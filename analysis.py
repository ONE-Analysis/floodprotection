"""
Analysis script for Flood Infrastructure Analysis
Handles flood exposure analysis and protection scenarios
"""
import os
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
from rasterio.features import rasterize, shapes
import shapely.geometry as sg
from shapely.geometry import Point, LineString, Polygon
from scipy.ndimage import binary_dilation
from config import DATASET_INFO, PROJECT_CRS, OVERWRITE, WEB_CRS
import sys
import warnings
import math
import shutil
from tqdm import tqdm

def extract_nsi_data(buildings_gdf, nsi_gdf):
    """
    Extract NSI data and join to buildings
    
    Parameters:
    buildings_gdf (GeoDataFrame): GeoDataFrame of buildings
    nsi_gdf (GeoDataFrame): GeoDataFrame of NSI data
    
    Returns:
    GeoDataFrame: Buildings with NSI data joined
    """
    print("Extracting NSI data and joining to buildings...")
    
    # Debug information
    print(f"Buildings GDF: {len(buildings_gdf)} features with CRS: {buildings_gdf.crs}")
    print(f"NSI GDF: {len(nsi_gdf)} features with CRS: {nsi_gdf.crs}")
    
    # Ensure CRS match
    if buildings_gdf.crs != nsi_gdf.crs:
        print(f"CRS mismatch: Buildings {buildings_gdf.crs}, NSI {nsi_gdf.crs}")
        print("Reprojecting NSI to match buildings CRS")
        nsi_gdf = nsi_gdf.to_crs(buildings_gdf.crs)
    
    # Print the bounding boxes of each dataset to verify overlap
    buildings_bounds = buildings_gdf.total_bounds
    nsi_bounds = nsi_gdf.total_bounds
    print(f"Buildings bounds: {buildings_bounds}")
    print(f"NSI bounds: {nsi_bounds}")
    
    # Check spatial bounds overlap
    bounds_overlap = (
        buildings_bounds[0] < nsi_bounds[2] and buildings_bounds[2] > nsi_bounds[0] and
        buildings_bounds[1] < nsi_bounds[3] and buildings_bounds[3] > nsi_bounds[1]
    )
    
    if not bounds_overlap:
        print("WARNING: Buildings and NSI datasets do not appear to overlap spatially!")
        print(f"Buildings bounds: {buildings_bounds}")
        print(f"NSI bounds: {nsi_bounds}")
    
    # Fields to extract from NSI
    nsi_fields = ['sqft', 'found_type', 'found_ht', 'med_yr_blt', 
                  'val_struct', 'val_cont', 'val_vehic', 'ground_elv']
    
    # Ensure all required fields exist in NSI data
    for field in nsi_fields:
        if field not in nsi_gdf.columns:
            print(f"Warning: Field '{field}' not found in NSI data. Creating empty column.")
            nsi_gdf[field] = None
    
    # Create centroids for spatial join
    buildings_centroids = buildings_gdf.copy()
    buildings_centroids['centroid'] = buildings_centroids.geometry.centroid
    buildings_centroids.set_geometry('centroid', inplace=True)
    
    # Spatial join NSI to buildings based on centroids
    joined = gpd.sjoin(buildings_centroids, nsi_gdf, how='left', predicate='within')
    
    # If no NSI points are within buildings, just leave as NULL (don't use nearest neighbor)
    if joined[nsi_fields].isna().all().all():
        print(f"No NSI points found within {len(buildings_gdf)} buildings, leaving fields as NULL values")
    
    # Reset geometry to original building polygons
    joined.set_geometry('geometry', inplace=True)
    
    # Get only the NSI fields and original building fields
    buildings_nsi = buildings_gdf.copy()
    for field in nsi_fields:
        buildings_nsi[field] = joined[field]
    
    return buildings_nsi

def analyze_flood_exposure(buildings_gdf, fema_flood_path):
    """
    Analyze flood exposure for buildings using a raster-based approach
    
    Parameters:
    buildings_gdf (GeoDataFrame): Buildings data
    fema_flood_path (str): Path to FEMA flood raster
    
    Returns:
    GeoDataFrame: Buildings with flood exposure
    """
    print("Analyzing flood exposure using raster-based approach...")
    
    # Create a copy of the buildings GeoDataFrame
    buildings_flood = buildings_gdf.copy()
    
    # Add flood exposure columns if they don't exist
    if 'flood_02pct' not in buildings_flood.columns:
        buildings_flood['flood_02pct'] = False
    if 'flood_01pct' not in buildings_flood.columns:
        buildings_flood['flood_01pct'] = False
    
    # Open the FEMA flood raster
    with rasterio.open(fema_flood_path) as src:
        # Print info about the raster
        print(f"FEMA raster shape: {src.shape}")
        fema_data = src.read(1)
        unique_values = np.unique(fema_data)
        print(f"FEMA raster unique values: {unique_values}")
        print(f"FEMA raster CRS: {src.crs}")
        
        # Check if buildings CRS matches raster CRS
        if buildings_flood.crs != src.crs:
            print(f"Converting buildings from {buildings_flood.crs} to {src.crs}")
            buildings_flood = buildings_flood.to_crs(src.crs)
        
        # Create masks for 1% and 0.2% flood zones
        # In FEMA raster: 1 = 0.2% annual chance, 2 = 1% annual chance
        mask_1pct = fema_data == 2
        mask_02pct = (fema_data == 1) | (fema_data == 2)
        
        # Process buildings in batches to prevent memory issues
        batch_size = 5000  # Increased from 1000 to 5000
        total_buildings = len(buildings_flood)
        num_batches = (total_buildings + batch_size - 1) // batch_size  # Ceiling division
        
        print(f"Processing {total_buildings} buildings in {num_batches} batches of {batch_size}")
        
        # Using tqdm for batch progress
        for batch_idx in tqdm(range(num_batches), desc="Processing building batches"):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, total_buildings)
            
            batch_buildings = buildings_flood.iloc[start_idx:end_idx]
            
            try:
                # Convert batch of building polygons to a raster
                # Create a dictionary to map building index to rasterized ID
                building_id_map = {}
                shapes_to_rasterize = []
                
                for i, (idx, building) in enumerate(batch_buildings.iterrows()):
                    # Use i+1 as raster value to avoid 0 (which is nodata)
                    raster_id = i + 1
                    building_id_map[raster_id] = idx
                    shapes_to_rasterize.append((building.geometry, raster_id))
                
                # Rasterize buildings
                buildings_raster = rasterize(
                    shapes_to_rasterize,
                    out_shape=src.shape,
                    transform=src.transform,
                    fill=0,
                    dtype=np.uint32,  # Use uint32 to handle large number of buildings
                    all_touched=True  # Consider all pixels touched by polygon
                )
                
                # Find buildings that overlap with each flood zone
                # For 1% flood zone
                buildings_in_1pct = np.unique(buildings_raster[mask_1pct])
                # Remove 0 value (background)
                buildings_in_1pct = buildings_in_1pct[buildings_in_1pct > 0]
                
                # For 0.2% flood zone
                buildings_in_02pct = np.unique(buildings_raster[mask_02pct])
                # Remove 0 value (background)
                buildings_in_02pct = buildings_in_02pct[buildings_in_02pct > 0]
                
                # Mark buildings as being in flood zones
                for raster_id in buildings_in_1pct:
                    idx = building_id_map[raster_id]
                    buildings_flood.at[idx, 'flood_01pct'] = True
                    buildings_flood.at[idx, 'flood_02pct'] = True  # 1% is also in 0.2%
                
                for raster_id in buildings_in_02pct:
                    idx = building_id_map[raster_id]
                    buildings_flood.at[idx, 'flood_02pct'] = True
                
            except Exception as e:
                print(f"Error processing batch {batch_idx+1}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        # Calculate total buildings in each flood zone
        total_02pct = buildings_flood['flood_02pct'].sum()
        total_01pct = buildings_flood['flood_01pct'].sum()
        print(f"Found {total_02pct} buildings in 0.2% flood zone and {total_01pct} in 1% flood zone")
    
    return buildings_flood

def create_trimmed_flood_raster(fema_flood_path, alignment_gdf, water_gdf, output_path):
    """
    Create a trimmed version of the FEMA flood raster based on the alignment
    Using a direct raster-based approach to determine which side of the alignment is protected
    
    Parameters:
    fema_flood_path (str): Path to FEMA flood raster
    alignment_gdf (GeoDataFrame): Alignment polyline
    water_gdf (GeoDataFrame): Water polygon
    output_path (str): Path to save trimmed raster
    
    Returns:
    str: Path to trimmed flood raster
    """
    print("Creating trimmed flood raster using raster-based approach...")
    
    # Create output directory
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if output already exists
    if os.path.exists(output_path) and not OVERWRITE["buildings_scenarios"]:
        print(f"Trimmed flood raster already exists at {output_path}")
        try:
            with rasterio.open(output_path) as src:
                data = src.read(1)
                if np.all(data == 0) or np.all(np.isnan(data)):
                    print("Existing trimmed raster appears to be empty, will recreate it")
                else:
                    return output_path
        except Exception as e:
            print(f"Error checking existing trimmed raster: {e}")
    
    # Check for valid inputs
    if alignment_gdf is None or alignment_gdf.empty:
        print("Error: Alignment GeoDataFrame is empty or None!")
        shutil.copy2(fema_flood_path, output_path)
        return output_path
        
    if water_gdf is None or water_gdf.empty:
        print("Error: Water GeoDataFrame is empty or None!")
        shutil.copy2(fema_flood_path, output_path)
        return output_path
    
    # Verify CRS
    if alignment_gdf.crs is None:
        print(f"ERROR: Alignment has no CRS. Setting to {PROJECT_CRS}")
        alignment_gdf.crs = PROJECT_CRS
    elif alignment_gdf.crs != PROJECT_CRS:
        print(f"Converting alignment from {alignment_gdf.crs} to {PROJECT_CRS}")
        alignment_gdf = alignment_gdf.to_crs(PROJECT_CRS)
        
    if water_gdf.crs is None:
        print(f"ERROR: Water has no CRS. Setting to {PROJECT_CRS}")
        water_gdf.crs = PROJECT_CRS
    elif water_gdf.crs != PROJECT_CRS:
        print(f"Converting water from {water_gdf.crs} to {PROJECT_CRS}")
        water_gdf = water_gdf.to_crs(PROJECT_CRS)
    
    try:
        # Get the alignment line
        alignment_line = alignment_gdf.geometry.iloc[0]
        
        # Save the alignment for reference
        alignment_path = os.path.join(output_dir, "alignment.geojson")
        alignment_gdf.to_file(alignment_path, driver="GeoJSON")
        print(f"Alignment saved to {alignment_path}")
        
        # Get the water geometry
        water_geom = water_gdf.geometry.unary_union
        
        # STEP 1: Split the water geometry with the alignment if they intersect
        if alignment_line.intersects(water_geom):
            print("Alignment intersects water - splitting water geometry")
            
            # Create a small buffer around the alignment to ensure clean split
            alignment_buffer = alignment_line.buffer(0.1)
            
            # Split the water geometry
            split_water = water_geom.difference(alignment_buffer)
            
            # Get the resulting geometries
            if split_water.geom_type == 'MultiPolygon':
                water_parts = list(split_water.geoms)
            elif split_water.geom_type == 'Polygon':
                water_parts = [split_water]
            else:
                print(f"Unexpected geometry type after split: {split_water.geom_type}")
                water_parts = [water_geom]  # Fall back to original water
            
            # Sort by area and keep the largest part
            if len(water_parts) > 1:
                water_parts.sort(key=lambda g: g.area, reverse=True)
                kept_water = water_parts[0]
                print(f"Split water into {len(water_parts)} parts. Kept largest with area {kept_water.area:.1f}")
            else:
                kept_water = water_parts[0]
                print("Water wasn't properly split, using original water geometry")
        else:
            print("Alignment doesn't intersect water - using full water geometry")
            kept_water = water_geom
        
        # Save kept water geometry for reference
        kept_water_gdf = gpd.GeoDataFrame(geometry=[kept_water], crs=PROJECT_CRS)
        kept_water_path = os.path.join(output_dir, "kept_water.geojson")
        kept_water_gdf.to_file(kept_water_path, driver="GeoJSON")
        
        # Get water centroid for reference
        water_centroid = kept_water.centroid
        
        # Open the FEMA flood raster
        with rasterio.open(fema_flood_path) as src:
            fema_data = src.read(1)
            transform = src.transform
            
            print(f"Creating side determination mask with shape {src.shape}")
            
            # STEP 2: Create a binary mask of the alignment with a buffer
            # This will help us visualize and debug the alignment
            alignment_raster = rasterize(
                [(alignment_line.buffer(10), 1)],  # Use a 10-unit buffer for visibility
                out_shape=src.shape,
                transform=transform,
                fill=0,
                dtype=np.uint8
            )
            
            # Save the alignment raster for debugging
            alignment_raster_path = os.path.join(output_dir, "alignment_raster.tif")
            with rasterio.open(alignment_raster_path, 'w', 
                              driver='GTiff', 
                              height=src.height, 
                              width=src.width, 
                              count=1, 
                              dtype=np.uint8, 
                              crs=src.crs, 
                              transform=transform) as dst:
                dst.write(alignment_raster, 1)
            print(f"Alignment raster saved to {alignment_raster_path}")
            
            # STEP 3: Determine which side of the alignment the water is on
            # Process EVERY pixel directly for maximum accuracy
            
            print("Creating full-resolution water side mask by checking every pixel...")
            
            # Create a full-resolution mask at the start
            full_water_side_mask = np.zeros(src.shape, dtype=bool)
            has_protected_side = False
            
            # Function to determine if a point is on the same side as water
            def is_on_water_side(point, alignment, water_centroid):
                """
                Check if a point is on the same side of the alignment as the water.
                
                Uses ray casting - draw a line from the point to the water and 
                count how many times it crosses the alignment. If even (or 0), 
                it's on the same side. If odd, it's on the opposite side.
                """
                ray = LineString([point, water_centroid])
                intersections = ray.intersection(alignment)
                
                if intersections.is_empty:
                    return True  # No crossings, same side
                elif intersections.geom_type == 'Point':
                    return False  # One crossing, opposite side
                elif intersections.geom_type == 'MultiPoint':
                    # Count number of intersection points
                    num_points = len(list(intersections.geoms))
                    return num_points % 2 == 0  # Even = same side, Odd = opposite side
                else:
                    # Unexpected result - default to water side to be safe
                    return True
            
            # Only process where there are actual flood pixels to improve performance
            flood_pixels = fema_data > 0
            flood_pixel_count = np.sum(flood_pixels)
            processed_count = 0
            
            # Get the total number of rows and cols
            rows, cols = src.shape
            
            # Process each pixel where there's flood data
            print(f"Processing {flood_pixel_count} flood pixels...")
            
            # Create a progress bar
            pbar = tqdm(total=flood_pixel_count, desc="Processing pixels")
            
            # Process in row-major order
            for row in range(rows):
                for col in range(cols):
                    # Only process flood pixels to save time
                    if fema_data[row, col] > 0:
                        # Skip points in the alignment buffer
                        if alignment_raster[row, col] > 0:
                            full_water_side_mask[row, col] = False
                        else:
                            # Convert to real-world coordinates
                            x, y = rasterio.transform.xy(transform, row, col)
                            point = Point(x, y)
                            
                            # Check if on water side
                            try:
                                on_water_side = is_on_water_side(point, alignment_line, water_centroid)
                                full_water_side_mask[row, col] = on_water_side
                                
                                # Track if we found protected areas
                                if not on_water_side:
                                    has_protected_side = True
                            except Exception as e:
                                # Silently handle errors and default to water side
                                full_water_side_mask[row, col] = True
                        
                        # Update progress
                        processed_count += 1
                        if processed_count % 1000 == 0:  # Update progress every 1000 pixels
                            pbar.update(1000)
                            
            # Close progress bar
            pbar.update(flood_pixel_count - processed_count)  # Update remaining count
            pbar.close()
            
            # Check if we found any protected areas
            if not has_protected_side:
                print("Warning: No protected side detected! The alignment may not be forming a barrier.")
                print("Using original flood data as there is no clear protected side.")
                shutil.copy2(fema_flood_path, output_path)
                return output_path
            
            # We now have a full-resolution mask - no need for STEP 4 interpolation
            
            # Save the water side mask for debugging
            water_side_mask_path = os.path.join(output_dir, "water_side_mask.tif")
            with rasterio.open(water_side_mask_path, 'w', 
                              driver='GTiff', 
                              height=src.height, 
                              width=src.width, 
                              count=1, 
                              dtype=np.uint8, 
                              crs=src.crs, 
                              transform=transform) as dst:
                dst.write(full_water_side_mask.astype(np.uint8), 1)
            print(f"Water side mask saved to {water_side_mask_path}")
            
            # STEP 5: Apply the water side mask to the flood raster
            print("Applying water side mask to flood raster...")
            
            # Create the output raster - only keep flood pixels on the water side
            trimmed_data = np.where(full_water_side_mask, fema_data, 0)
            
            # Check if result is empty (no flood pixels on water side)
            if np.all(trimmed_data == 0):
                print("Warning: No flood pixels on water side! Using original flood raster.")
                shutil.copy2(fema_flood_path, output_path)
                return output_path
            
            # Save the trimmed flood raster
            with rasterio.open(output_path, 'w', **src.meta) as dst:
                dst.write(trimmed_data, 1)
            
            print(f"Trimmed flood raster saved to {output_path}")
    
    except Exception as e:
        print(f"Error in create_trimmed_flood_raster: {e}")
        import traceback
        traceback.print_exc()
        
        # Fall back to original flood raster
        print("Using original flood raster due to error")
        shutil.copy2(fema_flood_path, output_path)
    
    return output_path

def analyze_protection_scenario(buildings_flood, trimmed_fema_path, flood_fraction_threshold=0.5):
    """
    Analyze protection scenario using the trimmed flood map by checking the fraction 
    of each building's area that falls in the trimmed flood zone.
    
    A building is considered to remain flooded in scenario1 (with protection) if it still has more
    than the threshold fraction of its area in the flood zone after the protective alignment.
    
    Parameters:
    buildings_flood (GeoDataFrame): Buildings with original flood exposure flags.
    trimmed_fema_path (str): Path to the trimmed FEMA flood raster.
    flood_fraction_threshold (float): Fraction threshold above which a building is considered still flooded.
    
    Returns:
    GeoDataFrame: Buildings with additional protection scenario fields.
    """
    import rasterio
    from rasterio.features import rasterize
    import numpy as np

    print("Analyzing protection scenario (per-building flood fraction method)...")
    
    # Make a copy of the buildings data and initialize scenario flags.
    buildings_scenarios = buildings_flood.copy()
    buildings_scenarios['flood_01pct_scenario1'] = False
    buildings_scenarios['flood_02pct_scenario1'] = False

    # Open the trimmed FEMA flood raster.
    with rasterio.open(trimmed_fema_path) as src:
        transform = src.transform
        trimmed_fema_data = src.read(1)
        # Define flood masks (FEMA convention: 1 = 0.2%, 2 = 1%).
        mask_1pct = (trimmed_fema_data == 2)
        mask_02pct = ((trimmed_fema_data == 1) | (trimmed_fema_data == 2))
        
        # Ensure building geometries are in the same CRS as the raster.
        if buildings_scenarios.crs != src.crs:
            print(f"Converting buildings from {buildings_scenarios.crs} to {src.crs}")
            buildings_scenarios = buildings_scenarios.to_crs(src.crs)
        
        # Evaluate each building individually.
        for idx, building in buildings_scenarios.iterrows():
            # Rasterize the building footprint onto the same grid.
            building_mask = rasterize(
                [(building.geometry, 1)],
                out_shape=src.shape,
                transform=transform,
                fill=0,
                dtype=np.uint8,
                all_touched=True
            ) == 1  # ensure boolean mask
            
            total_pixels = np.sum(building_mask)
            if total_pixels == 0:
                continue  # skip if no pixels covered
            
            # Count how many pixels fall into the flood zones.
            flooded_pixels_1pct = np.sum(np.logical_and(building_mask, mask_1pct))
            flooded_pixels_02pct = np.sum(np.logical_and(building_mask, mask_02pct))
            
            # Compute the fraction of the building that is flooded.
            fraction_1pct = flooded_pixels_1pct / total_pixels
            fraction_02pct = flooded_pixels_02pct / total_pixels
            
            # Mark building as still flooded if it has more than threshold fraction flooded
            if building['flood_01pct'] and fraction_1pct >= flood_fraction_threshold:
                buildings_scenarios.at[idx, 'flood_01pct_scenario1'] = True
            if building['flood_02pct'] and fraction_02pct >= flood_fraction_threshold:
                buildings_scenarios.at[idx, 'flood_02pct_scenario1'] = True

    # Report results.
    total_1pct = buildings_scenarios['flood_01pct'].sum()
    total_02pct = buildings_scenarios['flood_02pct'].sum()
    still_flooded_1pct = buildings_scenarios['flood_01pct_scenario1'].sum()
    still_flooded_02pct = buildings_scenarios['flood_02pct_scenario1'].sum()
    
    print(f"{still_flooded_1pct} of {total_1pct} buildings remain in 1% flood zone ({still_flooded_1pct/max(1,total_1pct)*100:.1f}%)")
    print(f"{still_flooded_02pct} of {total_02pct} buildings remain in 0.2% flood zone ({still_flooded_02pct/max(1,total_02pct)*100:.1f}%)")
    
    return buildings_scenarios

def analyze_buildings_in_floodplain(buildings_gdf, flood_raster_path, output_fields=None, base_fields=None):
    """
    Analyze buildings to determine which ones are in flood zones using a raster-based approach
    
    Parameters:
    buildings_gdf (GeoDataFrame): Buildings to analyze
    flood_raster_path (str): Path to flood raster (original or trimmed)
    output_fields (tuple): Field names to store results (0.2% field, 1% field)
    base_fields (tuple): Optional base fields to copy values from (for comparison)
    
    Returns:
    GeoDataFrame: Buildings with flood zone analysis results
    """
    print(f"Analyzing buildings using flood raster: {os.path.basename(flood_raster_path)}")
    
    # Define output field names if not provided
    if output_fields is None:
        output_fields = ('flood_02pct', 'flood_01pct')
    
    # Create a copy of the buildings GeoDataFrame
    result_gdf = buildings_gdf.copy()
    
    # Initialize output columns to False
    result_gdf[output_fields[0]] = False  # 0.2% annual chance field
    result_gdf[output_fields[1]] = False  # 1% annual chance field
    
    # Copy values from base fields if provided (for comparison scenarios)
    if base_fields is not None:
        result_gdf[output_fields[0]] = result_gdf[base_fields[0]]
        result_gdf[output_fields[1]] = result_gdf[base_fields[1]]
    
    # Open the flood raster
    with rasterio.open(flood_raster_path) as src:
        # Print info about the raster
        print(f"Flood raster shape: {src.shape}")
        flood_data = src.read(1)
        unique_values = np.unique(flood_data)
        print(f"Flood raster unique values: {unique_values}")
        
        # Check if buildings CRS matches raster CRS
        if result_gdf.crs != src.crs:
            print(f"Converting buildings from {result_gdf.crs} to {src.crs}")
            result_gdf = result_gdf.to_crs(src.crs)
        
        # Create masks for 1% and 0.2% flood zones
        # In FEMA raster: 1 = 0.2% annual chance, 2 = 1% annual chance
        mask_1pct = flood_data == 2
        mask_02pct = (flood_data == 1) | (flood_data == 2)
        
        # Count pixels in each mask to verify they exist
        print(f"Pixels in 1% flood zone: {np.sum(mask_1pct)}")
        print(f"Pixels in 0.2% flood zone: {np.sum(mask_02pct)}")
        
        # Process buildings in batches to prevent memory issues
        batch_size = 5000
        total_buildings = len(result_gdf)
        num_batches = (total_buildings + batch_size - 1) // batch_size  # Ceiling division
        
        print(f"Processing {total_buildings} buildings in {num_batches} batches of {batch_size}")
        
        # Count for tracking
        buildings_in_1pct_count = 0
        buildings_in_02pct_count = 0
        
        # Using tqdm for batch progress
        for batch_idx in tqdm(range(num_batches), desc="Processing building batches"):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, total_buildings)
            
            batch_buildings = result_gdf.iloc[start_idx:end_idx]
            
            try:
                # Create a dictionary to map building index to rasterized ID
                building_id_map = {}
                shapes_to_rasterize = []
                
                for i, (idx, building) in enumerate(batch_buildings.iterrows()):
                    # Use i+1 as raster value to avoid 0 (which is nodata)
                    raster_id = i + 1
                    building_id_map[raster_id] = idx
                    shapes_to_rasterize.append((building.geometry, raster_id))
                
                # Rasterize buildings
                buildings_raster = rasterize(
                    shapes_to_rasterize,
                    out_shape=src.shape,
                    transform=src.transform,
                    fill=0,
                    dtype=np.uint32,  # Use uint32 to handle large number of buildings
                    all_touched=True  # Consider all pixels touched by polygon
                )
                
                # Find buildings that intersect with each flood zone
                buildings_in_1pct = np.unique(buildings_raster[mask_1pct])
                buildings_in_1pct = buildings_in_1pct[buildings_in_1pct > 0]  # Remove 0 (background)
                
                buildings_in_02pct = np.unique(buildings_raster[mask_02pct])
                buildings_in_02pct = buildings_in_02pct[buildings_in_02pct > 0]  # Remove 0 (background)
                
                # Mark buildings that are in flood zones
                for raster_id in buildings_in_1pct:
                    idx = building_id_map[raster_id]
                    # For protection scenario, we set to False because building is still flooded
                    result_gdf.at[idx, output_fields[1]] = True if base_fields is None else False
                    buildings_in_1pct_count += 1
                
                for raster_id in buildings_in_02pct:
                    idx = building_id_map[raster_id]
                    # For protection scenario, we set to False because building is still flooded
                    result_gdf.at[idx, output_fields[0]] = True if base_fields is None else False
                    buildings_in_02pct_count += 1
                
            except Exception as e:
                print(f"Error processing batch {batch_idx+1}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        # Print summary of results
        print(f"Identified {buildings_in_1pct_count} buildings in 1% flood zone")
        print(f"Identified {buildings_in_02pct_count} buildings in 0.2% flood zone")
    
    return result_gdf

def calculate_flood_damage(buildings_df, nsi_data):
    """
    Calculate estimated flood damages using NSI data
    
    Parameters:
    buildings_df (GeoDataFrame): Buildings with flood exposure analysis
    nsi_data (GeoDataFrame): NSI data with damage estimates
    
    Returns:
    GeoDataFrame: Buildings with damage estimates added
    """
    print("Calculating flood damage estimates...")
    
    # Create a copy to avoid modifying the input
    buildings_damage = buildings_df.copy()
    
    # Add damage estimate columns
    buildings_damage['total_damage_01pct'] = 0.0
    buildings_damage['total_damage_02pct'] = 0.0
    
    # Check if NSI data is available
    if nsi_data is None or len(nsi_data) == 0:
        print("No NSI data available for damage calculations")
        return buildings_damage
    
    try:
        # Ensure the projections match
        if buildings_damage.crs != nsi_data.crs:
            nsi_data = nsi_data.to_crs(buildings_damage.crs)
        
        # Spatial join to connect buildings with NSI damage curves
        print("Joining buildings with NSI damage data...")
        buildings_nsi = gpd.sjoin_nearest(
            buildings_damage, 
            nsi_data,
            how='left',
            max_distance=50  # 50 units (usually meters) max distance
        )
        
        # Handle buildings without a match (use median damage values by occupancy type)
        missing_nsi = buildings_nsi['nsi_id'].isna()
        if missing_nsi.any():
            print(f"Note: {missing_nsi.sum()} buildings ({missing_nsi.sum()/len(buildings_nsi)*100:.1f}%) without NSI match")
            # Use building type and median damage by occupancy as fallback
            # This would require additional processing
        
        # Calculate damages for buildings in flood zones
        print("Calculating damages for buildings in flood zones...")
        
        # For 1% annual chance flood (100-year)
        mask_1pct = buildings_nsi['flood_01pct']
        buildings_damage.loc[mask_1pct, 'total_damage_01pct'] = buildings_nsi.loc[mask_1pct, 'damage_1pct']
        
        # For 0.2% annual chance flood (500-year)
        mask_02pct = buildings_nsi['flood_02pct']
        buildings_damage.loc[mask_02pct, 'total_damage_02pct'] = buildings_nsi.loc[mask_02pct, 'damage_02pct']
        
        # Calculate total potential damages
        total_damage_1pct = buildings_damage['total_damage_01pct'].sum()
        total_damage_02pct = buildings_damage['total_damage_02pct'].sum()
        
        print(f"Total potential damage (1% flood): ${total_damage_1pct:,.2f}")
        print(f"Total potential damage (0.2% flood): ${total_damage_02pct:,.2f}")
        
    except Exception as e:
        print(f"Error in damage calculation: {e}")
        import traceback
        traceback.print_exc()
    
    return buildings_damage

def run_analysis(preprocessed_data):
    """
    Run the analysis for both flood exposure scenarios:
      1. Original flood exposure using the full FEMA flood raster.
      2. Protection scenario exposure using the trimmed flood raster.
      
    Returns:
        dict: Dictionary with both analysis results.
    """
    print("Starting analysis...")

    # Unpack preprocessed data
    buildings_gdf = preprocessed_data["buildings"]
    fema_flood_crop_path = preprocessed_data["fema_flood_crop"]
    alignment_gdf = preprocessed_data["alignment"]
    water_gdf = preprocessed_data["water"]
    nsi_gdf = preprocessed_data["nsi"]

    print(f"Analysis using {len(buildings_gdf)} buildings")
    if buildings_gdf is None or buildings_gdf.empty:
        print("Error: No buildings data")
        sys.exit(1)
    if fema_flood_crop_path is None or not os.path.exists(fema_flood_crop_path):
        print("Error: No FEMA flood crop data")
        sys.exit(1)

    # --- Scenario 1: Original Flood Exposure ---
    # Optionally join NSI data if available
    if nsi_gdf is not None and not nsi_gdf.empty:
        buildings_gdf = extract_nsi_data(buildings_gdf, nsi_gdf)
    
    # Analyze flood exposure using the original FEMA flood raster.
    print("Running original flood exposure analysis...")
    buildings_flood = analyze_flood_exposure(buildings_gdf, fema_flood_crop_path)
    
    # Save the original flood exposure results.
    buildings_flood_path = DATASET_INFO["Output"]["Buildings_flood"]["path"]
    buildings_flood.to_file(buildings_flood_path, driver='GeoJSON')
    print(f"Buildings flood exposure saved to {buildings_flood_path}")

    # --- Create the Trimmed Flood Raster for Protection Scenario ---
    if alignment_gdf is not None and not alignment_gdf.empty:
        fema_flood_trim_path = DATASET_INFO["Output"]["FEMA_Flood_crop_trim"]["path"]
        trimmed_fema_path = create_trimmed_flood_raster(
            fema_flood_crop_path, 
            alignment_gdf, 
            water_gdf,
            fema_flood_trim_path
        )
    else:
        print("Warning: No alignment data available, using original flood raster for protection scenario.")
        trimmed_fema_path = fema_flood_crop_path

    # --- Scenario 2: Protection Scenario ---
    # This function compares the original exposure with the trimmed flood raster.
    print("Running protection scenario analysis...")
    buildings_protected = analyze_protection_scenario(buildings_flood, trimmed_fema_path)
    
    # Save the protection scenario results.
    buildings_scenarios_path = DATASET_INFO["Output"]["Buildings_scenarios"]["path"]
    buildings_protected.to_file(buildings_scenarios_path, driver='GeoJSON')
    print(f"Buildings protection scenario saved to {buildings_scenarios_path}")

    # Optionally, you can merge the results so that each building has both sets of attributes.
    # For example:
    # merged = buildings_flood.copy()
    # merged["flood_01pct_protected"] = buildings_protected["flood_01pct_protected"]
    # merged["flood_02pct_protected"] = buildings_protected["flood_02pct_protected"]
    # merged.to_file("path_to_merged_output.geojson", driver='GeoJSON')
    
    # Continue with web map generation, etc.
    return {
        "buildings_flood": buildings_flood,
        "buildings_protected": buildings_protected,
        "alignment": alignment_gdf
    }

def prepare_data_for_webmap(analysis_results):
    """
    Prepare analysis results for webmap visualization
    
    Parameters:
    analysis_results (dict): Analysis results
    
    Returns:
    dict: Data prepared for webmap
    """
    import os
    import numpy as np
    from config import DATASET_INFO, WEB_CRS
    import geopandas as gpd
    
    webmap_data = {}
    
    # Process FEMA flood raster for web display
    original_path = DATASET_INFO["Output"]["FEMA_Flood_crop"]["path"]
    web_path = DATASET_INFO["Output"]["FEMA_Flood_crop_web"]["path"]
    if os.path.exists(str(original_path)):
        from webmap import process_raster_for_web
        webmap_data["fema_web_path"], webmap_data["fema_bounds"] = process_raster_for_web(
            original_path, web_path
        )
    else:
        print(f"Warning: Original FEMA flood raster not found at {original_path}")
        webmap_data["fema_web_path"] = None
        webmap_data["fema_bounds"] = [0, 0, 0, 0]
    
    # Process trimmed FEMA flood raster
    trim_path = DATASET_INFO["Output"]["FEMA_Flood_crop_trim"]["path"]
    trim_web_path = DATASET_INFO["Output"]["FEMA_Flood_crop_trim_web"]["path"]
    if os.path.exists(str(trim_path)):
        from webmap import process_raster_for_web
        webmap_data["fema_trim_web_path"], webmap_data["fema_trim_bounds"] = process_raster_for_web(
            trim_path, trim_web_path
        )
    else:
        print(f"Warning: Trimmed FEMA flood raster not found at {trim_path}")
        webmap_data["fema_trim_web_path"] = None
        webmap_data["fema_trim_bounds"] = [0, 0, 0, 0]
    
    # Include alignment data
    if "alignment" in analysis_results and analysis_results["alignment"] is not None:
        alignment = analysis_results["alignment"]
        if alignment.crs != WEB_CRS:
            alignment = alignment.to_crs(WEB_CRS)
        webmap_data["alignment"] = alignment
    else:
        webmap_data["alignment"] = None
    
    # Include buildings data
    if "buildings_protected" in analysis_results:
        buildings = analysis_results["buildings_protected"]
        if buildings.crs != WEB_CRS:
            buildings = buildings.to_crs(WEB_CRS)
        webmap_data["buildings"] = buildings
    else:
        webmap_data["buildings"] = None
    
    return webmap_data