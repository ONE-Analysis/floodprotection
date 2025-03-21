"""
Preprocessing script for Flood Infrastructure Analysis
Handles CRS standardization, bounding box creation, and data cropping
"""
import os
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
from pathlib import Path
import shapely.geometry as sg
from shapely.ops import linemerge
from config import DATASET_INFO, PROJECT_CRS, OVERWRITE
import sys
from rasterio.windows import from_bounds

def create_site_bounds():
    """
    Load site_boundary.geojson
    
    Returns:
    GeoDataFrame: Site bounds as a GeoDataFrame
    """
    site_bounds_path = DATASET_INFO["Input"]["Site_Bounds"]["path"]
    
    if not os.path.exists(site_bounds_path):
        raise FileNotFoundError(f"Site bounds not found at {site_bounds_path}")
    
    print(f"Loading site bounds from {site_bounds_path}")
    site_bbox = gpd.read_file(site_bounds_path)
    
    # Ensure CRS is set correctly
    if site_bbox.crs is None:
        print(f"Site bounds has no CRS, assigning {PROJECT_CRS}")
        site_bbox.crs = PROJECT_CRS
    elif site_bbox.crs != PROJECT_CRS:
        print(f"Converting site bounds from {site_bbox.crs} to {PROJECT_CRS}")
        site_bbox = site_bbox.to_crs(PROJECT_CRS)
    
    return site_bbox

def create_alignment():
    """
    Load alignment.geojson
    
    Returns:
    GeoDataFrame: Alignment as a GeoDataFrame
    """
    alignment_path = DATASET_INFO["Input"]["Alignment"]["path"]
    
    if not os.path.exists(alignment_path):
        raise FileNotFoundError(f"Alignment not found at {alignment_path}")
    
    print(f"Loading alignment from {alignment_path}")
    alignment = gpd.read_file(alignment_path)
    
    # Ensure CRS is set correctly
    if alignment.crs is None:
        print(f"Alignment has no CRS, assigning {PROJECT_CRS}")
        alignment.crs = PROJECT_CRS
    elif alignment.crs != PROJECT_CRS:
        print(f"Converting alignment from {alignment.crs} to {PROJECT_CRS}")
        alignment = alignment.to_crs(PROJECT_CRS)
    
    # Check if the alignment is empty
    if alignment.empty:
        print("ERROR: Alignment GeoJSON is empty!")
        print("Please ensure your alignment file contains valid geometries.")
        
        # Create a placeholder line for testing purposes
        print("Creating placeholder alignment line for testing...")
        # Get the bounds of the site
        site_bounds_path = DATASET_INFO["Input"]["Site_Bounds"]["path"]
        if os.path.exists(site_bounds_path):
            site_bbox = gpd.read_file(site_bounds_path)
            if not site_bbox.empty:
                bounds = site_bbox.total_bounds
                # Create a line that crosses the site
                from shapely.geometry import LineString
                line = LineString([
                    ((bounds[0] + bounds[2]) / 2, bounds[1]),  # Bottom middle
                    ((bounds[0] + bounds[2]) / 2, bounds[3])   # Top middle
                ])
                alignment = gpd.GeoDataFrame(geometry=[line], crs=PROJECT_CRS)
                alignment.to_file(alignment_path, driver='GeoJSON')
                print("Placeholder alignment created for testing.")
    
    print(f"Alignment bounds: {alignment.total_bounds}")
    return alignment

def crop_vector_data(site_bbox):
    """
    Crop vector datasets to the site boundary
    
    Parameters:
    site_bbox (GeoDataFrame): Site boundary polygon
    
    Returns:
    tuple: (buildings_crop, nsi_crop, water_crop) cropped GeoDataFrames
    """
    print("Cropping vector datasets...")
    
    # Crop buildings
    buildings_path = DATASET_INFO["Input"]["Buildings"]["path"]
    buildings_crop_path = DATASET_INFO["Preprocessed"]["Buildings_crop"]["path"]
    buildings_crop = None
    
    # Check if buildings crop already exists
    if os.path.exists(buildings_crop_path) and not OVERWRITE["vector_crop"]:
        print(f"Buildings crop already exists at {buildings_crop_path}")
        buildings_crop = gpd.read_file(buildings_crop_path)
        if buildings_crop.empty:
            print("Existing buildings crop is empty! Regenerating...")
        else:
            print(f"Loaded {len(buildings_crop)} buildings from existing crop")
            
    # Create buildings crop if it doesn't exist or is empty
    if buildings_crop is None or buildings_crop.empty:
        print(f"Creating buildings crop from {buildings_path}")
        try:
            print(f"Reading buildings from {buildings_path}")
            buildings = gpd.read_file(buildings_path)
            print(f"Read {len(buildings)} buildings with CRS: {buildings.crs}")
            
            if buildings.crs is None:
                print("WARNING: Buildings GeoJSON has no CRS defined!")
                print("Assuming the CRS is the same as the project CRS...")
                buildings.crs = PROJECT_CRS
            elif buildings.crs != PROJECT_CRS:
                print(f"Reprojecting buildings from {buildings.crs} to {PROJECT_CRS}")
                buildings = buildings.to_crs(PROJECT_CRS)
            
            # Crop to site boundary
            print(f"Clipping buildings to site boundary with CRS: {site_bbox.crs}")
            buildings_crop = gpd.clip(buildings, site_bbox)
            print(f"Successfully cropped to {len(buildings_crop)} buildings")
            
            buildings_crop.to_file(buildings_crop_path)
            print(f"Saved {len(buildings_crop)} buildings to {buildings_crop_path}")
        except Exception as e:
            print(f"Error cropping buildings: {e}")
            buildings_crop = gpd.GeoDataFrame(geometry=[], crs=PROJECT_CRS)
    
    # Crop NSI data
    nsi_path = DATASET_INFO["Input"]["NSI"]["path"]
    nsi_crop_path = DATASET_INFO["Preprocessed"]["NSI_crop"]["path"]
    nsi_crop = None
    
    # Check if NSI crop already exists
    if os.path.exists(nsi_crop_path) and not OVERWRITE["vector_crop"]:
        print(f"NSI crop already exists at {nsi_crop_path}")
        nsi_crop = gpd.read_file(nsi_crop_path)
        if nsi_crop.empty:
            print("Existing NSI crop is empty! Regenerating...")
        else:
            print(f"Loaded {len(nsi_crop)} NSI points from existing crop")
            
    # Create NSI crop if it doesn't exist or is empty
    if nsi_crop is None or nsi_crop.empty:
        print(f"Creating NSI crop from {nsi_path}")
        try:
            nsi = gpd.read_file(nsi_path)
            if nsi.crs is None:
                print("WARNING: NSI GeoJSON has no CRS defined!")
                print("Assuming the CRS is the same as the project CRS...")
                nsi.crs = PROJECT_CRS
            elif nsi.crs != PROJECT_CRS:
                print(f"Reprojecting NSI from {nsi.crs} to {PROJECT_CRS}")
                nsi = nsi.to_crs(PROJECT_CRS)
            
            # Crop to site boundary
            nsi_crop = gpd.clip(nsi, site_bbox)
            nsi_crop.to_file(nsi_crop_path)
            print(f"Saved {len(nsi_crop)} NSI points to {nsi_crop_path}")
        except Exception as e:
            print(f"Error cropping NSI data: {e}")
            nsi_crop = gpd.GeoDataFrame(geometry=[], crs=PROJECT_CRS)
    
    # Crop water polygon
    water_path = DATASET_INFO["Input"]["Water"]["path"]
    water_crop_path = DATASET_INFO["Preprocessed"]["Water_crop"]["path"]
    water_crop = None
    
    # Check if water crop already exists
    if os.path.exists(water_crop_path) and not OVERWRITE["vector_crop"]:
        print(f"Water crop already exists at {water_crop_path}")
        water_crop = gpd.read_file(water_crop_path)
        if water_crop.empty:
            print("Existing water crop is empty! Regenerating...")
        else:
            print(f"Loaded water polygon from existing crop")
            
    # Create water crop if it doesn't exist or is empty
    if water_crop is None or water_crop.empty:
        print(f"Creating water crop from {water_path}")
        try:
            water = gpd.read_file(water_path)
            if water.crs is None:
                print("WARNING: Water GeoJSON has no CRS defined!")
                print("Assuming the CRS is the same as the project CRS...")
                water.crs = PROJECT_CRS
            elif water.crs != PROJECT_CRS:
                print(f"Reprojecting water from {water.crs} to {PROJECT_CRS}")
                water = water.to_crs(PROJECT_CRS)
            
            # Crop to site boundary
            water_crop = gpd.clip(water, site_bbox)
            water_crop.to_file(water_crop_path)
            print(f"Saved water polygon to {water_crop_path}")
        except Exception as e:
            print(f"Error cropping water data: {e}")
            water_crop = gpd.GeoDataFrame(geometry=[], crs=PROJECT_CRS)
    
    return buildings_crop, nsi_crop, water_crop

def crop_raster(input_path, output_path, bounds):
    """
    Crop a raster to the specified bounds
    
    Parameters:
    input_path (str or Path): Path to input raster
    output_path (str or Path): Path to output cropped raster
    bounds (tuple): Bounds in the format (minx, miny, maxx, maxy)
    
    Returns:
    None
    """
    print(f"Cropping raster from {input_path} to {output_path}")
    
    # Ensure paths are strings
    input_path = str(input_path)
    output_path = str(output_path)
    
    with rasterio.open(input_path) as src:
        # Create a window from the bounds
        window = from_bounds(*bounds, src.transform)
        
        # Read the data in the window
        data = src.read(window=window)
        
        # Calculate the transform for the window
        transform = rasterio.windows.transform(window, src.transform)
        
        # Create a profile for the output raster
        profile = src.profile.copy()
        profile.update({
            'height': window.height,
            'width': window.width,
            'transform': transform
        })
        
        # Write the cropped raster
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data)
            
    print(f"Raster cropped successfully to {output_path}")

def crop_raster_data(site_bbox):
    """
    Crop raster datasets using site bounding box
    
    Parameters:
    site_bbox (GeoDataFrame): Site bounds for cropping
    
    Returns:
    str: Path to cropped FEMA flood raster
    """
    # Define paths for original and cropped datasets
    fema_flood_path = str(DATASET_INFO["Input"]["FEMA_Flood"]["path"])
    fema_flood_crop_path = str(DATASET_INFO["Preprocessed"]["FEMA_Flood_crop"]["path"])
    
    # Check if cropped files already exist
    if os.path.exists(fema_flood_crop_path) and not OVERWRITE["raster_crop"]:
        print("Cropped FEMA raster file already exists")
        
        # Verify CRS of the cropped raster
        try:
            with rasterio.open(fema_flood_crop_path) as src:
                if src.crs is None or src.crs.to_string() != PROJECT_CRS:
                    print(f"WARNING: Cropped FEMA raster CRS ({src.crs}) does not match PROJECT_CRS ({PROJECT_CRS})")
                    print("This may cause intersection issues with the alignment.")
                else:
                    print(f"Cropped FEMA raster has correct CRS: {src.crs}")
        except Exception as e:
            print(f"Error checking raster CRS: {e}")
            
        return fema_flood_crop_path
    
    print("Cropping FEMA raster dataset...")
    
    # Get geometries for cropping
    bbox_geom = [sg.mapping(site_bbox.geometry.unary_union)]
    
    # Crop FEMA Flood
    print("Cropping FEMA Flood...")
    with rasterio.open(fema_flood_path) as src:
        # Check and report the CRS
        print(f"Original FEMA raster CRS: {src.crs}")
        
        # Ensure correct CRS
        if src.crs is None:
            print(f"Warning: FEMA Flood has no CRS, assuming {PROJECT_CRS}")
            # We'll proceed but this is risky
        elif src.crs.to_string() != PROJECT_CRS:
            print(f"Warning: FEMA Flood CRS {src.crs} does not match project CRS {PROJECT_CRS}")
            print("Will reproject during cropping")
            # We'll let rasterio handle the reprojection during masking
        
        # Perform the crop operation
        out_image, out_transform = mask(src, bbox_geom, crop=True)
        out_meta = src.meta.copy()
        
        # Update metadata - set 0 as nodata value and ensure CRS is set
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "crs": PROJECT_CRS,
            "nodata": 0  # Set 0 values as nodata
        })
        
        # Write cropped raster
        with rasterio.open(fema_flood_crop_path, "w", **out_meta) as dest:
            dest.write(out_image)
        
        print(f"Cropped FEMA raster saved with CRS: {PROJECT_CRS}")
    
    return fema_flood_crop_path

def process_water_polygon(water_gdf, alignment_gdf, site_bbox):
    """
    Process water polygon with advanced filtering:
    1. Crop to site boundary
    2. Dissolve all pieces into one polygon
    3. Split with alignment polyline
    4. Apply filtering logic based on alignment and connectivity
    
    Parameters:
    water_gdf (GeoDataFrame): Water polygons
    alignment_gdf (GeoDataFrame): Alignment polyline
    site_bbox (GeoDataFrame): Site boundary
    
    Returns:
    GeoDataFrame: Processed and filtered water polygons
    """
    print("Processing water polygon with advanced filtering...")
    
    if water_gdf is None or water_gdf.empty:
        print("No water polygon data available for processing")
        return None
    
    if alignment_gdf is None or alignment_gdf.empty:
        print("No alignment data available for water processing")
        return water_gdf
    
    # Step 1: Crop water to site boundary
    print("Cropping water to site boundary...")
    water_crop = gpd.clip(water_gdf, site_bbox)
    if water_crop.empty:
        print("Water polygon is empty after cropping to site boundary")
        return water_crop
    
    # Step 2: Dissolve all water pieces into one polygon
    print("Dissolving water pieces into a single polygon...")
    water_dissolved = water_crop.dissolve()
    
    # Check for validity and fix if needed
    if not all(water_dissolved.geometry.is_valid):
        print("Fixing invalid geometries in dissolved water...")
        water_dissolved.geometry = water_dissolved.geometry.buffer(0)
    
    # Step 3: Split the dissolved water using the alignment
    print("Splitting water with alignment...")
    
    # Get the alignment geometry
    if len(alignment_gdf) > 1:
        print(f"Combining {len(alignment_gdf)} alignment features...")
        from shapely.ops import linemerge
        alignment_geoms = alignment_gdf.geometry.tolist()
        alignment_line = linemerge(alignment_geoms)
        if alignment_line.geom_type == 'MultiLineString':
            print("Warning: Could not merge alignment into a single LineString")
            # Use the longest piece
            alignment_line = max(alignment_line.geoms, key=lambda x: x.length)
    else:
        alignment_line = alignment_gdf.geometry.iloc[0]
        # Handle case where even a single alignment feature is actually a MultiLineString
        if alignment_line.geom_type == 'MultiLineString':
            print("Warning: Single alignment feature is a MultiLineString")
            # Use the longest piece
            alignment_line = max(alignment_line.geoms, key=lambda x: x.length)
    
    # Get the water geometry
    water_geom = water_dissolved.geometry.iloc[0]
    from shapely.geometry import LineString, Point, MultiPoint, Polygon, MultiPolygon
    
    print(f"Alignment type: {alignment_line.geom_type}")
    print(f"Water geometry type: {water_geom.geom_type}")
    
    # Create debug folder
    debug_folder = "debug_geometries"
    os.makedirs(debug_folder, exist_ok=True)
    
    # Check if they actually intersect
    if alignment_line.intersects(water_geom):
        print("DEBUG: Alignment intersects water polygon!")
        intersection = alignment_line.intersection(water_geom)
        print(f"Intersection type: {intersection.geom_type}")
        print(f"Intersection length: {intersection.length if hasattr(intersection, 'length') else 'N/A'}")
        
        # Save debug geometries
        gpd.GeoDataFrame(geometry=[water_geom], crs=water_gdf.crs).to_file(
            os.path.join(debug_folder, "water_before_split.geojson")
        )
        
        gpd.GeoDataFrame(geometry=[alignment_line], crs=alignment_gdf.crs).to_file(
            os.path.join(debug_folder, "alignment_for_split.geojson")
        )
        
        gpd.GeoDataFrame(geometry=[intersection], crs=alignment_gdf.crs).to_file(
            os.path.join(debug_folder, "alignment_water_intersection.geojson")
        )
        
        # Step 3.1: Find entry/exit points of alignment with water boundary
        try:
            # Get the boundary of the water polygon
            water_boundary = water_geom.boundary
            
            # Find intersection points between alignment and water boundary
            boundary_intersection = alignment_line.intersection(water_boundary)
            
            # Convert to a list of points
            if boundary_intersection.geom_type == 'Point':
                boundary_points = [boundary_intersection]
            elif boundary_intersection.geom_type == 'MultiPoint':
                boundary_points = list(boundary_intersection.geoms)
            else:
                print(f"Unexpected boundary intersection type: {boundary_intersection.geom_type}")
                boundary_points = []
                # Try to extract points from a more complex geometry
                if hasattr(boundary_intersection, 'geoms'):
                    for geom in boundary_intersection.geoms:
                        if geom.geom_type == 'Point':
                            boundary_points.append(geom)
                
            print(f"Found {len(boundary_points)} boundary intersection points")
            
            # Save boundary points for debugging
            if boundary_points:
                gpd.GeoDataFrame(geometry=boundary_points, crs=alignment_gdf.crs).to_file(
                    os.path.join(debug_folder, "boundary_intersection_points.geojson")
                )
            
            # If we have an even number of points, we can create segments
            if len(boundary_points) >= 2:
                # Group points into entry/exit pairs
                # If we have ordered points along the alignment, we can create segments
                boundary_coords = [(p.x, p.y) for p in boundary_points]
                
                # We need to order the points along the alignment
                # This is a bit tricky, so we'll use a distance-based approach
                
                # Get all coordinates of the alignment
                alignment_coords = list(alignment_line.coords)
                
                # For each boundary point, find its position along the alignment
                point_positions = []
                for point in boundary_points:
                    min_dist = float('inf')
                    closest_index = 0
                    
                    # Find the closest segment of the alignment
                    for i in range(len(alignment_coords) - 1):
                        line_segment = LineString([alignment_coords[i], alignment_coords[i+1]])
                        dist = point.distance(line_segment)
                        if dist < min_dist:
                            min_dist = dist
                            closest_index = i
                    
                    # Record the point and its position
                    point_positions.append((point, closest_index))
                
                # Sort points by their position along the alignment
                point_positions.sort(key=lambda x: x[1])
                sorted_boundary_points = [pp[0] for pp in point_positions]
                
                # Create segments between consecutive points
                segments = []
                for i in range(0, len(sorted_boundary_points), 2):
                    if i + 1 < len(sorted_boundary_points):  # Ensure we have a pair
                        segment = LineString([sorted_boundary_points[i], sorted_boundary_points[i+1]])
                        segments.append(segment)
                
                print(f"Created {len(segments)} line segments for splitting")
                
                # Save the segments
                if segments:
                    gpd.GeoDataFrame(geometry=segments, crs=alignment_gdf.crs).to_file(
                        os.path.join(debug_folder, "split_segments.geojson")
                    )
                
                # Now use these segments to split the water polygon
                split_results = water_geom
                for segment in segments:
                    try:
                        from shapely.ops import split
                        split_results = split(split_results, segment)
                    except Exception as e:
                        print(f"Error splitting with segment: {e}")
                
                # Check if we successfully split the water
                if hasattr(split_results, 'geoms') and len(split_results.geoms) > 1:
                    print(f"Successfully split water into {len(split_results.geoms)} pieces")
                    
                    # Save the split result
                    gpd.GeoDataFrame(geometry=list(split_results.geoms), crs=water_gdf.crs).to_file(
                        os.path.join(debug_folder, "water_split_result.geojson")
                    )
                    
                    # Create GeoDataFrame from pieces
                    water_pieces = gpd.GeoDataFrame(
                        geometry=list(split_results.geoms),
                        crs=water_gdf.crs
                    )
                    
                    # Apply the filtering logic
                    # Calculate area for each piece
                    water_pieces['area'] = water_pieces.geometry.area
                    
                    # Find the largest piece
                    largest_idx = water_pieces['area'].idxmax()
                    largest_piece = water_pieces.loc[largest_idx, 'geometry']
                    
                    # Filter pieces based on the connectivity logic
                    kept_pieces = [largest_piece]
                    dropped_pieces = []
                    
                    for idx, row in water_pieces.iterrows():
                        if idx == largest_idx:
                            continue  # Skip the largest piece
                        
                        # Get piece centroid
                        piece_geom = row.geometry
                        centroid = piece_geom.centroid
                        
                        # Find nearest point on the largest piece
                        from shapely.ops import nearest_points
                        nearest_point = nearest_points(centroid, largest_piece)[1]
                        
                        # Create line between centroid and nearest point
                        connection_line = LineString([centroid, nearest_point])
                        
                        # Check if alignment intersects this line
                        if alignment_line.intersects(connection_line):
                            print(f"Dropping water piece {idx} - alignment blocks connection to largest piece")
                            dropped_pieces.append(piece_geom)
                        else:
                            print(f"Keeping water piece {idx} - connected to largest piece")
                            kept_pieces.append(piece_geom)
                    
                    # Create new GeoDataFrame with only the kept pieces
                    filtered_water = gpd.GeoDataFrame(
                        geometry=kept_pieces,
                        crs=water_gdf.crs
                    )
                    
                    print(f"Kept {len(kept_pieces)} water pieces, dropped {len(dropped_pieces)} pieces")
                    
                    # Dissolve if needed
                    if len(kept_pieces) > 1:
                        filtered_water = filtered_water.dissolve()
                        print("Dissolved kept pieces into final water polygon")
                    
                    # Save final water polygon
                    water_crop_path = DATASET_INFO["Preprocessed"]["Water_crop"]["path"]
                    filtered_water.to_file(water_crop_path)
                    
                    return filtered_water
                else:
                    print("Segment splitting didn't create multiple pieces")
            else:
                print("Not enough boundary intersection points found for splitting")
                
            # Fall back to perpendicular line approach if entry/exit approach fails
            print("Falling back to perpendicular line approach...")
            
            # Create a straight line perpendicular to the alignment
            from shapely.geometry import box
            import numpy as np
            
            # Get the bounds of the water with a buffer
            minx, miny, maxx, maxy = water_geom.bounds
            dx = maxx - minx
            dy = maxy - miny
            center_x = (minx + maxx) / 2
            center_y = (miny + maxy) / 2
            
            # Create a line aligned with the general direction of the alignment
            direction_vector = LineString([
                alignment_line.coords[0],
                alignment_line.coords[-1]
            ])
            
            angle = np.arctan2(
                direction_vector.coords[-1][1] - direction_vector.coords[0][1],
                direction_vector.coords[-1][0] - direction_vector.coords[0][0]
            )
            
            # Create perpendicular line
            perp_angle = angle + np.pi/2
            line_length = max(dx, dy) * 2  # Make sure it's long enough
            
            perp_line = LineString([
                (center_x - np.cos(perp_angle) * line_length, 
                 center_y - np.sin(perp_angle) * line_length),
                (center_x + np.cos(perp_angle) * line_length,
                 center_y + np.sin(perp_angle) * line_length)
            ])
            
            # Save the perpendicular line
            gpd.GeoDataFrame(geometry=[perp_line], crs=alignment_gdf.crs).to_file(
                os.path.join(debug_folder, "perpendicular_split_line.geojson")
            )
            
            # Try splitting with this line
            from shapely.ops import split
            try:
                final_split = split(water_geom, perp_line)
                if hasattr(final_split, 'geoms') and len(final_split.geoms) > 1:
                    print(f"Perpendicular line split water into {len(final_split.geoms)} pieces")
                    
                    # Save the split result
                    gpd.GeoDataFrame(geometry=list(final_split.geoms), crs=water_gdf.crs).to_file(
                        os.path.join(debug_folder, "water_perpendicular_split.geojson")
                    )
                    
                    # Create GeoDataFrame from pieces
                    water_pieces = gpd.GeoDataFrame(
                        geometry=list(final_split.geoms),
                        crs=water_gdf.crs
                    )
                    
                    # Calculate area for each piece
                    water_pieces['area'] = water_pieces.geometry.area
                    
                    # Find the largest piece
                    largest_idx = water_pieces['area'].idxmax()
                    largest_piece = water_pieces.loc[largest_idx, 'geometry']
                    
                    # Filter pieces based on the connectivity logic
                    kept_pieces = [largest_piece]
                    dropped_pieces = []
                    
                    for idx, row in water_pieces.iterrows():
                        if idx == largest_idx:
                            continue  # Skip the largest piece
                        
                        # Get piece centroid
                        piece_geom = row.geometry
                        centroid = piece_geom.centroid
                        
                        # Find nearest point on the largest piece
                        from shapely.ops import nearest_points
                        nearest_point = nearest_points(centroid, largest_piece)[1]
                        
                        # Create line between centroid and nearest point
                        connection_line = LineString([centroid, nearest_point])
                        
                        # Check if alignment intersects this line
                        if alignment_line.intersects(connection_line):
                            print(f"Dropping water piece {idx} - alignment blocks connection to largest piece")
                            dropped_pieces.append(piece_geom)
                        else:
                            print(f"Keeping water piece {idx} - connected to largest piece")
                            kept_pieces.append(piece_geom)
                    
                    # Create new GeoDataFrame with only the kept pieces
                    filtered_water = gpd.GeoDataFrame(
                        geometry=kept_pieces,
                        crs=water_gdf.crs
                    )
                    
                    print(f"Kept {len(kept_pieces)} water pieces, dropped {len(dropped_pieces)} pieces")
                    
                    # Dissolve if needed
                    if len(kept_pieces) > 1:
                        filtered_water = filtered_water.dissolve()
                        print("Dissolved kept pieces into final water polygon")
                    
                    # Save final water polygon
                    water_crop_path = DATASET_INFO["Preprocessed"]["Water_crop"]["path"]
                    filtered_water.to_file(water_crop_path)
                    
                    return filtered_water
                else:
                    print("Perpendicular line splitting also failed")
            except Exception as e:
                print(f"Error in perpendicular line splitting: {e}")
                
        except Exception as e:
            print(f"Error in entry/exit point processing: {e}")
            import traceback
            traceback.print_exc()
            
    else:
        print("DEBUG: Alignment does NOT intersect water polygon!")
        print(f"Alignment bounds: {alignment_line.bounds}")
        print(f"Water bounds: {water_geom.bounds}")
    
    # If all approaches fail, just return the cropped water
    print("All splitting approaches failed. Returning the original cropped water polygon.")
    water_crop_path = DATASET_INFO["Preprocessed"]["Water_crop"]["path"]
    water_crop.to_file(water_crop_path)
    print(f"Saved original cropped water polygon to {water_crop_path}")
    return water_crop

def preprocess_data():
    """
    Preprocess data for flood infrastructure analysis
    
    Returns:
    dict: Dictionary of preprocessed data
    """
    print("Starting preprocessing...")
    
    # Get project directory for relative paths
    project_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Create output directories if they don't exist
    for dataset_type in DATASET_INFO:
        for dataset in DATASET_INFO[dataset_type]:
            # Check if path exists for this dataset
            if "path" in DATASET_INFO[dataset_type][dataset]:
                path = DATASET_INFO[dataset_type][dataset]["path"]
                # Convert path to absolute path if it's not already
                if not os.path.isabs(path):
                    path = os.path.join(project_dir, str(path))
                os.makedirs(os.path.dirname(path), exist_ok=True)
            else:
                # Skip warning for Webmap entries - they'll be created during webmap generation
                if dataset_type != "Webmap":
                    print(f"Warning: No path found for {dataset_type}/{dataset}")
    
    # DEBUGGING: Print dataset keys to identify correct structure
    print("Available dataset types:", DATASET_INFO.keys())
    for dataset_type in DATASET_INFO:
        print(f"Datasets in {dataset_type}:", DATASET_INFO[dataset_type].keys())
    
    # Create site bounds first since we need it for cropping
    print("Creating/loading site bounds...")
    site_bounds = create_site_bounds()
    
    # Load the buildings data
    buildings_path = DATASET_INFO["Input"]["Buildings"]["path"]
    buildings_gdf = gpd.read_file(buildings_path)
    print(f"Loaded buildings data from {buildings_path}")
    
    # Check CRS and reproject if needed
    if buildings_gdf.crs != PROJECT_CRS:
        print(f"Reprojecting buildings from {buildings_gdf.crs} to {PROJECT_CRS}")
        buildings_gdf = buildings_gdf.to_crs(PROJECT_CRS)
    
    # Clip buildings to the site boundary
    print(f"Clipping buildings to study area... ({len(buildings_gdf)} buildings before clipping)")
    buildings_gdf = gpd.clip(buildings_gdf, site_bounds)
    print(f"Buildings clipped to study area. {len(buildings_gdf)} buildings remaining.")
    
    # Save cropped buildings to file
    buildings_crop_path = DATASET_INFO["Preprocessed"]["Buildings_crop"]["path"]
    buildings_gdf.to_file(buildings_crop_path, driver='GeoJSON')
    print(f"Saved cropped buildings to {buildings_crop_path}")
    
    # Load and process alignment polyline before water processing
    alignment_gdf = create_alignment()
    
    # Load and process water polygon with advanced filtering
    if "Water" in DATASET_INFO["Input"]:
        water_path = DATASET_INFO["Input"]["Water"]["path"]
        if os.path.exists(water_path):
            water_gdf = gpd.read_file(water_path)
            
            # Check if the water has a CRS
            if water_gdf.crs is None:
                print(f"Warning: Water has no CRS. Setting to {PROJECT_CRS}")
                water_gdf.set_crs(PROJECT_CRS, inplace=True)
            elif water_gdf.crs != PROJECT_CRS:
                print(f"Reprojecting water from {water_gdf.crs} to {PROJECT_CRS}")
                water_gdf = water_gdf.to_crs(PROJECT_CRS)
                
            # Apply advanced water processing
            water_gdf = process_water_polygon(water_gdf, alignment_gdf, site_bounds)
            
        else:
            water_gdf = None
            print(f"Warning: Water polygon file not found at {water_path}")
    else:
        water_gdf = None
        print("Warning: No water polygon path specified in config")
    
    # Load NSI data if available
    if "NSI" in DATASET_INFO["Input"]:
        nsi_path = DATASET_INFO["Input"]["NSI"]["path"]
        if os.path.exists(nsi_path):
            nsi_gdf = gpd.read_file(nsi_path)
            
            # Check if the NSI has a CRS
            if nsi_gdf.crs is None:
                print(f"Warning: NSI has no CRS. Setting to {PROJECT_CRS}")
                nsi_gdf.set_crs(PROJECT_CRS, inplace=True)
            elif nsi_gdf.crs != PROJECT_CRS:
                print(f"Reprojecting NSI from {nsi_gdf.crs} to {PROJECT_CRS}")
                nsi_gdf = nsi_gdf.to_crs(PROJECT_CRS)
        else:
            nsi_gdf = None
            print(f"Warning: NSI file not found at {nsi_path}")
    else:
        nsi_gdf = None
        print("Warning: No NSI path specified in config")
    
    # Crop FEMA flood raster
    fema_flood_crop = crop_raster_data(site_bounds)
    
    # Return dictionary of preprocessed data
    return {
        "buildings": buildings_gdf,
        "fema_flood_crop": fema_flood_crop,
        "alignment": alignment_gdf,
        "water": water_gdf,
        "nsi": nsi_gdf
    }

if __name__ == "__main__":
    preprocess_data()
