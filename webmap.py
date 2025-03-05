"""
Webmap script for Flood Infrastructure Analysis
Creates an interactive Folium map with flood protection scenarios
"""
import os
import folium
from folium import plugins
import geopandas as gpd
import numpy as np
import rasterio
from rasterio import warp
from rasterio.enums import Resampling
from affine import Affine
from PIL import Image
import webbrowser
from pathlib import Path
from config import (
    DATASET_INFO, WEB_CRS, INITIAL_ZOOM, RESOLUTION,
    LEGEND_STYLES, TITLE_STYLE, INFO_STYLE, LOGO_STYLE,
    hex_to_rgb, hex_to_rgba, interpolate_color, interpolate_color_with_alpha,
    OVERWRITE
)

def process_raster_for_web(input_raster, output_png, target_crs=WEB_CRS, colormap="flood"):
    """
    Reprojects, downsamples, and applies a colormap to a raster,
    then saves it as a PNG for use as an ImageOverlay.
    """
    try:
        # Ensure input_raster and output_png are strings, not Path objects
        input_raster = str(input_raster)
        output_png = str(output_png)
        
        # Verify the input file exists
        if not os.path.exists(input_raster):
            print(f"ERROR: Input raster file does not exist: {input_raster}")
            raise FileNotFoundError(f"Input file not found: {input_raster}")
            
        # Reproject the raster
        with rasterio.open(input_raster) as src:
            # Store original bounds for accurate placement
            src_bounds = src.bounds
            
            # Get reprojected bounds for exact positioning
            dst_bounds = warp.transform_bounds(src.crs, target_crs, *src.bounds)
            
            print(f"Processing {input_raster}")
            print(f"  Original bounds (source CRS): {src_bounds}")
            print(f"  Transformed bounds (target CRS): {dst_bounds}")
            print(f"  Using resolution: {RESOLUTION} feet")
            
            # Verify the raster has valid bounds
            if not all(np.isfinite(b) for b in dst_bounds):
                print(f"ERROR: Invalid bounds after reprojection: {dst_bounds}")
                print("Attempting to use original bounds as fallback")
                try:
                    dst_bounds = src_bounds
                    if not all(np.isfinite(b) for b in dst_bounds):
                        raise ValueError("Source bounds are also invalid")
                except:
                    print("Failed to use original bounds. Using NYC extent as fallback.")
                    dst_bounds = (-74.26, 40.49, -73.69, 40.91)
            
            # Calculate dimensions at specified resolution
            try:
                if target_crs == "EPSG:4326":
                    # Approximate conversion at 40.7° N:
                    # 1 degree latitude ≈ 364,320 feet
                    # 1 degree longitude ≈ 364,320 * cos(40.7) feet
                    deg_per_foot_lat = 1 / 364320
                    deg_per_foot_lng = 1 / (364320 * np.cos(np.radians(40.7)))
                    
                    # Compute pixel size in degrees for the given resolution (in feet)
                    pixel_size_deg_lat = RESOLUTION * deg_per_foot_lat
                    pixel_size_deg_lng = RESOLUTION * deg_per_foot_lng
                    
                    # Calculate target dimensions based on the converted resolution
                    width = max(100, int((dst_bounds[2] - dst_bounds[0]) / pixel_size_deg_lng))
                    height = max(100, int((dst_bounds[3] - dst_bounds[1]) / pixel_size_deg_lat))
                else:
                    width = max(100, int((dst_bounds[2] - dst_bounds[0]) / RESOLUTION))
                    height = max(100, int((dst_bounds[3] - dst_bounds[1]) / RESOLUTION))
                
                print(f"  Target dimensions: {width}x{height} pixels")
            except Exception as e:
                print(f"Error calculating dimensions: {e}")
                width, height = 1000, 1000
                print(f"  Using fallback dimensions: {width}x{height} pixels")
    
            # Calculate reprojection transform
            transform = Affine.translation(dst_bounds[0], dst_bounds[3]) * Affine.scale(
                (dst_bounds[2] - dst_bounds[0]) / width, (dst_bounds[1] - dst_bounds[3]) / height
            )
            
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': target_crs,
                'transform': transform,
                'width': width,
                'height': height
            })
            
            data = np.empty((src.count, height, width), dtype=src.dtypes[0])
            
            try:
                for i in range(1, src.count + 1):
                    warp.reproject(
                        source=rasterio.band(src, i),
                        destination=data[i-1],
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=target_crs,
                        resampling=Resampling.bilinear
                    )
                nodata = src.nodata
            except Exception as e:
                print(f"Error during reprojection: {e}")
                data = np.zeros((src.count, height, width), dtype=np.float32)
                nodata = None

        # Use the first band for visualization
        band = data[0]
        
        # Handle nodata values
        if nodata is not None:
            valid = (band != nodata) & (~np.isnan(band))
        else:
            valid = ~np.isnan(band)
        
        # Debug: Print out the range of valid data values
        valid_data = band[valid]
        
        # Create an empty RGBA array
        rgba = np.zeros((height, width, 4), dtype=np.uint8)
        
        # Select colors for FEMA flood data
        if "FEMA" in input_raster:
            hex_1pct = DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_1pct"]
            hex_02pct = DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_0_2pct"]
            opacity = DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["opacity"]
            
            for i in range(height):
                for j in range(width):
                    if not valid[i, j] or band[i, j] == 0:
                        rgba[i, j] = [0, 0, 0, 0]  # Fully transparent
                    elif band[i, j] == 2:  # 1% annual chance (100-year flood)
                        r, g, b = hex_to_rgb(hex_1pct)
                        alpha = int(opacity * 255)
                        rgba[i, j] = [r, g, b, alpha]
                    elif band[i, j] == 1:  # 0.2% annual chance (500-year flood)
                        r, g, b = hex_to_rgb(hex_02pct)
                        alpha = int(opacity * 255)
                        rgba[i, j] = [r, g, b, alpha]
                    else:
                        rgba[i, j] = [0, 0, 0, 0]
        
        img = Image.fromarray(rgba)
        img.save(output_png)
        print(f"  Successfully saved {output_png}")
        
        if dst_bounds is None or not all(np.isfinite(b) for b in dst_bounds):
            print("WARNING: Invalid bounds detected. Using NYC fallback bounds.")
            dst_bounds = (-74.26, 40.49, -73.69, 40.91)
        
        return output_png, dst_bounds
        
    except Exception as e:
        print(f"Error processing raster {input_raster}: {e}")
        print("Creating a fallback image to ensure visibility.")
        
        # Create a distinguishable fallback image
        img_size = 500
        fallback_img = np.zeros((img_size, img_size, 4), dtype=np.uint8)
        
        # Create a pattern that will be clearly visible and identifiable
        # Red border to show it's a fallback
        fallback_img[0:10, :] = [255, 0, 0, 255]
        fallback_img[-10:, :] = [255, 0, 0, 255]
        fallback_img[:, 0:10] = [255, 0, 0, 255]
        fallback_img[:, -10:] = [255, 0, 0, 255]
        
        # Add some blue patterns inside
        if "trim" in input_raster:
            # For trimmed raster, make a different pattern
            for i in range(50, 450, 100):
                fallback_img[i:i+50, 50:450] = [0, 0, 255, 128]
        else:
            # For original raster
            for i in range(50, 450, 50):
                fallback_img[i:i+25, 50:450] = [0, 0, 255, 200]
        
        # Save the test image
        Image.fromarray(fallback_img).save(output_png)
        print(f"Fallback image saved to {output_png}")
        
        print("Using NYC fallback bounds.")
        return output_png, (-74.26, 40.49, -73.69, 40.91)

def style_function(feature, scenario="existing"):
    """
    Style function for GeoJSON layers
    
    Parameters:
    feature: GeoJSON feature
    scenario: 'existing' or 'alignment'
    
    Returns:
    dict: Style dictionary for the feature
    """
    properties = feature.get('properties', {})
    
    # Print debugging info for the first feature
    if hasattr(style_function, 'first_call') and not style_function.first_call:
        print(f"Feature properties: {properties}")
        style_function.first_call = True
    
    # Determine if the building is flooded
    if scenario == "existing":
        is_flooded = properties.get('0_2_PctFlood', False)
    else:
        is_flooded = properties.get('0_2_PctFlood_Alignment', False)
    
    # Set style based on flood status
    if is_flooded:
        return {
            'fillColor': DATASET_INFO["Webmap"]["Buildings"]["flooded"],
            'color': '#000000',
            'weight': 1,
            'fillOpacity': DATASET_INFO["Webmap"]["Buildings"]["opacity"]
        }
    else:
        return {
            'fillColor': DATASET_INFO["Webmap"]["Buildings"]["not_flooded"],
            'color': '#000000',
            'weight': 0.5,
            'fillOpacity': DATASET_INFO["Webmap"]["Buildings"]["opacity"] * 0.7
        }

# Initialize the first_call flag
style_function.first_call = False
def create_legend(m):
    """
    Create a legend for the map
    
    Parameters:
    m: Folium map object
    """
    legend_html = f'''
    <div style="{LEGEND_STYLES['container']}">
        <h4 style="{LEGEND_STYLES['header']}">Legend</h4>
        
        <h5 style="{LEGEND_STYLES['sectionHeader']}">Buildings</h5>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["Buildings"]["flooded"]}"></div>
            <span style="{LEGEND_STYLES['label']}">Flooded Buildings</span>
        </div>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["Buildings"]["not_flooded"]}"></div>
            <span style="{LEGEND_STYLES['label']}">Non-Flooded Buildings</span>
        </div>
        
        <h5 style="{LEGEND_STYLES['sectionHeader']}">FEMA Flood Zones</h5>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_1pct"]}"></div>
            <span style="{LEGEND_STYLES['label']}">1% Annual Chance (100-yr)</span>
        </div>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_0_2pct"]}"></div>
            <span style="{LEGEND_STYLES['label']}">0.2% Annual Chance (500-yr)</span>
        </div>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(legend_html))

def add_title(m):
    """
    Add a title to the map
    
    Parameters:
    m: Folium map object
    """
    title_html = f'''
    <div style="{TITLE_STYLE['container']}">
        <h1 style="{TITLE_STYLE['title']}">Flood Protection Scenarios</h1>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(title_html))

def add_info_box(m, count_existing, count_alignment):
    """
    Add info box with flooded building counts
    
    Parameters:
    m: Folium map object
    count_existing: Count of flooded buildings in existing scenario
    count_alignment: Count of flooded buildings in alignment scenario
    """
    info_html = f'''
    <div id="info-box" style="{INFO_STYLE['container']}">
        <p style="{INFO_STYLE['text']}">Flooded buildings: <span id="flood-count">{count_existing}</span></p>
    </div>
    
    <script>
        // Update count when layer visibility changes
        document.addEventListener('DOMContentLoaded', function() {{
            // Wait for map layers to be initialized
            setTimeout(function() {{
                var existingButton = document.querySelector('input[type="radio"][name="leaflet-base-layers"][value="Existing Scenario"]');
                var alignmentButton = document.querySelector('input[type="radio"][name="leaflet-base-layers"][value="Alignment Scenario"]');
                var countElement = document.getElementById('flood-count');
                
                if (existingButton && alignmentButton && countElement) {{
                    // Set up event listeners
                    existingButton.addEventListener('change', function() {{
                        if (this.checked) {{
                            countElement.textContent = "{count_existing}";
                        }}
                    }});
                    
                    alignmentButton.addEventListener('change', function() {{
                        if (this.checked) {{
                            countElement.textContent = "{count_alignment}";
                        }}
                    }});
                }}
            }}, 1000);
        }});
    </script>
    '''
    
    m.get_root().html.add_child(folium.Element(info_html))

def add_logo(m):
    """
    Add ONE Analysis logo to the map
    
    Parameters:
    m: Folium map object
    """
    logo_html = f'''
    <div style="{LOGO_STYLE['container']}" class="analysis-text">
        <span style="font-family: 'Futura', sans-serif; font-weight: bold; color: #4c5da4;">one</span>
        <span style="font-family: 'Futura', sans-serif; font-weight: 300; color: #4c5da4;"> analysis</span>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(logo_html))

def create_webmap(analysis_results):
    """
    Create a Folium webmap with flood protection scenarios
    
    Parameters:
    analysis_results (dict): Dictionary containing analysis results
    
    Returns:
    str: Path to the generated HTML file
    """
    buildings_scenarios = analysis_results["buildings_scenarios"]
    
    # Print debugging info
    print(f"Number of buildings: {len(buildings_scenarios)}")
    print(f"CRS of buildings data: {buildings_scenarios.crs}")
    print(f"Bounds of buildings data: {buildings_scenarios.total_bounds}")
    
    # Get paths to raster data
    fema_flood_path = str(DATASET_INFO["Output"]["FEMA_Flood_crop"]["path"])
    fema_flood_trim_path = str(DATASET_INFO["Output"]["FEMA_Flood_crop_trim"]["path"])
    
    # Check if the raster files exist
    print(f"FEMA flood path exists: {os.path.exists(fema_flood_path)}")
    print(f"FEMA flood trim path exists: {os.path.exists(fema_flood_trim_path)}")
    
    # Process rasters for web with more verbose logging
    fema_web_path = str(DATASET_INFO["Output"]["FEMA_Flood_crop_web"]["path"])
    fema_trim_web_path = str(DATASET_INFO["Output"]["FEMA_Flood_crop_trim_web"]["path"])
    
    print("Processing original FEMA flood raster...")
    fema_web_path, fema_bounds = process_raster_for_web(str(fema_flood_path), fema_web_path)
    print(f"FEMA bounds: {fema_bounds}")
    
    print("Processing trimmed FEMA flood raster...")
    fema_trim_web_path, fema_trim_bounds = process_raster_for_web(str(fema_flood_trim_path), fema_trim_web_path)
    print(f"FEMA trim bounds: {fema_trim_bounds}")
    
    # Ensure the web images exist
    print(f"FEMA web image exists: {os.path.exists(fema_web_path)}")
    print(f"FEMA trim web image exists: {os.path.exists(fema_trim_web_path)}")
    
    # Reproject buildings to WGS84 for web mapping
    print("Reprojecting buildings to WGS84...")
    buildings_wgs84 = buildings_scenarios.to_crs(WEB_CRS)
    
    # Get center point for map initialization
    if not buildings_wgs84.empty:
        bounds = buildings_wgs84.total_bounds
        center = [
            (bounds[1] + bounds[3]) / 2,  # average of min and max y (latitude)
            (bounds[0] + bounds[2]) / 2   # average of min and max x (longitude)
        ]
    else:
        # Default to a New York City location if no buildings data
        center = [40.7128, -74.0060]
        
    print(f"Map center: {center}")
        
    # Create Folium map
    m = folium.Map(
        location=center,
        zoom_start=INITIAL_ZOOM,
        tiles='CartoDB Positron'
    )
    
    # If we have valid bounds, make sure the map is zoomed to show all data
    if not buildings_wgs84.empty:
        # Convert shapely/geopandas bounds to folium bounds [[lat_min, lon_min], [lat_max, lon_max]]
        folium_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
        m.fit_bounds(folium_bounds)
        print(f"Setting map bounds to: {folium_bounds}")
    
    # Add layer control
    folium.LayerControl(position='topright', collapsed=False).add_to(m)
    
    # Calculate flooded building counts
    if not buildings_wgs84.empty:
        count_existing = int(buildings_wgs84['0_2_PctFlood'].sum())
        count_alignment = int(buildings_wgs84['0_2_PctFlood_Alignment'].sum())
    else:
        count_existing = 0
        count_alignment = 0
    
    # Create layer groups for scenarios
    existing_fg = folium.FeatureGroup(name='Existing Scenario', overlay=False, show=True)
    alignment_fg = folium.FeatureGroup(name='Alignment Scenario', overlay=False, show=False)
    
    # Add FEMA rasters to layer groups if the images exist and bounds are valid
    if os.path.exists(fema_web_path) and all(np.isfinite(b) for b in fema_bounds):
        print(f"Adding FEMA flood image overlay with bounds: {fema_bounds}")
        folium.raster_layers.ImageOverlay(
            image=fema_web_path,
            bounds=[[fema_bounds[1], fema_bounds[0]], [fema_bounds[3], fema_bounds[2]]],
            opacity=0.7,  # Increase opacity
            name='FEMA Flood Zones'
        ).add_to(existing_fg)
    else:
        print("Warning: Could not add FEMA flood image overlay - missing image or invalid bounds")
    
    if os.path.exists(fema_trim_web_path) and all(np.isfinite(b) for b in fema_trim_bounds):
        print(f"Adding FEMA trim flood image overlay with bounds: {fema_trim_bounds}")
        folium.raster_layers.ImageOverlay(
            image=fema_trim_web_path,
            bounds=[[fema_trim_bounds[1], fema_trim_bounds[0]], [fema_trim_bounds[3], fema_trim_bounds[2]]],
            opacity=0.7,  # Increase opacity
            name='FEMA Flood Zones (with Alignment)'
        ).add_to(alignment_fg)
    else:
        print("Warning: Could not add FEMA trim flood image overlay - missing image or invalid bounds")
    
    # Add buildings to layer groups using style functions
    if not buildings_wgs84.empty:
        # First, check for required flood columns
        has_flood_cols = '0_2_PctFlood' in buildings_wgs84.columns and '0_2_PctFlood_Alignment' in buildings_wgs84.columns
        print(f"Building data has required flood columns: {has_flood_cols}")
        if not has_flood_cols:
            print(f"Available columns: {buildings_wgs84.columns.tolist()}")
            
        # Add missing columns if needed
        if '0_2_PctFlood' not in buildings_wgs84.columns:
            buildings_wgs84['0_2_PctFlood'] = False
        if '0_2_PctFlood_Alignment' not in buildings_wgs84.columns:
            buildings_wgs84['0_2_PctFlood_Alignment'] = False
            
        # Check if any buildings are flooded
        print(f"Buildings in flood zone (existing): {buildings_wgs84['0_2_PctFlood'].sum()}")
        print(f"Buildings in flood zone (alignment): {buildings_wgs84['0_2_PctFlood_Alignment'].sum()}")
        
        buildings_geojson = buildings_wgs84.to_json()
        
        folium.GeoJson(
            buildings_geojson,
            name='Buildings (Existing)',
            style_function=lambda x: style_function(x, "existing"),
            tooltip=folium.GeoJsonTooltip(
                fields=['1_PctFlood', '0_2_PctFlood'],
                aliases=['1% Annual Chance Flood:', '0.2% Annual Chance Flood:'],
                style=("background-color: white; color: #333333; font-family: arial; font-size: 12px; padding: 10px;")
            )
        ).add_to(existing_fg)
        
        folium.GeoJson(
            buildings_geojson,
            name='Buildings (with Alignment)',
            style_function=lambda x: style_function(x, "alignment"),
            tooltip=folium.GeoJsonTooltip(
                fields=['1_PctFlood_Alignment', '0_2_PctFlood_Alignment'],
                aliases=['1% Annual Chance Flood:', '0.2% Annual Chance Flood:'],
                style=("background-color: white; color: #333333; font-family: arial; font-size: 12px; padding: 10px;")
            )
        ).add_to(alignment_fg)
    else:
        print("Warning: No building data to display")
    
    # Add layer groups to map
    existing_fg.add_to(m)
    alignment_fg.add_to(m)
    
    # Add legend, title, info box, and logo
    create_legend(m)
    add_title(m)
    add_info_box(m, count_existing, count_alignment)
    add_logo(m)
    
    # Save map to HTML file
    output_path = str(DATASET_INFO["Output"]["Webmap"]["path"])
    m.save(output_path)
    print(f"Web map saved to {output_path}")
    
    # Add a debug map to help diagnose issues
    debug_path = os.path.join(os.path.dirname(output_path), "debug_map.html")
    
    # Create a simplified debug map just showing the buildings
    debug_map = folium.Map(location=center, zoom_start=INITIAL_ZOOM)
    
    # First, add a tile layer with markers at the corners of the bounding box
    if not buildings_wgs84.empty:
        bounds = buildings_wgs84.total_bounds
        folium.Marker([bounds[1], bounds[0]], popup="SW Corner").add_to(debug_map)
        folium.Marker([bounds[1], bounds[2]], popup="SE Corner").add_to(debug_map)
        folium.Marker([bounds[3], bounds[0]], popup="NW Corner").add_to(debug_map)
        folium.Marker([bounds[3], bounds[2]], popup="NE Corner").add_to(debug_map)
        
        # Add flood raster bounding box markers
        if all(np.isfinite(b) for b in fema_bounds):
            folium.Marker([fema_bounds[1], fema_bounds[0]], popup="FEMA SW", icon=folium.Icon(color='blue')).add_to(debug_map)
            folium.Marker([fema_bounds[3], fema_bounds[2]], popup="FEMA NE", icon=folium.Icon(color='blue')).add_to(debug_map)
        
        # Add all buildings as red polygons for visibility
        folium.GeoJson(
            buildings_wgs84,
            name='All Buildings',
            style_function=lambda x: {'fillColor': 'red', 'color': 'black', 'weight': 1, 'fillOpacity': 0.8}
        ).add_to(debug_map)
        
        # Fit to bounds
        folium_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
        debug_map.fit_bounds(folium_bounds)
    
    # Add a special marker at the center for reference
    folium.Marker(center, popup="Map Center", icon=folium.Icon(color='green')).add_to(debug_map)
    
    # Create a visual test overlay that should definitely be visible
    test_bounds = [[center[0] - 0.01, center[1] - 0.01], [center[0] + 0.01, center[1] + 0.01]]
    img = np.zeros((100, 100, 4), dtype=np.uint8)
    img[25:75, 25:75] = [255, 0, 0, 200]  # Red square in the center
    test_img_path = os.path.join(os.path.dirname(output_path), "test_overlay.png")
    Image.fromarray(img).save(test_img_path)
    
    folium.raster_layers.ImageOverlay(
        image=test_img_path,
        bounds=test_bounds,
        opacity=0.8,
        name='Test Overlay (Red Square)'
    ).add_to(debug_map)
    
    debug_map.save(debug_path)
    print(f"Debug map saved to {debug_path}")
    
    return output_path

def generate_webmap(analysis_results):
    """
    Generate the webmap if it doesn't exist or if overwrite is enabled
    
    Parameters:
    analysis_results (dict): Dictionary containing analysis results
    
    Returns:
    str: Path to the generated HTML file
    """
    output_path = str(DATASET_INFO["Output"]["Webmap"]["path"])
    
    # Check if webmap already exists and respect OVERWRITE setting
    if os.path.exists(output_path) and not OVERWRITE["webmap"]:
        print(f"Web map already exists at {output_path}")
        return output_path
    
    return create_webmap(analysis_results)

if __name__ == "__main__":
    # This would normally be called from main.py
    from preprocessing import preprocess_data
    from analysis import run_analysis
    
    preprocessed_data = preprocess_data()
    analysis_results = run_analysis(preprocessed_data)
    webmap_path = generate_webmap(analysis_results)
    
    # Open the webmap in a browser
    webbrowser.open('file://' + os.path.abspath(webmap_path))