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
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_png), exist_ok=True)
        
        # Verify the input file exists
        if not os.path.exists(input_raster):
            print(f"ERROR: Input raster file does not exist: {input_raster}")
            raise FileNotFoundError(f"Input file not found: {input_raster}")
            
        # Reproject the raster - cap the size to prevent huge files
        with rasterio.open(input_raster) as src:
            # Store original bounds for accurate placement
            src_bounds = src.bounds
            
            # Get reprojected bounds for exact positioning
            dst_bounds = warp.transform_bounds(src.crs, target_crs, *src.bounds)
            
            print(f"Processing {input_raster}")
            print(f"  Original CRS: {src.crs}")
            print(f"  Target CRS: {target_crs}")
            print(f"  Original bounds (source CRS): {src_bounds}")
            print(f"  Transformed bounds (target CRS): {dst_bounds}")
            print(f"  Using resolution: {RESOLUTION} feet")
            
            # Calculate target width and height with max size limit
            MAX_SIZE = 2000  # Max dimension in pixels
            width = int((dst_bounds[2] - dst_bounds[0]) / (RESOLUTION / 3.28084 / 111000) * 10)
            height = int((dst_bounds[3] - dst_bounds[1]) / (RESOLUTION / 3.28084 / 111000) * 10)
            
            # Scale down if needed
            if width > MAX_SIZE or height > MAX_SIZE:
                scale = max(width / MAX_SIZE, height / MAX_SIZE)
                width = int(width / scale)
                height = int(height / scale)
            
            print(f"  Target dimensions: {width}x{height} pixels")
            
            # Reproject to WGS84
            dst_crs = {'init': target_crs} if target_crs.startswith('EPSG:') else target_crs
            data, transform = warp.reproject(
                source=src.read(),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_crs=dst_crs,
                dst_shape=(height, width),
                resampling=warp.Resampling.nearest
            )
            
            # Squeeze data to get 2D array
            data = np.squeeze(data)
            
            # Apply colormap based on the data type
            valid_min = int(np.min(data[data > 0])) if np.any(data > 0) else 1
            valid_max = int(np.max(data)) if np.any(data > 0) else 2
            print(f"  Valid data range: {valid_min} to {valid_max}")
            unique_values = np.unique(data[data > 0])
            print(f"  Unique values: {unique_values}")
            
            # Convert data to RGBA
            rgba = np.zeros((height, width, 4), dtype=np.uint8)
            
            # Color scheme for flood zones
            if colormap == "flood":
                # 0: Transparent (no data)
                # 1: Blue (1% annual chance flood zone / 100-year)
                # 2: Light blue (0.2% annual chance flood zone / 500-year)
                data_colors = {
                    0: [0, 0, 0, 0],     # Transparent
                    1: [0, 0, 255, 150], # Blue with 60% opacity
                    2: [150, 200, 255, 125]  # Light blue with 50% opacity
                }
                
                # Apply the colors
                for val, color in data_colors.items():
                    rgba[data == val] = color
            else:
                # Generic colormap
                for val in np.unique(data):
                    if val == 0:
                        rgba[data == val] = [0, 0, 0, 0]  # Transparent
                    else:
                        # Linear interpolation from blue to red
                        t = (val - valid_min) / (valid_max - valid_min) if valid_max > valid_min else 0
                        r = int(0 + t * 255)
                        g = int(0 + (1-t) * 255)
                        b = int(255 - t * 255)
                        rgba[data == val] = [r, g, b, 200]
            
            # Save as PNG
            Image.fromarray(rgba).save(output_png)
            print(f"  Successfully saved {output_png}")
            
            # Return the path to the PNG and the bounds
            return output_png, dst_bounds
            
    except Exception as e:
        print(f"Error processing raster for web: {e}")
        
        # Create a fallback image - reduced in size
        width, height = 500, 500
        
        # Create a blank RGBA image, small grid pattern
        fallback_img = np.zeros((height, width, 4), dtype=np.uint8)
        fallback_img[::20, :, :] = [150, 150, 150, 100]  # Horizontal lines
        fallback_img[:, ::20, :] = [150, 150, 150, 100]  # Vertical lines
        
        # Add some colored regions to indicate an error
        for i in range(50, 450, 50):
            fallback_img[i:i+25, 50:450] = [255, 0, 0, 200]  # Red bars
        
        # Save the test image
        os.makedirs(os.path.dirname(output_png), exist_ok=True)
        Image.fromarray(fallback_img).save(output_png)
        print(f"Fallback image saved to {output_png}")
        
        print("Using NYC fallback bounds.")
        return output_png, (-74.05, 40.61, -73.85, 40.81)  # Smaller area for testing

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
        print(f"Using flooded color: {DATASET_INFO['Webmap']['Buildings']['flooded']}")
        print(f"Using non-flooded color: {DATASET_INFO['Webmap']['Buildings']['not_flooded']}")
        print(f"Using opacity: {DATASET_INFO['Webmap']['Buildings']['opacity']}")
        style_function.first_call = True
    
    # Determine if the building is flooded
    if scenario == "existing":
        # Look for different possible column names
        for col in ['flood_02pct', '0_2_PctFlood', 'is_flooded']:
            if col in properties:
                is_flooded = properties.get(col, False)
                if is_flooded:
                    break
        else:
            is_flooded = False
    else:
        # Look for different possible column names for scenario
        for col in ['flood_02pct_scenario1', 'scenario1', 'is_protected']:
            if col in properties:
                is_flooded = properties.get(col, False)
                if is_flooded:
                    break
        else:
            is_flooded = False
    
    # Set style based on flood status - use colors directly from config
    if is_flooded:
        return {
            'fillColor': DATASET_INFO["Webmap"]["Buildings"]["flooded"],  # Red for flooded buildings
            'color': '#000000',
            'weight': 1,
            'fillOpacity': DATASET_INFO["Webmap"]["Buildings"]["opacity"]
        }
    else:
        return {
            'fillColor': DATASET_INFO["Webmap"]["Buildings"]["not_flooded"],  # Dark grey for non-flooded buildings
            'color': '#000000',
            'weight': 0.5,
            'fillOpacity': DATASET_INFO["Webmap"]["Buildings"]["opacity"] * 0.7
        }

def alignment_style_function(feature):
    """Style function for alignment polyline"""
    return {
        'color': '#000000',  # Solid black
        'weight': 3,         # Moderate thickness
        'opacity': 1.0,      # Full opacity
        'dashArray': None,   # Solid line
        'fillOpacity': 0
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
        
        <h5 style="{LEGEND_STYLES['sectionHeader']}">FEMA Flood Zones</h5>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_1pct"]}"></div>
            <span style="{LEGEND_STYLES['label']}">1% Annual Chance (100-yr)</span>
        </div>
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: {DATASET_INFO["Webmap"]["FEMA_FloodHaz"]["hex_0_2pct"]}"></div>
            <span style="{LEGEND_STYLES['label']}">0.2% Annual Chance (500-yr)</span>
        </div>
        
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: #000000; height: 3px; margin-top: 8px;"></div>
            <span style="{LEGEND_STYLES['label']}">Flood Protection Alignment</span>
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

def add_info_box(m, count_existing, count_protected):
    """
    Add an info box with protection statistics
    
    Parameters:
    m: Folium map object
    count_existing: Number of buildings flooded in existing scenario
    count_protected: Number of buildings remaining flooded with protection
    """
    # Calculate the percent reduction
    percent_reduction = 0
    if count_existing > 0:
        percent_reduction = 100 - (count_protected / count_existing * 100)
    
    info_html = f'''
    <div style="{INFO_STYLE['container']}">
        <p style="{INFO_STYLE['text']}">
            {count_existing} buildings flooded in existing conditions<br>
            {count_protected} buildings remain flooded with protection<br>
            {percent_reduction:.1f}% reduction in flooded buildings
        </p>
    </div>
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

def rasterize_buildings(buildings_gdf, output_png, bounds, width=1000, height=1000, flooded_column='0_2_PctFlood'):
    """
    Rasterize buildings for web visualization
    
    Parameters:
    buildings_gdf: GeoDataFrame of buildings
    output_png: Path to save the PNG output
    bounds: Bounds of the area to rasterize [minx, miny, maxx, maxy]
    width, height: Dimensions of the output raster
    flooded_column: Column indicating flood status
    
    Returns:
    str: Path to the PNG file
    """
    print(f"Rasterizing buildings with {flooded_column} to {output_png}")
    
    # Create an empty RGBA array
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    
    # Get colors from config
    flooded_color = hex_to_rgb(DATASET_INFO["Webmap"]["Buildings"]["flooded"])
    non_flooded_color = hex_to_rgb(DATASET_INFO["Webmap"]["Buildings"]["not_flooded"])
    opacity = int(DATASET_INFO["Webmap"]["Buildings"]["opacity"] * 255)
    
    try:
        # Convert bounds to pixel coordinates
        def world_to_pixel(x, y):
            px = int((x - bounds[0]) / (bounds[2] - bounds[0]) * width)
            py = int((bounds[3] - y) / (bounds[3] - bounds[1]) * height)
            return max(0, min(px, width-1)), max(0, min(py, height-1))
        
        # Process each building
        for idx, row in buildings_gdf.iterrows():
            # Get simplified geometry to improve performance
            geom = row.geometry.simplify(0.00001)
            
            # Skip invalid geometries
            if not geom.is_valid or geom.is_empty:
                continue
                
            # Check if building is flooded
            is_flooded = row[flooded_column] if flooded_column in buildings_gdf.columns else False
            
            # Set the color based on flood status
            color = flooded_color if is_flooded else non_flooded_color
            
            # For polygons, we'll rasterize them
            if geom.geom_type == 'Polygon':
                exterior = list(geom.exterior.coords)
                # Convert to pixel coordinates
                pixels = [world_to_pixel(x, y) for x, y in exterior]
                
                # Draw filled polygon (simple approach)
                from PIL import Image, ImageDraw
                img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
                draw = ImageDraw.Draw(img)
                draw.polygon(pixels, fill=(*color, opacity))
                
                # Merge with our rgba array
                temp_array = np.array(img)
                mask = temp_array[:, :, 3] > 0
                rgba[mask] = temp_array[mask]
            
            elif geom.geom_type == 'MultiPolygon':
                # Process each polygon in the multipolygon
                for poly in geom.geoms:
                    exterior = list(poly.exterior.coords)
                    pixels = [world_to_pixel(x, y) for x, y in exterior]
                    
                    # Draw filled polygon
                    from PIL import Image, ImageDraw
                    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))
                    draw = ImageDraw.Draw(img)
                    draw.polygon(pixels, fill=(*color, opacity))
                    
                    # Merge with our rgba array
                    temp_array = np.array(img)
                    mask = temp_array[:, :, 3] > 0
                    rgba[mask] = temp_array[mask]
        
        # Save the final image
        Image.fromarray(rgba).save(output_png)
        print(f"Successfully saved rasterized buildings to {output_png}")
        return output_png
        
    except Exception as e:
        print(f"Error rasterizing buildings: {e}")
        
        # Create a fallback image
        fallback = np.zeros((height, width, 4), dtype=np.uint8)
        fallback[10:20, :] = [255, 0, 0, 255]  # Red strip at top
        Image.fromarray(fallback).save(output_png)
        print(f"Saved fallback image to {output_png}")
        return output_png

def create_webmap(webmap_data):
    """
    Create an interactive webmap showing the flood protection scenarios
    
    Parameters:
    webmap_data (dict): Dictionary containing data prepared for the webmap
    
    Returns:
    str: Path to the generated webmap
    """
    # Extract data
    fema_web_path = webmap_data["fema_web_path"]
    fema_bounds = webmap_data["fema_bounds"]
    fema_trim_web_path = webmap_data["fema_trim_web_path"]
    fema_trim_bounds = webmap_data["fema_trim_bounds"]
    alignment = webmap_data["alignment"]
    buildings = webmap_data.get("buildings", None)
    
    # Print debugging information 
    print("\n==== WEBMAP CREATION: DATA CHECK ====")
    print(f"FEMA overlay path: {fema_web_path} (exists: {os.path.exists(fema_web_path)})")
    print(f"FEMA bounds: {fema_bounds}")
    print(f"Buildings data available: {buildings is not None}")
    print(f"Alignment data available: {alignment is not None}")
    
    # Force everything to EPSG:4326 (WGS84)
    if alignment is not None and str(alignment.crs) != "EPSG:4326":
        print(f"Converting alignment from {alignment.crs} to EPSG:4326")
        alignment = alignment.to_crs("EPSG:4326")
    
    if buildings is not None and str(buildings.crs) != "EPSG:4326":
        print(f"Converting buildings from {buildings.crs} to EPSG:4326")
        buildings = buildings.to_crs("EPSG:4326")
    
    # Find center and bounds for the map
    if buildings is not None and not buildings.empty:
        try:
            total_bounds = buildings.total_bounds
            center_lat = (total_bounds[1] + total_bounds[3]) / 2
            center_lng = (total_bounds[0] + total_bounds[2]) / 2
            map_bounds = [[total_bounds[1], total_bounds[0]], [total_bounds[3], total_bounds[2]]]
            print(f"Using buildings for map center: ({center_lng}, {center_lat})")
            print(f"Map bounds: {map_bounds}")
        except Exception as e:
            print(f"Error getting buildings bounds: {e}")
            center_lat, center_lng = 40.7, -74.0  # New York City
            map_bounds = None
    else:
        print("No buildings data, using default center")
        center_lat, center_lng = 40.7, -74.0  # New York City
        map_bounds = None
    
    # Create a basic map - using Leaflet directly
    m = folium.Map(
        location=[center_lat, center_lng],
        zoom_start=12,
        tiles='CartoDB positron'
    )
    
    # Create TWO SEPARATE MAPS instead of layer groups
    # This bypasses potential layer visibility issues
    
    # Add a custom HTML control for switching between maps
    html = """
    <div style="position: fixed; top: 10px; right: 10px; z-index: 1000; background: white; padding: 10px; border-radius: 5px;">
        <strong>Scenario: </strong>
        <select id="scenario-selector" onchange="switchScenario()">
            <option value="existing">Existing Conditions</option>
            <option value="protected">With Protection</option>
        </select>
    </div>
    
    <div id="map-existing" style="width: 100%; height: 100%; position: absolute; top: 0; left: 0; z-index: 500;"></div>
    <div id="map-protected" style="width: 100%; height: 100%; position: absolute; top: 0; left: 0; z-index: 400; display: none;"></div>
    
    <script>
        function switchScenario() {
            var scenario = document.getElementById('scenario-selector').value;
            if (scenario === 'existing') {
                document.getElementById('map-existing').style.display = 'block';
                document.getElementById('map-protected').style.display = 'none';
            } else {
                document.getElementById('map-existing').style.display = 'none';
                document.getElementById('map-protected').style.display = 'block';
            }
        }
    </script>
    """
    
    m.get_root().html.add_child(folium.Element(html))
    
    # Add buildings - EXISTING CONDITIONS
    if buildings is not None:
        print("Converting buildings to simplified GeoJSON...")
        
        # First, create a simplified version with only needed columns
        simple_buildings = buildings.copy()
        keep_columns = ['geometry']
        
        if 'flood_02pct' in buildings.columns:
            keep_columns.append('flood_02pct')
        if 'flood_02pct_scenario1' in buildings.columns:
            keep_columns.append('flood_02pct_scenario1')
            
        # Add any other columns for tooltips
        for col in ['name', 'bin', 'heightroof']:
            if col in buildings.columns:
                keep_columns.append(col)
                
        # Keep only necessary columns and convert booleans to integers for geojson
        for col in keep_columns:
            if col != 'geometry' and buildings[col].dtype == bool:
                simple_buildings[col] = simple_buildings[col].astype(int)
                
        # Direct GeoJSON creation with simplified styling
        print("Adding buildings to map...")
        
        # Custom style function as a JavaScript function string
        style_function_js = """
        function(feature) {
            var isFlooded = feature.properties.flood_02pct;
            if (isFlooded) {
                return {
                    fillColor: '#FF0000',
                    color: '#000000',
                    weight: 1,
                    fillOpacity: 0.8
                };
            } else {
                return {
                    fillColor: '#444444',
                    color: '#000000',
                    weight: 0.5,
                    fillOpacity: 0.6
                };
            }
        }
        """
        
        # Create custom GeoJSON script element
        geojson_data = simple_buildings.to_json()
        
        # Add existing conditions map
        buildings_script = f"""
        <script>
            var existingMap = L.map('map-existing').setView([{center_lat}, {center_lng}], 12);
            L.tileLayer('https://cartodb-basemaps-{{s}}.global.ssl.fastly.net/light_all/{{z}}/{{x}}/{{y}}.png', {{
                attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="http://cartodb.com/attributions">CartoDB</a>'
            }}).addTo(existingMap);
            
            var buildingsData = {geojson_data};
            
            L.geoJSON(buildingsData, {{
                style: {style_function_js}
            }}).addTo(existingMap);
            
            // Add FEMA overlay if exists
            """
        
        if os.path.exists(fema_web_path):
            buildings_script += f"""
            L.imageOverlay(
                '{fema_web_path}',
                [[{fema_bounds[1]}, {fema_bounds[0]}], [{fema_bounds[3]}, {fema_bounds[2]}]],
                {{opacity: 0.6}}
            ).addTo(existingMap);
            """
            
        # Add alignment if exists
        if alignment is not None:
            alignment_json = alignment.to_json()
            buildings_script += f"""
            var alignmentData = {alignment_json};
            L.geoJSON(alignmentData, {{
                style: function() {{
                    return {{
                        color: '#000000',
                        weight: 3,
                        opacity: 1.0
                    }};
                }}
            }}).addTo(existingMap);
            """
            
        # Add protected conditions map
        buildings_script += f"""
            // PROTECTED CONDITIONS MAP
            var protectedMap = L.map('map-protected').setView([{center_lat}, {center_lng}], 12);
            L.tileLayer('https://cartodb-basemaps-{{s}}.global.ssl.fastly.net/light_all/{{z}}/{{x}}/{{y}}.png', {{
                attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="http://cartodb.com/attributions">CartoDB</a>'
            }}).addTo(protectedMap);
            
            // Protected buildings
            L.geoJSON(buildingsData, {{
                style: function(feature) {{
                    var isFlooded = feature.properties.flood_02pct_scenario1;
                    if (isFlooded) {{
                        return {{
                            fillColor: '#FF0000',
                            color: '#000000',
                            weight: 1,
                            fillOpacity: 0.8
                        }};
                    }} else {{
                        return {{
                            fillColor: '#444444',
                            color: '#000000',
                            weight: 0.5,
                            fillOpacity: 0.6
                        }};
                    }}
                }}
            }}).addTo(protectedMap);
        """
        
        # Add trimmed FEMA overlay
        if os.path.exists(fema_trim_web_path):
            buildings_script += f"""
            L.imageOverlay(
                '{fema_trim_web_path}',
                [[{fema_trim_bounds[1]}, {fema_trim_bounds[0]}], [{fema_trim_bounds[3]}, {fema_trim_bounds[2]}]],
                {{opacity: 0.6}}
            ).addTo(protectedMap);
            """
            
        # Add alignment to protected view also
        if alignment is not None:
            buildings_script += f"""
            L.geoJSON(alignmentData, {{
                style: function() {{
                    return {{
                        color: '#000000',
                        weight: 3,
                        opacity: 1.0
                    }};
                }}
            }}).addTo(protectedMap);
            """
            
        # Set bounds for both maps
        if map_bounds:
            bounds_str = f"[[{map_bounds[0][0]}, {map_bounds[0][1]}], [{map_bounds[1][0]}, {map_bounds[1][1]}]]"
            buildings_script += f"""
            existingMap.fitBounds({bounds_str});
            protectedMap.fitBounds({bounds_str});
            """
            
        buildings_script += """
        </script>
        """
        
        m.get_root().html.add_child(folium.Element(buildings_script))
    
    # Add legend - simplified for clarity
    legend_html = f'''
    <div style="{LEGEND_STYLES['container']}">
        <h4 style="{LEGEND_STYLES['header']}">Legend</h4>
        
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: #FF0000;"></div>
            <span style="{LEGEND_STYLES['label']}">Flooded Buildings</span>
        </div>
        
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: #444444;"></div>
            <span style="{LEGEND_STYLES['label']}">Non-Flooded Buildings</span>
        </div>
        
        <div style="{LEGEND_STYLES['itemContainer']}">
            <div style="{LEGEND_STYLES['colorBox']}; background-color: #000000; height: 3px; margin-top: 8px;"></div>
            <span style="{LEGEND_STYLES['label']}">Flood Protection Alignment</span>
        </div>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(legend_html))
    
    # Add title
    title_html = f'''
    <div style="{TITLE_STYLE['container']}">
        <h1 style="{TITLE_STYLE['title']}">Flood Protection Scenarios</h1>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(title_html))
    
    # Add info box with statistics
    # Calculate the percent reduction
    percent_reduction = 0
    if count_existing > 0:
        percent_reduction = 100 - (count_protected / count_existing * 100)
    
    info_html = f'''
    <div style="{INFO_STYLE['container']}">
        <p style="{INFO_STYLE['text']}">
            {count_existing} buildings flooded in existing conditions<br>
            {count_protected} buildings remain flooded with protection<br>
            {percent_reduction:.1f}% reduction in flooded buildings
        </p>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(info_html))
    
    # Add logo
    logo_html = f'''
    <div style="{LOGO_STYLE['container']}" class="analysis-text">
        <span style="font-family: 'Futura', sans-serif; font-weight: bold; color: #4c5da4;">one</span>
        <span style="font-family: 'Futura', sans-serif; font-weight: 300; color: #4c5da4;"> analysis</span>
    </div>
    '''
    
    m.get_root().html.add_child(folium.Element(logo_html))
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(DATASET_INFO["Output"]["Webmap"]["path"])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save the map to an HTML file
    output_path = str(DATASET_INFO["Output"]["Webmap"]["path"])
    m.save(output_path)
    print(f"Web map saved to {output_path}")
    
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
    
    # Print the full analysis_results to help debug
    print("DEBUG: Analysis results keys:", list(analysis_results.keys()))
    for key, value in analysis_results.items():
        if isinstance(value, (int, float, str, bool)):
            print(f"DEBUG: {key} = {value}")
        elif hasattr(value, "shape"):  # For numpy arrays
            print(f"DEBUG: {key} shape = {value.shape}")
        elif hasattr(value, "__len__"):  # For lists, dicts, etc.
            print(f"DEBUG: {key} length = {len(value)}")
        else:
            print(f"DEBUG: {key} type = {type(value)}")
    
    # Ensure all required keys exist in the analysis_results dictionary
    webmap_data = {}
    
    # Required paths - check if they exist in analysis_results, otherwise use defaults
    if "fema_web_path" not in analysis_results:
        print("WARNING: fema_web_path not found in analysis results. Using processed FEMA data.")
        # Fix: Use the correct key "FEMA_Flood" instead of "FEMA_FloodHaz"
        fema_path = DATASET_INFO["Input"]["FEMA_Flood"]["path"]
        fema_web_dir = os.path.join(os.path.dirname(output_path), "web_data")
        os.makedirs(fema_web_dir, exist_ok=True)
        fema_web_path = os.path.join(fema_web_dir, "fema_flood.png")
        
        # Process the FEMA raster for web if it exists
        if os.path.exists(fema_path):
            print(f"Processing FEMA raster: {fema_path}")
            fema_web_path, fema_bounds = process_raster_for_web(fema_path, fema_web_path)
        else:
            print(f"ERROR: FEMA flood hazard file not found: {fema_path}")
            # Use fallback values
            fema_web_path = ""
            fema_bounds = (-74.26, 40.49, -73.69, 40.91)  # NYC fallback bounds
        
        webmap_data["fema_web_path"] = fema_web_path
        webmap_data["fema_bounds"] = fema_bounds
    else:
        webmap_data["fema_web_path"] = analysis_results["fema_web_path"]
        webmap_data["fema_bounds"] = analysis_results["fema_bounds"]
    
    # Same for trimmed FEMA data
    if "fema_trim_web_path" not in analysis_results:
        print("WARNING: fema_trim_web_path not found in analysis results. Using fallback.")
        webmap_data["fema_trim_web_path"] = webmap_data["fema_web_path"]  # Fallback to regular FEMA
        webmap_data["fema_trim_bounds"] = webmap_data["fema_bounds"]
    else:
        webmap_data["fema_trim_web_path"] = analysis_results["fema_trim_web_path"]
        webmap_data["fema_trim_bounds"] = analysis_results["fema_trim_bounds"]
    
    # Process other data
    if "alignment" in analysis_results:
        print(f"Found alignment data with type: {type(analysis_results['alignment'])}")
        webmap_data["alignment"] = analysis_results["alignment"]
    else:
        print("WARNING: No alignment data found. Looking for alignment file in config...")
        # Try to load alignment from file
        alignment_path = DATASET_INFO["Input"]["Alignment"]["path"]
        if os.path.exists(alignment_path):
            try:
                print(f"Loading alignment from: {alignment_path}")
                webmap_data["alignment"] = gpd.read_file(alignment_path)
            except Exception as e:
                print(f"Error loading alignment: {e}")
                webmap_data["alignment"] = None
        else:
            print(f"Alignment file not found: {alignment_path}")
            webmap_data["alignment"] = None
    
    if "buildings" in analysis_results:
        print(f"Found buildings data with type: {type(analysis_results['buildings'])}")
        webmap_data["buildings"] = analysis_results["buildings"]
    else:
        print("WARNING: No buildings data found. Looking for buildings file in config...")
        # Try to load buildings from output file
        buildings_path = DATASET_INFO["Output"]["Buildings_scenarios"]["path"]
        if os.path.exists(buildings_path):
            try:
                print(f"Loading buildings from: {buildings_path}")
                buildings_gdf = gpd.read_file(buildings_path)
                webmap_data["buildings"] = buildings_gdf
            except Exception as e:
                print(f"Error loading buildings: {e}")
                webmap_data["buildings"] = None
        else:
            print(f"Buildings file not found: {buildings_path}")
            webmap_data["buildings"] = None

    # Add these to globals so create_webmap can access them
    global count_existing, count_protected

    if webmap_data["buildings"] is not None:
        buildings = webmap_data["buildings"]
        if 'flood_02pct' in buildings.columns:
            count_existing = int(buildings['flood_02pct'].sum())
            print(f"Calculated count_existing: {count_existing}")
        else:
            print("WARNING: 'flood_02pct' column not found in buildings data")
        
        if 'flood_02pct_scenario1' in buildings.columns:
            count_protected = int(buildings['flood_02pct_scenario1'].sum())
            print(f"Calculated count_protected: {count_protected}")
        else:
            print("WARNING: 'flood_02pct_scenario1' column not found in buildings data")
    else:
        # Use values from analysis_results if available
        count_existing = analysis_results.get("count_existing", 0)
        count_protected = analysis_results.get("count_protected", 0)
    
    print(f"Using count values: existing={count_existing}, protected={count_protected}")
    
    return create_webmap(webmap_data)


if __name__ == "__main__":
    # This is only for directly running webmap.py
    # Import here to avoid circular imports when called from main.py
    import sys
    # Add the project directory to the path if needed
    import os
    project_dir = os.path.dirname(os.path.abspath(__file__))
    if project_dir not in sys.path:
        sys.path.append(project_dir)
        
    from preprocessing import preprocess_data
    from analysis import run_analysis
    
    preprocessed_data = preprocess_data()
    analysis_results = run_analysis(preprocessed_data)
    webmap_path = generate_webmap(analysis_results)
    
    # Open the webmap in a browser
    webbrowser.open('file://' + os.path.abspath(webmap_path))