"""
Flood Protection Analysis - Side-by-Side Map Visualization
Creates two maps comparing current flood conditions with protection scenario
"""

import os
import geopandas as gpd
import numpy as np
import rasterio
from rasterio import warp
from PIL import Image
from pathlib import Path
import webbrowser
import json

# Constants and styling
WEB_CRS = "EPSG:4326"
OUTPUT_PATH = "./output/multi_flood_scenario_map.html"

# Color definitions - centralized for consistency
FLOOD_COLORS = {
    "annual_0_2_percent": "#75E1FF",  # 0.2% annual chance (raster value 1)
    "annual_1_percent": "#BEE7FF",    # 1% annual chance (raster value 2)
    "buildings_flooded": "#0040A0",   # Flooded buildings
    "alignment": "#000000",           # Flood protection alignment
    "boundary": "#555555"             # Area boundary
}

# Convert hex to RGB for raster processing (with 70% opacity)
def hex_to_rgba(hex_color, alpha=179):  # 179 is ~70% of 255
    h = hex_color.lstrip('#')
    return [int(h[i:i+2], 16) for i in (0, 2, 4)] + [alpha]

# Style dictionaries
MAP_STYLES = {
    "width": "45%",
    "height": "500px",
    "display": "inline-block",
    "border": "1px solid #ccc",
    "border-radius": "5px",
    "margin": "10px 2%",  # Increased margin for better spacing
    "vertical-align": "top"  # Ensure consistent vertical alignment
}

TITLE_STYLE = {
    "container": "background-color: rgba(255,255,255,0.8); padding: 5px 10px; border-radius: 4px; text-align: left;",
    "title": "margin: 0; font-size: 16px; font-weight: bold; font-family: Arial, sans-serif;"
}

LEGEND_STYLE = {
    "container": "background-color: white; padding: 10px; border-radius: 5px; box-shadow: 0 0 5px rgba(0,0,0,0.5);",
    "header": "margin: 0; font-size: 14px; font-weight: bold;",
    "item": "margin: 5px 0; font-size: 12px;",
    "color_box": "display: inline-block; width: 15px; height: 15px; margin-right: 5px; vertical-align: middle;"
}

STATS_BOX_STYLE = {
    "container": "background-color: white; padding: 10px; border-radius: 5px; box-shadow: 0 0 5px rgba(0,0,0,0.5);",
    "header": "margin: 0; font-size: 14px; font-weight: bold;",
    "text": "font-size: 12px; margin: 0; line-height: 1.4;"
}

# Style definitions for vector features
BUILDING_STYLE = {
    "flooded": {
        "fillColor": FLOOD_COLORS["buildings_flooded"],
        "color": "#000000",
        "weight": 0.5,
        "fillOpacity": 0.7
    }
}

ALIGNMENT_STYLE = {
    "color": FLOOD_COLORS["alignment"],
    "weight": 2,
    "opacity": 1
}

BOUNDARY_STYLE = {
    "color": FLOOD_COLORS["boundary"],
    "weight": 2,
    "opacity": 1,
    "fillOpacity": 0,
    "dashArray": "5, 5"  # Dashed line
}

def create_maps(buildings_path, fema_raster_path, fema_trim_raster_path, alignment_path, boundary_path):
    """
    Create side-by-side maps comparing flood scenarios
    
    Parameters:
    buildings_path (str): Path to buildings GeoJSON with flood attributes
    fema_raster_path (str): Path to FEMA flood raster
    fema_trim_raster_path (str): Path to trimmed FEMA flood raster
    alignment_path (str): Path to flood protection alignment GeoJSON
    boundary_path (str): Path to site boundary GeoJSON
    
    Returns:
    str: Path to the generated HTML file
    """
    print("\n=== Loading Data ===")
    
    # Load buildings data
    buildings_flooded_existing = None
    buildings_flooded_protected = None
    buildings = None
    try:
        print(f"Loading buildings from: {buildings_path}")
        buildings = gpd.read_file(buildings_path)
        if buildings.crs != WEB_CRS:
            buildings = buildings.to_crs(WEB_CRS)
        
        # Check column existence
        required_columns = ['flood_02pct', 'flood_02pct_scenario1']
        missing_columns = [col for col in required_columns if col not in buildings.columns]
        if missing_columns:
            print(f"Warning: Missing columns in buildings data: {missing_columns}")
            # Create default columns if missing
            for col in missing_columns:
                buildings[col] = False
        
        # Filter for flooded buildings only
        buildings_flooded_existing = buildings[buildings['flood_02pct']].copy()
        buildings_flooded_protected = buildings[buildings['flood_02pct_scenario1']].copy()
        
        print(f"Total buildings: {len(buildings)}")
        print(f"Flooded buildings (existing): {len(buildings_flooded_existing)}")
        print(f"Flooded buildings (protected): {len(buildings_flooded_protected)}")
    except Exception as e:
        print(f"Error loading buildings data: {e}")
        buildings = None
    
    # Load site boundary
    boundary = None
    try:
        print(f"Loading boundary from: {boundary_path}")
        boundary = gpd.read_file(boundary_path)
        if boundary.crs != WEB_CRS:
            boundary = boundary.to_crs(WEB_CRS)
        print(f"Boundary loaded with {len(boundary)} features")
    except Exception as e:
        print(f"Error loading site boundary: {e}")
        boundary = None
    
    # Load alignment
    alignment = None
    try:
        print(f"Loading alignment from: {alignment_path}")
        alignment = gpd.read_file(alignment_path)
        if alignment.crs != WEB_CRS:
            alignment = alignment.to_crs(WEB_CRS)
        print(f"Alignment loaded with {len(alignment)} features")
    except Exception as e:
        print(f"Error loading alignment data: {e}")
        alignment = None
    
    # Count flooded buildings
    count_existing = 0
    count_protected = 0
    buildings_protected = 0
    percent_reduction = 0.0
    
    if buildings is not None:
        count_existing = int(buildings['flood_02pct'].sum())
        count_protected = int(buildings['flood_02pct_scenario1'].sum())
        buildings_protected = count_existing - count_protected
        if count_existing > 0:
            percent_reduction = (buildings_protected / count_existing) * 100
    
    # Get map center and bounds - use the buildings data if available
    if buildings is not None and not buildings.empty:
        # Get the bounds of all buildings
        bounds = buildings.total_bounds  # [minx, miny, maxx, maxy]
        # Convert to Leaflet format [[miny, minx], [maxy, maxx]] = [[south, west], [north, east]]
        map_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
        # Use centroid for initial view
        center = buildings.unary_union.centroid
        center = (center.y, center.x)  # lat, lon
        print(f"Using buildings data for map center: {center}")
    elif boundary is not None and not boundary.empty:
        # If no buildings, use boundary
        bounds = boundary.total_bounds
        map_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
        center = boundary.unary_union.centroid
        center = (center.y, center.x)
        print(f"Using boundary data for map center: {center}")
    else:
        # Default to NYC
        center = (40.7128, -74.0060)
        map_bounds = [[40.6, -74.1], [40.8, -73.9]]  # Rough NYC bounds
        print("No buildings or boundary data available, using default NYC center and bounds")
    
    print("\n=== Processing Rasters ===")
    
    # Process FEMA raster to PNG
    web_data_dir = os.path.join(os.path.dirname(OUTPUT_PATH), "web_data")
    os.makedirs(web_data_dir, exist_ok=True)
    
    # Process original FEMA flood raster
    fema_png = os.path.join(web_data_dir, "fema_flood.png")
    fema_bounds = None
    
    try:
        print(f"Processing FEMA flood raster: {fema_raster_path}")
        with rasterio.open(fema_raster_path) as src:
            # Read data
            data = src.read(1)
            raster_bounds = src.bounds
            
            # Convert to WGS84 bounds if needed
            if src.crs != WEB_CRS:
                fema_bounds = warp.transform_bounds(src.crs, WEB_CRS, *raster_bounds)
            else:
                fema_bounds = raster_bounds
            
            # Create RGBA image
            height, width = data.shape
            rgba = np.zeros((height, width, 4), dtype=np.uint8)
            
            # Apply colors from our central color definitions
            rgba[data == 1] = hex_to_rgba(FLOOD_COLORS["annual_0_2_percent"])
            rgba[data == 2] = hex_to_rgba(FLOOD_COLORS["annual_1_percent"])
            
            # Save as PNG
            Image.fromarray(rgba).save(fema_png)
            print(f"Processed FEMA flood raster to {fema_png}")
            print(f"FEMA bounds: {fema_bounds}")
    except Exception as e:
        print(f"Error processing FEMA raster: {e}")
        fema_png = None
        fema_bounds = None
    
    # Process trimmed FEMA flood raster
    fema_trim_png = os.path.join(web_data_dir, "fema_flood_trim.png")
    fema_trim_bounds = None
    
    try:
        print(f"Processing trimmed FEMA flood raster: {fema_trim_raster_path}")
        with rasterio.open(fema_trim_raster_path) as src:
            # Read data
            data = src.read(1)
            raster_bounds = src.bounds
            
            # Convert to WGS84 bounds if needed
            if src.crs != WEB_CRS:
                fema_trim_bounds = warp.transform_bounds(src.crs, WEB_CRS, *raster_bounds)
            else:
                fema_trim_bounds = raster_bounds
            
            # Create RGBA image
            height, width = data.shape
            rgba = np.zeros((height, width, 4), dtype=np.uint8)
            
            # Apply colors from our central color definitions
            rgba[data == 1] = hex_to_rgba(FLOOD_COLORS["annual_0_2_percent"])
            rgba[data == 2] = hex_to_rgba(FLOOD_COLORS["annual_1_percent"])
            
            # Save as PNG
            Image.fromarray(rgba).save(fema_trim_png)
            print(f"Processed trimmed FEMA flood raster to {fema_trim_png}")
            print(f"Trimmed FEMA bounds: {fema_trim_bounds}")
    except Exception as e:
        print(f"Error processing trimmed FEMA raster: {e}")
        fema_trim_png = None
        fema_trim_bounds = None
    
    print("\n=== Creating HTML ===")
    
    # Create HTML content
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Flood Protection Scenarios Comparison</title>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <link rel="stylesheet" href="https://unpkg.com/leaflet@1.7.1/dist/leaflet.css" />
        <script src="https://unpkg.com/leaflet@1.7.1/dist/leaflet.js"></script>
        <style>
            body {{ 
                margin: 0;
                padding: 20px;
                font-family: Arial, sans-serif;
                background-color: #f5f5f5;
            }}
            h1 {{
                font-family: Arial, sans-serif;
                font-size: 24px;
                margin-bottom: 20px;
                text-align: left;
                color: #333;
            }}
            .map-container {{
                display: flex;
                flex-wrap: wrap;
                justify-content: center;
                width: 100%;
            }}
            .map {{
                {'; '.join([f"{k}: {v}" for k, v in MAP_STYLES.items()])};
            }}
        </style>
    </head>
    <body>
        <h1>Flood Protection Scenarios Comparison</h1>
        <div class="map-container">
            <div id="map1" class="map"></div>
            <div id="map2" class="map"></div>
        </div>
        
        <script>
        // Initialize the maps
        var map1 = L.map('map1', {{
            scrollWheelZoom: true,
            dragging: true,
            doubleClickZoom: true
        }}).setView([{center[0]}, {center[1]}], 12.5);
        
        var map2 = L.map('map2', {{
            scrollWheelZoom: true,
            dragging: true,
            doubleClickZoom: true
        }}).setView([{center[0]}, {center[1]}], 12.5);
        
        // Remove zoom controls
        map1.removeControl(map1.zoomControl);
        map2.removeControl(map2.zoomControl);
        
        // Add basemap tiles
        L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png', {{
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
            maxZoom: 19
        }}).addTo(map1);
        
        L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png', {{
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>',
            maxZoom: 19
        }}).addTo(map2);
        
        // Add titles
        var title1 = L.control({{position: 'topleft'}});
        title1.onAdd = function() {{
            var div = L.DomUtil.create('div', 'info title');
            div.innerHTML = '<div style="{TITLE_STYLE["container"]}"><h2 style="{TITLE_STYLE["title"]}">Current Flood Hazard</h2></div>';
            return div;
        }};
        title1.addTo(map1);
        
        var title2 = L.control({{position: 'topleft'}});
        title2.onAdd = function() {{
            var div = L.DomUtil.create('div', 'info title');
            div.innerHTML = '<div style="{TITLE_STYLE["container"]}"><h2 style="{TITLE_STYLE["title"]}">With Flood Protection</h2></div>';
            return div;
        }};
        title2.addTo(map2);
    """
    
    # Add FEMA flood rasters if available
    if fema_png and fema_bounds:
        rel_path = os.path.relpath(fema_png, os.path.dirname(OUTPUT_PATH))
        html += f"""
                // Add FEMA flood raster to Map 1
                L.imageOverlay(
                    '{rel_path.replace(os.sep, "/")}',
                    [[{fema_bounds[1]}, {fema_bounds[0]}], [{fema_bounds[3]}, {fema_bounds[2]}]],
                    {{opacity: 0.7}}
                ).addTo(map1);
        """
    
    if fema_trim_png and fema_trim_bounds:
        rel_path = os.path.relpath(fema_trim_png, os.path.dirname(OUTPUT_PATH))
        html += f"""
                // Add trimmed FEMA flood raster to Map 2
                L.imageOverlay(
                    '{rel_path.replace(os.sep, "/")}',
                    [[{fema_trim_bounds[1]}, {fema_trim_bounds[0]}], [{fema_trim_bounds[3]}, {fema_trim_bounds[2]}]],
                    {{opacity: 0.7}}
                ).addTo(map2);
        """
    
    # Add boundary if available - ensure valid GeoJSON
    if boundary is not None:
        try:
            # Test serialization to ensure valid GeoJSON
            boundary_json = boundary.to_json()
            json.loads(boundary_json)  # Test if valid JSON
            html += f"""
                // Add boundary to both maps
                var boundaryStyle = {BOUNDARY_STYLE};
                var boundaryData = {boundary_json};
                
                L.geoJSON(boundaryData, {{
                    style: function() {{ return boundaryStyle; }}
                }}).addTo(map1);
                
                L.geoJSON(boundaryData, {{
                    style: function() {{ return boundaryStyle; }}
                }}).addTo(map2);
            """
        except Exception as e:
            print(f"Error serializing boundary to GeoJSON: {e}")
    
    # Add alignment if available - ensure valid GeoJSON
    if alignment is not None:
        try:
            alignment_json = alignment.to_json()
            json.loads(alignment_json)  # Test if valid JSON
            html += f"""
                // Add alignment to Map 2
                var alignmentStyle = {ALIGNMENT_STYLE};
                var alignmentData = {alignment_json};
                
                L.geoJSON(alignmentData, {{
                    style: function() {{ return alignmentStyle; }}
                }}).addTo(map2);
            """
        except Exception as e:
            print(f"Error serializing alignment to GeoJSON: {e}")
    
    # Add flooded buildings to maps - ensure valid GeoJSON
    if buildings_flooded_existing is not None and not buildings_flooded_existing.empty:
        try:
            # Convert to GeoJSON and test
            buildings_existing_json = buildings_flooded_existing.to_json()
            json.loads(buildings_existing_json)  # Test if valid JSON
            
            html += f"""
                // Add flooded buildings to Map 1 (existing scenario)
                var buildingStyle = {BUILDING_STYLE["flooded"]};
                var floodedBuildingsData = {buildings_existing_json};
                
                L.geoJSON(floodedBuildingsData, {{
                    style: function() {{ return buildingStyle; }}
                }}).addTo(map1);
            """
        except Exception as e:
            print(f"Error serializing existing flooded buildings to GeoJSON: {e}")
    
    if buildings_flooded_protected is not None and not buildings_flooded_protected.empty:
        try:
            # Convert to GeoJSON and test
            buildings_protected_json = buildings_flooded_protected.to_json()
            json.loads(buildings_protected_json)  # Test if valid JSON
            
            html += f"""
                // Add flooded buildings to Map 2 (protected scenario)
                var buildingStyle = {BUILDING_STYLE["flooded"]};
                var floodedBuildingsProtectedData = {buildings_protected_json};
                
                L.geoJSON(floodedBuildingsProtectedData, {{
                    style: function() {{ return buildingStyle; }}
                }}).addTo(map2);
            """
        except Exception as e:
            print(f"Error serializing protected flooded buildings to GeoJSON: {e}")
    
    # Add legends and finish HTML
    html += f"""
                // Add legend to Map 1
                var legend1 = L.control({{position: 'bottomright'}});
                legend1.onAdd = function() {{
                    var div = L.DomUtil.create('div', 'info legend');
                    div.innerHTML = `
                        <div style="{LEGEND_STYLE["container"]}">
                            <h4 style="{LEGEND_STYLE["header"]}">Legend</h4>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; border: 2px dashed {FLOOD_COLORS["boundary"]}; background-color: transparent;"></span>
                                Area of Interest
                            </div>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["buildings_flooded"]};"></span>
                                Flooded Buildings
                            </div>
                            <h5 style="margin: 5px 0 2px 0; font-size: 12px;">FEMA Floodmap</h5>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["annual_0_2_percent"]};"></span>
                                0.2% Annual Chance
                            </div>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["annual_1_percent"]};"></span>
                                1% Annual Chance
                            </div>
                        </div>
                    `;
                    return div;
                }};
                legend1.addTo(map1);
                
                // Add legend to Map 2
                var legend2 = L.control({{position: 'bottomright'}});
                legend2.onAdd = function() {{
                    var div = L.DomUtil.create('div', 'info legend');
                    div.innerHTML = `
                        <div style="{LEGEND_STYLE["container"]}">
                            <h4 style="{LEGEND_STYLE["header"]}">Legend</h4>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; border: 2px dashed {FLOOD_COLORS["boundary"]}; background-color: transparent;"></span>
                                Area of Interest
                            </div>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; border: 2px solid {FLOOD_COLORS["alignment"]};"></span>
                                Flood Protection
                            </div>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["buildings_flooded"]};"></span>
                                Flooded Buildings
                            </div>
                            <h5 style="margin: 5px 0 2px 0; font-size: 12px;">FEMA Floodmap</h5>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["annual_0_2_percent"]};"></span>
                                0.2% Annual Chance
                            </div>
                            <div style="{LEGEND_STYLE["item"]}">
                                <span style="{LEGEND_STYLE["color_box"]}; background-color: {FLOOD_COLORS["annual_1_percent"]};"></span>
                                1% Annual Chance
                            </div>
                        </div>
                    `;
                    return div;
                }};
                legend2.addTo(map2);
                // Add stats box to Map 2
                var statsBox = L.control({{position: 'topright'}});
                statsBox.onAdd = function() {{
                    var div = L.DomUtil.create('div', 'info stats');
                    div.innerHTML = '<div style="{STATS_BOX_STYLE["container"]}">' +
                        '<h4 style="{STATS_BOX_STYLE["header"]}">Protection Statistics</h4>' +
                        '<p style="{STATS_BOX_STYLE["text"]}">' +
                        '{count_existing} buildings flooded without protection<br>' +
                        '{count_protected} buildings flooded with protection<br>' +
                        '<b>{buildings_protected} buildings protected</b><br>' +
                        '<b>{percent_reduction:.1f}% reduction</b>' +
                        '</p>' +
                        '</div>';
                    return div;
                }};
                statsBox.addTo(map2);
                // Sync the maps
                map1.on('move', function() {{
                    map2.setView(map1.getCenter(), map1.getZoom(), {{animate: false}});
                }});
                
                map2.on('move', function() {{
                    map1.setView(map2.getCenter(), map2.getZoom(), {{animate: false}});
                }});
                
                // Fit to buildings data bounds
                var mapBounds = {map_bounds};
                map1.fitBounds(mapBounds);
                
                // Add a slight delay to ensure maps are properly loaded before fitting bounds
                setTimeout(function() {{
                    map1.invalidateSize();
                    map2.invalidateSize();
                    map1.fitBounds(mapBounds);
                }}, 500);
        </script>
    </body>
    </html>
    """
    
    # Write the HTML to file
    with open(OUTPUT_PATH, 'w') as f:
        f.write(html)
    
    print(f"Side-by-side map created at: {OUTPUT_PATH}")
    return OUTPUT_PATH

def main():
    # Set file paths
    buildings_path = "./output/buildings_scenarios.geojson"
    fema_raster_path = "./output/FEMA_flood_crop.tif"
    fema_trim_raster_path = "./output/intermediate/fema_flood_crop_trim.tif"
    alignment_path = "./output/intermediate/alignment.geojson"
    boundary_path = "./input/site_boundary.geojson"
    
    # Check if files exist
    for path in [buildings_path, fema_raster_path, fema_trim_raster_path, alignment_path, boundary_path]:
        if not os.path.exists(path):
            print(f"Warning: File not found: {path}")
    
    # Create the maps
    output_path = create_maps(buildings_path, fema_raster_path, fema_trim_raster_path, alignment_path, boundary_path)
    print(f"Open {output_path} in a web browser to view the maps")
    
    # Automatically open the HTML file in the default browser
    try:
        absolute_path = os.path.abspath(output_path)
        print(f"Opening {absolute_path} in web browser...")
        webbrowser.open('file://' + absolute_path)
    except Exception as e:
        print(f"Could not automatically open the map: {e}")
        print(f"Please manually open {output_path} in your browser")

if __name__ == "__main__":
    main()