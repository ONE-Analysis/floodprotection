"""
Configuration file for Flood Infrastructure Analysis
Contains paths, settings, and parameters used throughout the analysis
"""
import os
from pathlib import Path

# Project structure and paths
ROOT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = ROOT_DIR / "input"
OUTPUT_DIR = ROOT_DIR / "output"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# CRS settings
PROJECT_CRS = "EPSG:6539"  # New York State Plane Long Island zone, meters
WEB_CRS = "EPSG:4326"      # WGS 84, for web mapping

# Processing toggles
OVERWRITE = {
    "site_bounds": True,
    "alignment": True,
    "vector_crop": True,
    "raster_crop": True,
    "fema_flood_crop": True,
    "buildings_flood": True,
    "buildings_scenarios": True,
    "webmap": True
}

# Map settings
RESOLUTION = 10  # Resolution in feet for web raster processing
INITIAL_ZOOM = 16  # Initial zoom level for the webmap

# Data dictionary for input and output datasets
DATASET_INFO = {
    # Input datasets
    "Input": {
        "Site_Bounds": {
            "path": INPUT_DIR / "site_boundary.geojson",
            "description": "Site boundary polygon"
        },
        "Alignment": {
            "path": INPUT_DIR / "alignment.geojson",
            "description": "Flood protection alignment"
        },
        "Buildings": {
            "path": INPUT_DIR / "buildings.geojson",
            "description": "Building footprints"
        },
        "NSI": {
            "path": INPUT_DIR / "nsi.geojson",
            "description": "National Structure Inventory points"
        },
        "Water": {
            "path": INPUT_DIR / "water.geojson",
            "description": "Water bodies"
        },
        "FEMA_Flood": {
            "path": INPUT_DIR / "fema_flood.tif",
            "description": "FEMA flood zones"
        },
        "DEM": {
            "path": INPUT_DIR / "dem.tif",
            "description": "Digital Elevation Model"
        }
    },
    
    # Preprocessed datasets
    "Preprocessed": {
        "Buildings_crop": {
            "path": INPUT_DIR / "NYC_Buildings_crop.geojson",
            "description": "Cropped NYC building footprints"
        },
        "NSI_crop": {
            "path": INPUT_DIR / "NSI_crop.geojson",
            "description": "Cropped National Structure Inventory data"
        },
        "Water_crop": {
            "path": INPUT_DIR / "Water_crop.geojson",
            "description": "Cropped water polygon for determining flood direction"
        },
        "FEMA_Flood_crop": {
            "path": INPUT_DIR / "FEMA_Flood_crop.tif",
            "description": "Clipped FEMA flood zones"
        }
    },
    
    # Output datasets with unique keys
    "Output": {
        "Buildings_Scenarios": {
            "path": OUTPUT_DIR / "buildings_scenarios.geojson",
            "description": "Buildings with flood exposure data for all scenarios"
        },
        "FEMA_Flood_crop": {
            "path": OUTPUT_DIR / "FEMA_flood_crop.tif",
            "description": "Clipped FEMA flood zones for existing scenario"
        },
        "FEMA_Flood_crop_trim": {
            "path": OUTPUT_DIR / "intermediate/fema_flood_crop_trim.tif",
            "description": "Trimmed FEMA flood raster (one side of alignment)"
        },
        "FEMA_Flood_crop_web": {
            "path": OUTPUT_DIR / "web/fema_flood.png",
            "description": "FEMA flood raster for web (PNG)"
        },
        "FEMA_Flood_crop_trim_web": {
            "path": OUTPUT_DIR / "web/fema_flood_trim.png",
            "description": "Trimmed FEMA flood raster for web (PNG)"
        },
        "Webmap": {
            "path": OUTPUT_DIR / "Flood_Protection.html",
            "description": "Interactive web map showing flood protection scenarios"
        },
        "DEM_crop": {
            "path": OUTPUT_DIR / "intermediate/dem_crop.tif",
            "description": "Cropped DEM raster"
        },
        "Buildings_flood": {
            "path": OUTPUT_DIR / "buildings_flood.geojson",
            "description": "Buildings with flood exposure"
        },
        "Buildings_scenarios": {
            "path": OUTPUT_DIR / "buildings_scenarios.geojson",
            "description": "Buildings with flood protection scenarios"
        }
    },
        
        # Webmap visualization settings
        "Webmap": {
            "FEMA_FloodHaz": {
                "hex_1pct": "#0000FF",       # Blue for 1% annual chance flood
                "hex_0_2pct": "#96C4FF",     # Light blue for 0.2% annual chance flood
                "opacity": 0.6
            },
            "Buildings": {
                "flooded": "#FF0000",        # Red for flooded buildings
                "not_flooded": "#444444",    # Dark grey for non-flooded buildings
                "opacity": 0.8
            }
        }
    }

# Legend styles for webmap
LEGEND_STYLES = {
    "container": "position: fixed; bottom: 50px; right: 50px; z-index:9999; background: white; padding: 15px; border-radius: 15px; font-family: Arial, sans-serif; box-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);",
    "header": "margin-top:0; margin-bottom: 12px; font-size: 16px; font-weight: bold;",
    "sectionHeader": "margin-top:10px; margin-bottom: 5px; font-size: 14px; font-weight: bold;",
    "itemContainer": "display: flex; align-items: center; margin-bottom: 5px;",
    "colorBox": "width: 20px; height: 20px; margin-right: 5px;",
    "label": "font-size: 13px;"
}

# Title box style
TITLE_STYLE = {
    "container": "position: fixed; top: 10px; left: 10px; z-index:9999; background: white; padding: 10px 20px; border-radius: 10px; font-family: Arial, sans-serif; box-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);",
    "title": "margin: 0; font-size: 20px; font-weight: bold; color: #333;"
}

# Info box style
INFO_STYLE = {
    "container": "position: fixed; top: 10px; right: 10px; z-index:9999; background: white; padding: 10px 15px; border-radius: 10px; font-family: Arial, sans-serif; box-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);",
    "text": "margin: 0; font-size: 14px; font-weight: bold; color: #333;"
}

# Logo style
LOGO_STYLE = {
    "container": "position: fixed; top: 100px; left: 10px; z-index:9999; background: white; padding: 5px 10px; border-radius: 5px; font-family: 'Futura', Arial, sans-serif; box-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);"
}

# Helper functions for working with colors
def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple"""
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 8:  # With alpha
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    else:
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def hex_to_rgba(hex_color):
    """Convert hex color with alpha to RGBA tuple"""
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 8:  # With alpha
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4, 6))
    else:
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4)) + (255,)

def interpolate_color(t, start_hex, end_hex):
    """Interpolate between two hex colors"""
    # Convert hex to RGB
    r_1, g_1, b_1 = hex_to_rgb(start_hex)
    r_2, g_2, b_2 = hex_to_rgb(end_hex)
    
    # Interpolate
    r = int(r_1 + (r_2 - r_1) * t)
    g = int(g_1 + (g_2 - g_1) * t)
    b = int(b_1 + (b_2 - b_1) * t)
    
    # Convert back to hex
    return f"#{r:02x}{g:02x}{b:02x}"

def interpolate_color_with_alpha(t, start_hex, end_hex):
    """Interpolate between two hex colors with alpha values"""
    # Convert hex to RGBA
    r_1, g_1, b_1, a_1 = hex_to_rgba(start_hex)
    r_2, g_2, b_2, a_2 = hex_to_rgba(end_hex)
    
    # Interpolate
    r = int(r_1 + (r_2 - r_1) * t)
    g = int(g_1 + (g_2 - g_1) * t)
    b = int(b_1 + (b_2 - b_1) * t)
    a = int(a_1 + (a_2 - a_1) * t)
    
    # Convert back to hex with alpha
    return f"#{r:02x}{g:02x}{b:02x}{a:02x}"

# Large polygon handling
LARGE_POLYGON_THRESHOLD = 10 * 4046.86  # 30 acres in square meters
POINTS_PER_AREA = 1 / (5 * 4046.86)    # 1 point per 15 acres