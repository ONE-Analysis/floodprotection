"""
Main execution script for Flood Infrastructure Analysis
Orchestrates the entire workflow: preprocessing, analysis, and webmap generation
"""
import os
import sys
import time
import argparse
import webbrowser
from pathlib import Path

# Add the project directory to the path
project_dir = os.path.dirname(os.path.abspath(__file__))
if project_dir not in sys.path:
    sys.path.append(project_dir)

# Import project modules
from config import DATASET_INFO, OVERWRITE
from preprocessing import preprocess_data
from analysis import run_analysis, prepare_data_for_webmap
from webmap import create_webmap

def parse_arguments():
    """Parse command line arguments to determine which steps to run"""
    parser = argparse.ArgumentParser(description='Run flood infrastructure analysis with toggleable steps')
    parser.add_argument('--skip-preprocess', action='store_true', help='Skip preprocessing step')
    parser.add_argument('--skip-analysis', action='store_true', help='Skip analysis step')
    parser.add_argument('--skip-webmap', action='store_true', help='Skip webmap generation step')
    parser.add_argument('--open-map', action='store_true', help='Open the webmap in browser when done')
    return parser.parse_args()

def load_saved_preprocessed_data():
    """Load preprocessed data from saved files"""
    print("Loading preprocessed data from saved files...")
    preprocessed_data = {}
    
    # Load site bounds
    site_bounds_path = DATASET_INFO["Output"]["Site_Bounds"]["path"]
    if os.path.exists(site_bounds_path):
        import geopandas as gpd
        preprocessed_data["site_bounds"] = gpd.read_file(site_bounds_path)
        print(f"Loaded site bounds from {site_bounds_path}")
    else:
        print(f"Warning: Site bounds file not found at {site_bounds_path}")
    
    # Load alignment
    alignment_path = DATASET_INFO["Output"]["Alignment"]["path"]
    if os.path.exists(alignment_path):
        import geopandas as gpd
        preprocessed_data["alignment"] = gpd.read_file(alignment_path)
        print(f"Loaded alignment from {alignment_path}")
    else:
        print(f"Warning: Alignment file not found at {alignment_path}")
        
    return preprocessed_data

def load_saved_analysis_results():
    """Load analysis results from saved files"""
    print("Loading analysis results from saved files...")
    analysis_results = {}
    
    # Load buildings with flood data
    buildings_flood_path = DATASET_INFO["Output"]["Buildings_flood"]["path"]
    if os.path.exists(buildings_flood_path):
        import geopandas as gpd
        analysis_results["buildings_flood"] = gpd.read_file(buildings_flood_path)
        print(f"Loaded buildings flood data from {buildings_flood_path}")
    else:
        print(f"Warning: Buildings flood data not found at {buildings_flood_path}")
    
    # Load buildings with protection scenario data
    buildings_scenarios_path = DATASET_INFO["Output"]["Buildings_scenarios"]["path"]
    if os.path.exists(buildings_scenarios_path):
        import geopandas as gpd
        analysis_results["buildings_protected"] = gpd.read_file(buildings_scenarios_path)
        print(f"Loaded buildings protection scenario from {buildings_scenarios_path}")
    else:
        print(f"Warning: Buildings protection scenario not found at {buildings_scenarios_path}")
    
    # Load alignment (if not already loaded in preprocessing)
    alignment_path = DATASET_INFO["Output"]["Alignment"]["path"]
    if "alignment" not in analysis_results and os.path.exists(alignment_path):
        import geopandas as gpd
        analysis_results["alignment"] = gpd.read_file(alignment_path)
        print(f"Loaded alignment from {alignment_path}")
    
    return analysis_results

def main():
    """Main execution function"""
    start_time = time.time()
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directories
    for dataset_type in DATASET_INFO:
        if dataset_type == "Input":
            continue
        for dataset in DATASET_INFO[dataset_type]:
            if "path" in DATASET_INFO[dataset_type][dataset]:
                path = DATASET_INFO[dataset_type][dataset]["path"]
                os.makedirs(os.path.dirname(str(path)), exist_ok=True)
    
    # Step 1: Preprocess data (if not skipped)
    if not args.skip_preprocess:
        print("\n=== STEP 1: PREPROCESSING ===")
        preprocessed_data = preprocess_data()
    else:
        print("\n=== SKIPPING PREPROCESSING ===")
        preprocessed_data = load_saved_preprocessed_data()
    
    # Step 2: Run analysis (if not skipped)
    if not args.skip_analysis:
        print("\n=== STEP 2: ANALYSIS ===")
        analysis_results = run_analysis(preprocessed_data)
    else:
        print("\n=== SKIPPING ANALYSIS ===")
        analysis_results = load_saved_analysis_results()
    
    # Step 3: Prepare data for webmap and generate webmap (if not skipped)
    if not args.skip_webmap:
        print("\n=== STEP 3: WEBMAP GENERATION ===")
        webmap_data = prepare_data_for_webmap(analysis_results)
        webmap_path = create_webmap(webmap_data)
        print(f"Web map saved to {webmap_path}")
        
        # Open webmap in browser if requested
        if args.open_map and os.path.exists(webmap_path):
            print(f"Opening web map in default browser...")
            webbrowser.open('file://' + os.path.abspath(webmap_path))
    else:
        print("\n=== SKIPPING WEBMAP GENERATION ===")
        webmap_path = str(DATASET_INFO["Output"]["Webmap"]["path"])
        if args.open_map and os.path.exists(webmap_path):
            print(f"Opening existing web map in default browser...")
            webbrowser.open('file://' + os.path.abspath(webmap_path))
    
    # Calculate and display elapsed time
    elapsed_time = time.time() - start_time
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")

if __name__ == "__main__":
    main()