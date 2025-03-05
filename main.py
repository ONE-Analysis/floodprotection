"""
Main execution script for Flood Infrastructure Analysis
Orchestrates the entire workflow: preprocessing, analysis, and webmap generation
"""
import os
import sys
import time
import webbrowser
from pathlib import Path

# Add the project directory to the path
project_dir = os.path.dirname(os.path.abspath(__file__))
if project_dir not in sys.path:
    sys.path.append(project_dir)

# Import project modules
from config import DATASET_INFO, OVERWRITE
from preprocessing import preprocess_data
from analysis import run_analysis
from webmap import generate_webmap

def check_input_files():
    """
    Check if required input files exist
    
    Returns:
    bool: True if all required files exist, False otherwise
    """
    required_files = [
        DATASET_INFO["Input"]["Site_Bounds"]["path"],
        DATASET_INFO["Input"]["Alignment"]["path"],
        DATASET_INFO["Input"]["Buildings"]["path"],
        DATASET_INFO["Input"]["NSI"]["path"],
        DATASET_INFO["Input"]["FEMA_Flood"]["path"]
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(str(file_path)):
            # Check if it's a CAD file, in which case we can convert from GeoJSON if it exists
            if file_path.suffix.lower() in ['.dwg', '.dxf']:
                geojson_path = file_path.with_suffix('.geojson')
                if os.path.exists(str(geojson_path)):
                    continue
            missing_files.append(file_path)
    
    if missing_files:
        print("The following required input files are missing:")
        for file_path in missing_files:
            print(f"  - {file_path}")
        return False
    
    return True

def main():
    """
    Main execution function
    """
    start_time = time.time()
    print("Starting Flood Infrastructure Analysis")
    
    # Check if required input files exist
    if not check_input_files():
        print("Error: Missing required input files. Please provide the missing files and try again.")
        return
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(str(DATASET_INFO["Output"]["Webmap"]["path"])), exist_ok=True)
    
    try:
        print("\n----- PREPROCESSING -----")
        preprocessed_data = preprocess_data()
        
        print("\n----- ANALYSIS -----")
        analysis_results = run_analysis(preprocessed_data)
        
        print("\n----- WEBMAP GENERATION -----")
        webmap_path = generate_webmap(analysis_results)
        
        # Open the webmap in a browser
        print(f"\nOpening webmap: {webmap_path}")
        try:
            webbrowser.open('file://' + os.path.abspath(webmap_path))
        except Exception as e:
            print(f"Error opening webmap in browser: {e}")
            print(f"Please manually open the file: {webmap_path}")
        
        # Print summary
        print("\n----- SUMMARY -----")
        
        # Count buildings impacted in each scenario
        buildings = analysis_results["buildings_scenarios"]
        existing_count = int(buildings['0_2_PctFlood'].sum())
        alignment_count = int(buildings['0_2_PctFlood_Alignment'].sum())
        total_buildings = len(buildings)
        
        print(f"Total buildings analyzed: {total_buildings}")
        print(f"Buildings in 0.2% flood zone (existing): {existing_count} ({existing_count/total_buildings*100:.1f}%)")
        print(f"Buildings in 0.2% flood zone (with alignment): {alignment_count} ({alignment_count/total_buildings*100:.1f}%)")
        
        # Calculate reduction
        if existing_count > 0:
            reduction = existing_count - alignment_count
            reduction_pct = reduction / existing_count * 100
            print(f"Reduction in flooded buildings: {reduction} ({reduction_pct:.1f}% reduction)")
        else:
            print("No buildings were flooded in the existing scenario.")
        
        # Print execution time
        end_time = time.time()
        print(f"\nTotal execution time: {end_time - start_time:.2f} seconds")
        print("Analysis complete! Results are available in the output directory.")
    
    except Exception as e:
        print(f"\nERROR: Analysis failed with error: {e}")
        print("Please fix the issues and try again.")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()