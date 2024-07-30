#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <zip_directory> <template_json>"
    exit 1
fi

# Input arguments
ZIP_DIR="$1"
TEMPLATE_JSON="$2"

# Define the path to the WDL file
WDL_FILE="/home/alex/DataspellProjects/DRAP/workflows/DroppedReadsAnalysis.wdl"

# Iterate over each zip file in the directory
# sort the zip files by their size in ascending order

for zip_file in $(ls -S "$ZIP_DIR"/*.zip); do
    echo "Processing $zip_file"

    # Extract the base name of the zip file (without directory and extension)
    zip_basename=$(basename "$zip_file" .zip)

    # Create an output directory named after the zip file
    OUTPUT_DIR="$zip_basename"
    mkdir -p "$OUTPUT_DIR"

    # Create a new JSON file for this zip file
    JSON_FILE="$OUTPUT_DIR/input.json"

    # Populate the JSON file with the template, replacing placeholders
    jq --arg input_zip "$zip_file" --arg run_name "$zip_basename" \
       '.input_zip = $input_zip | .run_name = $run_name' \
       "$TEMPLATE_JSON" > "$JSON_FILE"

    # Run miniwdl with the WDL file and the generated JSON file
    miniwdl run --no-quant-check "$WDL_FILE" --input "$JSON_FILE" -d "$OUTPUT_DIR"

done
