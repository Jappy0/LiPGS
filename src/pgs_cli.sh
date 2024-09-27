#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <prscs_path> <bfile> <weight_file> <output_prefix> <pheno_file> <covar_file>"
    echo "This script calculates PGS using PRScs."
    exit 1
fi

# Set input parameters from command line arguments
PRScs_PATH="$1"     # Path to the PRScs.py script
BFILE="$2"          # Base filename without extension
WEIGHT_FILE="$3"    # PGS weight file
OUTPUT_PREFIX="$4"  # Output prefix
PHENO_FILE="$5"     # Phenotype file
COVAR_FILE="$6"     # Covariates file

# Function to run PRScs for calculating Polygenic Scores
run_prscs() {
    # Run PRScs command directly with the full path to PRScs.py
    python "$PRScs_PATH" \
        --bfile "$BFILE" \
        --weight "$WEIGHT_FILE" \
        --out "$OUTPUT_PREFIX" \
        --pheno "$PHENO_FILE" \
        --covar "$COVAR_FILE" \
        --allow-no-sex

    # Check if the analysis was successful
    if [ $? -eq 0 ]; then
        echo "PGS calculation completed successfully. Output saved to: $OUTPUT_PREFIX"
    else
        echo "An error occurred during the PGS calculation."
    fi
}

# Main logic to run the PGS calculation
run_prscs
