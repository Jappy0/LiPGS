#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <mode> <bfile> <pheno_file> <covar_file> <covar_names> <output_prefix>"
    echo "Modes: 'binary' for binary phenotype analysis, 'alda' for ALDA analysis"
    exit 1
fi

# Set input parameters from command line arguments
MODE="$1"            # Analysis mode ('binary' or 'alda')
BFILE="$2"          # Base filename without extension
PHENO_FILE="$3"     # Phenotype file
COVAR_FILE="$4"     # Covariates file
COVAR_NAMES="$5"    # Comma-separated list of covariate names
OUTPUT_PREFIX="$6"  # Output prefix

# Function to run PLINK GWAS analysis for binary phenotypes
run_binary_analysis() {
    plink --bfile "$BFILE" \
          --pheno "$PHENO_FILE" \
          --covar "$COVAR_FILE" \
          --covar-name "$COVAR_NAMES" \
          --allow-no-sex \
          --logistic \
          --adjust \
          --out "$OUTPUT_PREFIX"

    # Check if the analysis was successful
    if [ $? -eq 0 ]; then
        echo "Binary GWAS analysis completed successfully. Output saved to: $OUTPUT_PREFIX"
    else
        echo "An error occurred during the binary GWAS analysis."
    fi
}

# Function to run PLINK GWAS analysis for ALDA (assuming different parameters)
run_alda_analysis() {
    plink --bfile "$BFILE" \
          --pheno "$PHENO_FILE" \
          --covar "$COVAR_FILE" \
          --covar-name "$COVAR_NAMES" \
          --allow-no-sex \
          --linear \  # Example for ALDA analysis; modify as needed
          --out "$OUTPUT_PREFIX"

    # Check if the analysis was successful
    if [ $? -eq 0 ]; then
        echo "ALDA analysis completed successfully. Output saved to: $OUTPUT_PREFIX"
    else
        echo "An error occurred during the ALDA analysis."
    fi
}

# Main logic to select analysis mode
case "$MODE" in
    "binary")
        run_binary_analysis
        ;;
    "alda")
        run_alda_analysis
        ;;
    *)
        echo "Invalid mode: $MODE. Use 'binary' or 'alda'."
        exit 1
        ;;
esac
