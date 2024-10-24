#!/bin/bash

# Function to perform quality control using PLINK
plink_qc() {
    # Parameters
    local bfile_prefix="$1"
    
    # Step 1: Initial Quality Control
    local output="${bfile_prefix}_qc"
    echo "Running PLINK QC Step 1..."
    plink --bfile "$bfile_prefix" \
        --maf 0.01 \
        --hwe 1e-6 \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out "$output"

    # # Step 2: Exclude individuals based on sex inconsistencies
    # local output2="${bfile_prefix}_qc_step2"
    # echo "Running PLINK QC Step 2: Check Sex..."
    # plink --bfile "${output}" \
    #     --check-sex \
    #     --out "$output2"

    # # Read the sex check results and create a list of individuals to exclude
    # local exclude_file="${bfile_prefix}_exclude.txt"
    # echo "Creating exclude file..."
    # awk 'NR > 1 && $2 != $3 { print $1 }' "${output2}.sexcheck" > "$exclude_file"

    # # Step 3: Perform final QC
    # local output_qc="${bfile_prefix}_final_qc"
    # echo "Running PLINK Final QC..."
    # plink --bfile "${output}" \
    #     --remove "$exclude_file" \
    #     --make-bed \
    #     --out "$output_qc"

    # echo "Quality control completed. Final dataset: ${output_qc}"
}

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bfile_prefix> <output_dir>"
    exit 1
fi

# Run the QC function with provided arguments
plink_qc "$1"
