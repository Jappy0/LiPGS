#!/bin/bash

plink_merge_with_exclude() {
    d1_prefix=$1
    d2_prefix=$2
    output_prefix=$3

    # First attempt to merge the datasets
    echo "Attempting first merge..."
    
    plink \
        --bfile "$d1_prefix" \
        --bmerge "${d2_prefix}.bed" "${d2_prefix}.bim" "${d2_prefix}.fam" \
        --make-bed \
        --allow-no-sex \
        --out "$output_prefix"

    missnp_file="${output_prefix}-merge.missnp"
    
    if [[ -f "$missnp_file" ]]; then
        echo "Merge failed due to problematic SNPs. Excluding SNPs listed in $missnp_file."

        # Exclude problematic SNPs from both datasets
        d1_cleaned_prefix="${d1_prefix}_snp_excluded"
        d2_cleaned_prefix="${d2_prefix}_snp_excluded"
        
        echo "Excluding problematic SNPs from dataset 1..."
        plink \
            --bfile "$d1_prefix" \
            --exclude "$missnp_file" \
            --make-bed \
            --allow-no-sex \
            --out "$d1_cleaned_prefix"

        echo "Excluding problematic SNPs from dataset 2..."
        plink \
            --bfile "$d2_prefix" \
            --exclude "$missnp_file" \
            --make-bed \
            --allow-no-sex \
            --out "$d2_cleaned_prefix"

        # Second attempt to merge
        echo "Attempting second merge after excluding problematic SNPs..."
        plink \
            --bfile "$d1_cleaned_prefix" \
            --bmerge "${d2_cleaned_prefix}.bed" "${d2_cleaned_prefix}.bim" "${d2_cleaned_prefix}.fam" \
            --make-bed \
            --allow-no-sex \
            --out "$output_prefix"
    fi
}

merge_bfiles() {
    dataset_file=$1
    input_dir=$2
    output_prefix=$3

    # Step 1: Read dataset prefixes from the specified file
    mapfile -t dataset_list < "$dataset_file"

    # Ensure there are at least two datasets to merge
    if [[ ${#dataset_list[@]} -lt 2 ]]; then
        echo "At least two datasets are required to perform the merge."
        exit 1
    fi

    # Step 2: Merge the first two datasets
    d0_prefix="${input_dir}/${dataset_list[0]}"
    d1_prefix="${input_dir}/${dataset_list[1]}"
    d2_prefix="${input_dir}/${dataset_list[2]}"
    d3_prefix="${input_dir}/${dataset_list[3]}"
    d4_prefix="${input_dir}/${dataset_list[4]}"
    d5_prefix="${input_dir}/${dataset_list[5]}"
    d6_prefix="${input_dir}/${dataset_list[6]}"

    merge_1="${input_dir}/${output_prefix}_1"
    plink_merge_with_exclude "$d0_prefix" "$d1_prefix" "$merge_1"

    merge_2="${output_prefix}_2"
    plink_merge_with_exclude "$merge_1" "$d2_prefix" "$merge_2"

    merge_3="${output_prefix}_3"
    plink_merge_with_exclude "$merge_2" "$d3_prefix" "$merge_3"

    merge_4="${output_prefix}_4"
    plink_merge_with_exclude "$merge_3" "$d4_prefix" "$merge_4"

    merge_5="${output_prefix}_5"
    plink_merge_with_exclude "$merge_4" "$d5_prefix" "$merge_5"

    merge_6="${output_prefix}_6"
    plink_merge_with_exclude "$merge_5" "$d6_prefix" "$merge_6"
}

# Check for correct usage
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <dataset_file> <input_dir> <output_prefix>"
    exit 1
fi

# Call the merge_bfiles function with the provided arguments
merge_bfiles "$1" "$2" "$3"
