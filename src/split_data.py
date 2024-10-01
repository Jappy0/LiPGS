import pandas as pd
import numpy as np
import os
import argparse

def extract_fid_iid(phenotype_file):
    """Extract FID and IID from phenotype file."""
    phenotype_df = pd.read_csv(phenotype_file, sep=r'\s+')  # Updated from delim_whitespace=True
    fid_iid_df = phenotype_df[['FID', 'IID']]
    return fid_iid_df

def split_into_groups(fid_iid_df, num_groups=5):
    """Split the FID/IID data into specified number of groups."""
    groups = np.array_split(fid_iid_df.sample(frac=1), num_groups)  # Shuffle before splitting
    return groups

def generate_files(groups, phenotype_file, plink_prefix, output_dir):
    """Generate target and discovery phenotype data and PLINK binary files."""
    os.makedirs(output_dir, exist_ok=True)

    for i, group in enumerate(groups):
        group_file = os.path.join(output_dir, f"group_{i + 1}.txt")
        group.to_csv(group_file, sep='\t', header=False, index=False)

        # Generate the target and discovery datasets
        generate_plink_files(group_file, plink_prefix, output_dir, i + 1, is_target=True)

        # Discovery data: individuals not in the current group
        discovery_file = os.path.join(output_dir, f"discovery_{i + 1}.txt")
        generate_discovery_file(group_file, plink_prefix, discovery_file)
        generate_plink_files(discovery_file, plink_prefix, output_dir, i + 1, is_target=False)

def generate_discovery_file(group_file, plink_prefix, discovery_file):
    """Create discovery group by excluding the target individuals."""
    target_data = pd.read_csv(group_file, sep='\t', header=None, names=['FID', 'IID'])

    # Read FIDs and IIDs from the .fam file
    fam_data = pd.read_csv(f"{plink_prefix}.fam", sep=r'\s+', header=None, names=['FID', 'IID', 'p1', 'p2', 'sex', 'pheno'])
    
    # Get individuals not in the target group (i.e., discovery group)
    discovery_data = fam_data[~fam_data[['FID', 'IID']].apply(tuple, axis=1).isin(target_data[['FID', 'IID']].apply(tuple, axis=1))]

    # Write the discovery FID/IID to a file
    discovery_data[['FID', 'IID']].to_csv(discovery_file, sep='\t', header=False, index=False)

def generate_plink_files(group_file, plink_prefix, output_dir, group_index, is_target):
    """Generate PLINK binary files for either the target or discovery groups."""
    if is_target:
        output_file_prefix = os.path.join(output_dir, f"target_{group_index}")
    else:
        output_file_prefix = os.path.join(output_dir, f"discovery_{group_index}")

    # Use PLINK to create the target or discovery binary files
    os.system(f"plink --bfile {plink_prefix} --keep {group_file} --make-bed --out {output_file_prefix}")

def main():
    parser = argparse.ArgumentParser(description="Split phenotype and genotype data into groups.")
    parser.add_argument("--phenotype_file", type=str, help="Path to the phenotype data file.")
    parser.add_argument("--plink_prefix", type=str, help="Prefix for the PLINK binary files (without extensions).")
    parser.add_argument("--output_dir", type=str, help="Directory where output files will be saved.")
    parser.add_argument("--num_groups", type=int, default=5, help="Number of groups to split into (default: 5).")

    args = parser.parse_args()

    fid_iid_df = extract_fid_iid(args.phenotype_file)
    groups = split_into_groups(fid_iid_df, num_groups=args.num_groups)
    generate_files(groups, args.phenotype_file, args.plink_prefix, args.output_dir)

if __name__ == "__main__":
    main()
