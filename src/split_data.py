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
    groups = np.array_split(fid_iid_df.sample(frac=1, random_state=42), num_groups)  # Shuffle before splitting
    return groups

def generate_files(groups, plink_prefix, pheno_file, output_dir):
    """Generate target and discovery phenotype data and PLINK binary files."""
    os.makedirs(output_dir, exist_ok=True)

    for i, group in enumerate(groups):
        group_file = os.path.join(output_dir, f"group_{i + 1}.txt")
        group.to_csv(group_file, sep='\t', header=False, index=False)

        # Generate the target and discovery datasets
        generate_plink_files(group_file, plink_prefix, output_dir, i + 1, is_target=True)

        # Discovery data: individuals not in the current group
        discovery_ID_file = os.path.join(output_dir, f"discovery_FID_IID_{i + 1}.txt")
        discovery_pheno_file = os.path.join(output_dir, f"discovery_pheno_{i + 1}.txt")
        generate_discovery_file(group_file, pheno_file, discovery_ID_file, discovery_pheno_file)
        generate_plink_files(discovery_ID_file, plink_prefix, output_dir, i + 1, is_target=False)

def generate_discovery_file(group_file, pheno_file, discovery_ID_file, discovery_pheno_file):
    """Create discovery group by excluding the target individuals and generate discovery phenotype data as integers."""
    # Read the target individuals' FID and IID from the group file
    target_data = pd.read_csv(group_file, sep='\t', header=None, names=['FID', 'IID'])

    # Read the phenotype file with the actual header
    pheno_data = pd.read_csv(pheno_file, sep='\t')

    # Filter out the target individuals from the phenotype data
    discovery_data = pheno_data[~pheno_data[['FID', 'IID']].apply(tuple, axis=1).isin(target_data[['FID', 'IID']].apply(tuple, axis=1))]

    # Convert all phenotype columns to integer type, excluding the first two columns (FID and IID)
    discovery_data.iloc[:, 2:] = discovery_data.iloc[:, 2:].fillna(0).astype(int)

    # Write the discovery FID/IID to the discovery ID file
    discovery_data[['FID', 'IID']].to_csv(discovery_ID_file, sep='\t', header=False, index=False)
    
    # Write the discovery phenotype data (excluding target individuals) to the discovery phenotype file
    discovery_data.to_csv(discovery_pheno_file, sep='\t', header=True, index=False)

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
    generate_files(groups, args.plink_prefix, args.phenotype_file, args.output_dir)

if __name__ == "__main__":
    main()
