import os
import argparse
import subprocess
import pandas as pd
import glob

def create_root_dir(root_dir):
    """Create the root output directory."""
    os.makedirs(root_dir, exist_ok=True)
    print(f"Created directory: {root_dir}")

def extract_phenotype(pheno_file, root_dir, pheno_type):
    """Extract phenotype based on the phenotype type without changing the header, adding 1 to the phenotype values, and printing only the first three columns of the header."""
    output_file = os.path.join(root_dir, "European_phenotype.txt")
    
    if pheno_type == 'binary':
        # Add 1 to the third column but only print first three columns of the header
        command = f"awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} NR==1 {{print $1, $2, $3}} NR>1 && $8 == \"European\" {{print $1, $2, $3+1}}' {pheno_file} > {output_file}"
    else:
        # Add 1 to the fourth column but only print first three columns of the header
        command = f"awk 'BEGIN {{FS=\"\\t\"; OFS=\"\\t\"}} NR==1 {{print $1, $2, $3}} NR>1 && $8 == \"European\" {{print $1, $2, $4+1}}' {pheno_file} > {output_file}"
    
    subprocess.run(command, shell=True, check=True)
    print(f"Extracted phenotype to: {output_file}")
    return output_file


def split_data(phenotype_file, bfile, output_dir, num_groups):
    """Split data for cross-validation."""
    command = f"python split_data.py --phenotype_file {phenotype_file} --plink_prefix {bfile} --output_dir {output_dir} --num_groups {num_groups}"
    subprocess.run(command, shell=True, check=True)
    print(f"Data split into {num_groups} folds.")

def run_gwas(fold, root_dir, covar, pheno_type):
    """Run GWAS for each fold."""
    # bed_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.bed")
    # bim_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.bim")
    # fam_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.fam")
    b_f = os.path.join(root_dir, f"fold_data/discovery_{fold}")
    pheno = os.path.join(root_dir, f"fold_data/discovery_pheno_{fold}.txt")
    output_prefix = os.path.join(root_dir, f"gwas_fold_{fold}")

    analysis_type = "--logistic" if pheno_type == 'binary' else "--linear"
    command = f"plink --bfile {b_f} --pheno {pheno} --covar {covar} --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --keep-allele-order --real-ref-alleles {analysis_type} --adjust --all-pheno --ci 0.95 --make-bed --out {output_prefix}"
    # command = f"plink --bed {bed_f} --bim {bim_f} --fam {fam_f} --pheno {pheno} --covar {covar} --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --keep-allele-order --real-ref-alleles {analysis_type} --adjust --all-pheno --ci 0.95 --make-bed --out {output_prefix}"
    subprocess.run(command, shell=True, check=True)
    print(f"GWAS completed for fold {fold}. Output: {output_prefix}.assoc.{pheno_type}")

def filter_gwas(gwas_file, bim_file, output_prefix):
    """Generate four different GWAS summary statistics files based on BIM data."""
    
    # Read the BIM file (no header; second column contains SNP IDs, sixth column contains A2 values)
    bim_data = pd.read_csv(bim_file, sep=r'\s+', header=None, usecols=[1, 5], names=['SNP', 'A2'])
    
    # Read the GWAS summary statistics file with headers
    gwas_data = pd.read_csv(gwas_file, sep=r'\s+')
    
    # Check if 'SNP' and necessary columns are in the GWAS data
    required_cols = ['SNP', 'A1', 'SE', 'P', 'OR'] #'BETA', 
    for col in required_cols:
        if col not in gwas_data.columns:
            print(f"GWAS Data does not contain '{col}' column. Current columns are:", gwas_data.columns.tolist())
            return

    # Create a new column 'A2' in the GWAS data by mapping SNPs to their A2 values from BIM data
    gwas_data['A2'] = gwas_data['SNP'].map(bim_data.set_index('SNP')['A2'])

    # Generate four different files based on the requested columns
    
    # # 1. beta_se file (SNP, A1, A2, BETA, SE)
    # beta_se = gwas_data[['SNP', 'A1', 'A2', 'BETA', 'SE']]
    # beta_se_file = f"{output_prefix}_beta_se.txt"
    # beta_se.to_csv(beta_se_file, sep='\t', index=False)
    # print(f"Generated beta_se file: {beta_se_file}")

    # 2. or_se file (SNP, A1, A2, OR, SE)
    or_se = gwas_data[['SNP', 'A1', 'A2', 'OR', 'SE']]
    or_se_file = f"{output_prefix}_or_se.txt"
    or_se.to_csv(or_se_file, sep='\t', index=False)
    print(f"Generated or_se file: {or_se_file}")

    # # 3. beta_p file (SNP, A1, A2, BETA, P)
    # beta_p = gwas_data[['SNP', 'A1', 'A2', 'BETA', 'P']]
    # beta_p_file = f"{output_prefix}_beta_p.txt"
    # beta_p.to_csv(beta_p_file, sep='\t', index=False)
    # print(f"Generated beta_p file: {beta_p_file}")

    # 4. or_p file (SNP, A1, A2, OR, P)
    or_p = gwas_data[['SNP', 'A1', 'A2', 'OR', 'P']]
    or_p_file = f"{output_prefix}_or_p.txt"
    or_p.to_csv(or_p_file, sep='\t', index=False)
    print(f"Generated or_p file: {or_p_file}")
    # Return the number of GWAS data points (rows in gwas_data)
    num_gwas = gwas_data.shape[0]
    print(f"Number of GWAS data points: {num_gwas}")
    
    return num_gwas


# Example usage:
# filter_gwas("gwas_summary.txt", "example.bim", "output_gwas")


# def filter_gwas(gwas_file, bim_file, output_file):
#     """Add A2 values to the GWAS summary statistics file based on SNP values from the BIM file."""
    
#     # Read the BIM file (no header; second column contains SNP IDs, sixth column contains A2 values)
#     bim_data = pd.read_csv(bim_file, sep=r'\s+', header=None, usecols=[1, 5], names=['SNP', 'A2'])
    
#     # Read the GWAS summary statistics file with headers
#     gwas_data = pd.read_csv(gwas_file, sep=r'\s+')
    
#     # Check if 'SNP' is in the columns of the GWAS data
#     if 'SNP' not in gwas_data.columns:
#         print("GWAS Data does not contain 'SNP' column. Current columns are:", gwas_data.columns.tolist())
#         return

#     # Create a new column 'A2' in the GWAS data by mapping SNPs to their A2 values from BIM data
#     gwas_data['A2'] = gwas_data['SNP'].map(bim_data.set_index('SNP')['A2'])

#     # Reorder columns to insert A2 after A1
#     cols = list(gwas_data.columns)
#     a1_index = cols.index('A1')
#     cols.insert(a1_index + 1, 'A2')  # Insert A2 after A1
#     gwas_data = gwas_data[cols]

#     # Write the updated GWAS data to the output file
#     gwas_data.to_csv(output_file, sep='\t', index=False)

#     print(f"Added A2 values to GWAS data. Output saved to: {output_file}")

def run_prscs(root_dir, path_to_prscs, ref_dir, summary_stats, bfile, n_gwas, fold):
    """Run PRScs to calculate PRS."""
    output_prs = os.path.join(root_dir, f"prs_fold_{fold}")
    plink_prs = os.path.join(root_dir, f"plink_prs_fold_{fold}.txt")
    b_file = bfile + ".bim"
    command = f"python {path_to_prscs} --ref_dir {ref_dir} --bim_file {b_file} --sst_file {summary_stats} --n_gwas {n_gwas} --out {output_prs}"
    subprocess.run(command, shell=True, check=True)

    prs_eff = combine_plink_prs(root_dir, fold)

    # Step 3: Run PLINK to score the PRS
    command = f"plink --bfile {bfile} --score {prs_eff} 2 4 6 --out {plink_prs}"
    subprocess.run(command, shell=True, check=True)

    print(f"PRS calculated for fold {fold}. Output: {plink_prs}")
    

def combine_plink_prs(root_dir, fold):
    """Combine PLINK PRS files for each fold into a single file."""
    combined_prs_file = os.path.join(root_dir, f"combined_plink_prs_fold_{fold}.txt")
    
    # Find all PLINK PRS files for the current fold
    plink_prs_files = glob.glob(os.path.join(root_dir, f"prs_fold_{fold}*.txt"))

    # Read and combine data from all found files
    combined_data = []
    for file in plink_prs_files:
        df = pd.read_csv(file, sep='\t', header=0)  # Adjust sep if needed
        combined_data.append(df)

    # Concatenate all DataFrames into a single DataFrame
    if combined_data:
        combined_df = pd.concat(combined_data, ignore_index=True)
        combined_df.to_csv(combined_prs_file, sep='\t', index=False)
        print(f"Combined PLINK PRS results saved to: {combined_prs_file}")
    return combined_prs_file

def combine_prs(root_dir, num_folds):
    """Combine all PRS scores into a single output."""
    prs_files = [os.path.join(root_dir, f"combined_plink_prs_fold_{fold}.txt") for fold in range(1, num_folds+1)]
    combined_prs_file = os.path.join(root_dir, "combined_prs.txt")

    with open(combined_prs_file, 'w') as outfile:
        for prs_file in prs_files:
            with open(prs_file) as infile:
                outfile.write(infile.read())
    print(f"Combined PRS scores into: {combined_prs_file}")

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run the GWAS and PRS pipeline")
    parser.add_argument('--bfile', required=True, help="Path to PLINK bfile prefix")
    parser.add_argument('--pheno_file', required=True, help="Path to phenotype file")
    parser.add_argument('--covar_file', required=True, help="Path to covariates file")
    parser.add_argument('--pheno_type', choices=['binary', 'continuous'], required=True, help="Type of phenotype")
    parser.add_argument('--ref_dir', required=True, help="Path to PRScs reference directory")
    parser.add_argument('--path_to_prscs', required=True, help="Path to PRScs.py script")
    parser.add_argument('--num_folds', type=int, default=5, help="Number of folds for cross-validation")
    parser.add_argument('--root_dir', default='output', help="Root output directory")

    args = parser.parse_args()

    # Step 1: Create root directory
    create_root_dir(args.root_dir)
    # Step 2: Extract phenotype
    phenotype_file = extract_phenotype(args.pheno_file, args.root_dir, args.pheno_type)
    # Step 3: Split data
    split_data(phenotype_file, args.bfile, os.path.join(args.root_dir, "fold_data"), args.num_folds)

    # Step 4: Run GWAS for each fold
    for fold in range(1, args.num_folds+1):
        run_gwas(fold, args.root_dir, args.covar_file, args.pheno_type)
        # Step 5: Calculate PRS for each fold
        summary_stats = os.path.join(args.root_dir, f"gwas_fold_{fold}.Response.assoc.logistic")
        target_bfile = os.path.join(args.root_dir, f"fold_data/target_{fold}")
        summary_stats_prefix = os.path.join(args.root_dir, f"gwas_fold_{fold}.assoc.logistic")
        num_gwas = filter_gwas(summary_stats, args.bfile+".bim", summary_stats_prefix)
        run_prscs(args.root_dir, args.path_to_prscs, args.ref_dir, summary_stats_prefix+"_or_se.txt", target_bfile, num_gwas, fold)

    # Step 6: Combine PRS scores
    combine_prs(args.root_dir, args.num_folds)

if __name__ == "__main__":
    main()
