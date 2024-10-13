import os
import argparse
import subprocess
import pandas as pd
import glob
from mpire import WorkerPool

def create_root_dir(args):
    """Create the root output directory."""
    os.makedirs(args.root_dir, exist_ok=True)
    os.makedirs(args.gwas_dir, exist_ok=True)
    os.makedirs(args.fold_data_dir, exist_ok=True)
    os.makedirs(args.prscs_dir, exist_ok=True)
    print(f"Output directory created!")

def extract_phenotype(args):
    """Extract phenotype based on the phenotype type without changing the header, adding 1 to the binary phenotype values, and standardizing the third column for non-binary data with 'z' added to the column name for the standardized column.
    """
    output_file = os.path.join(args.fold_data_dir, "European_phenotype.txt")
    
    # Load the phenotype file as a DataFrame
    phenotype_df = pd.read_csv(args.pheno_file, sep='\t')
    
    # Filter rows where the population is 'European'
    european_df = phenotype_df[phenotype_df.iloc[:, 7] == "European"]  # Assuming 8th column is population info (index 7)
    
    # For binary phenotypes, add 1 to the third column
    if args.pheno_type == 'binary':
        european_df.iloc[:, 2] += 1  # Assuming the third column is at index 2 (0-based index)
    else:
        # Calculate the mean and standard deviation for the column
        phenotype_values = european_df.iloc[:, 3].astype(float)  # Ensure column is float
        mean = phenotype_values.mean()
        std = phenotype_values.std()
        # Perform Z-score normalization on the 3rd column
        european_df.loc[:, european_df.columns[3]] = (phenotype_values - mean) / std
        # Rename the normalized column by adding 'z' prefix
        european_df.rename(columns={european_df.columns[3]: 'z' + european_df.columns[3]}, inplace=True)

    # Select only the first three columns
    output_df = european_df.iloc[:, [0, 1, 2] if args.pheno_type == 'binary' else [0, 1, 3]]
    
    # Save the result to the output file
    output_df.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"Extracted phenotype to: {output_file}")
    return output_file

def split_data(args):
    """Split data for cross-validation."""
    path2split_data = f"{args.src_path}/split_data.py"
    Phenotype_data = extract_phenotype(args)

    command = f"python {path2split_data} --phenotype_file {Phenotype_data} --plink_prefix {args.bfile} --output_dir {args.fold_data_dir} --num_groups {args.num_folds}"
    subprocess.run(command, shell=True, check=True)
    print(f"Data split into {args.num_folds} folds.")

# def extract_ref_allele(bim_file, output_file):
#     # Check if the input .bim file exists
#     if not os.path.exists(bim_file):
#         print(f"Error: {bim_file} does not exist.")
#         return
    
#     # Construct the awk command
#     awk_command = f"awk '{{print $2, $6}}' {bim_file} > {output_file}"
    
#     try:
#         # Run the awk command using subprocess
#         subprocess.run(awk_command, shell=True, check=True)
#         print(f"Reference alleles successfully extracted to {output_file}")
#     except subprocess.CalledProcessError as e:
#         print(f"Error occurred while extracting reference alleles: {e}")

def run_gwas(fold, args):
    """Run GWAS for each fold."""
    b_f = os.path.join(args.fold_data_dir, f"discovery_{fold}")
    # b_f = args.bfile
    pheno = os.path.join(args.fold_data_dir, f"discovery_pheno_{fold}.txt")

    output_prefix = os.path.join(args.gwas_dir, f"gwas_fold_{fold}")

    pgen_output = os.path.join(args.fold_data_dir, f"discovery_pgen_ref_fixed_{fold}") 
    command1 = f"plink2 --bfile {b_f} --ref-allele {b_f}.bim 6 2 --make-pgen --threads {args.threads} --out {pgen_output}"
    try:
        print(f"Running Step 1: Adjust REF/ALT Alleles with command: {command1}")
        subprocess.run(command1, shell=True, check=True)
        print("Step 1 completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during Step 1: {e}")
        return

    # command2 = f"plink2 --pfile {pgen_output} --pheno {pheno} --pheno-name {args.pheno_name} --covar {args.covar_file} --covar-name Age,Sex,chip1,zPC1,zPC2,zPC3,zPC4 --covar-quantile-normalize --glm hide-covar no-firth omit-ref --ci 0.95 --threads {args.threads} --out {output_prefix}"
    command2 = f"plink2 --pfile {pgen_output} --pheno {pheno} --pheno-name {args.pheno_name} --covar {args.covar_file} --covar-name zAge,Sex,chip1,zPC1,zPC2,zPC3,zPC4 --covar-quantile-normalize --glm hide-covar --ci 0.95 --threads {args.threads} --out {output_prefix}"

    try:
        print(f"Running Step 2: Run GWAS with the Adjusted File using command: {command2}")
        subprocess.run(command2, shell=True, check=True)
        print("Step 2 completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during Step 2: {e}")
        return

def extract_gwas_columns(fold, args):
    """
    Extracts columns (SNP, A1, A2, BETA, SE), (SNP, A1, A2, OR, SE), 
    (SNP, A1, A2, BETA, P), or (SNP, A1, A2, OR, P) from a GWAS summary statistics file.
    
    The 'ID' column will be used as SNP, 'REF' as A2, and 'LOG(OR)_SE' as SE.
    """
    gwas_file = os.path.join(args.gwas_dir, f"gwas_fold_{fold}.{args.pheno_name}.glm.{args.model}")
    output_prefix = os.path.join(args.gwas_dir, f"gwas_fold_{fold}.glm.{args.model}") 

    # Read the GWAS summary statistics file with headers
    gwas_data = pd.read_csv(gwas_file, sep=r'\s+')

    # Rename columns to match expected output
    gwas_data.rename(columns={'ID': 'SNP', 'REF': 'A2', 'LOG(OR)_SE': 'SE'}, inplace=True)

    # Check for the presence of required columns
    required_cols = ['SNP', 'A1', 'A2', 'SE']
    if not all(col in gwas_data.columns for col in required_cols):
        print(f"Required columns not found in the file. Current columns are: {gwas_data.columns.tolist()}")
        return

    # Generate four different output files based on available columns
    # 1. SNP, A1, A2, BETA, SE (if BETA column exists)
    if 'BETA' in gwas_data.columns:
        beta_se = gwas_data[['SNP', 'A1', 'A2', 'BETA', 'SE']].dropna(subset=['BETA', 'SE'])
        beta_se_file = f"{output_prefix}_beta_se.txt"
        beta_se.to_csv(beta_se_file, sep='\t', index=False)
        print(f"Generated beta_se file: {beta_se_file}")

    # 2. SNP, A1, A2, OR, SE (if OR column exists)
    if 'OR' in gwas_data.columns:
        or_se = gwas_data[['SNP', 'A1', 'A2', 'OR', 'SE']].dropna(subset=['OR', 'SE'])
        or_se_file = f"{output_prefix}_or_se.txt"
        or_se.to_csv(or_se_file, sep='\t', index=False)
        print(f"Generated or_se file: {or_se_file}")

    # 3. SNP, A1, A2, BETA, P (if BETA and P columns exist)
    if 'BETA' in gwas_data.columns and 'P' in gwas_data.columns:
        beta_p = gwas_data[['SNP', 'A1', 'A2', 'BETA', 'P']].dropna(subset=['BETA', 'P'])
        beta_p_file = f"{output_prefix}_beta_p.txt"
        beta_p.to_csv(beta_p_file, sep='\t', index=False)
        print(f"Generated beta_p file: {beta_p_file}")

    # 4. SNP, A1, A2, OR, P (if OR and P columns exist)
    if 'OR' in gwas_data.columns and 'P' in gwas_data.columns:
        or_p = gwas_data[['SNP', 'A1', 'A2', 'OR', 'P']].dropna(subset=['OR', 'P'])
        or_p_file = f"{output_prefix}_or_p.txt"
        or_p.to_csv(or_p_file, sep='\t', index=False)
        print(f"Generated or_p file: {or_p_file}")

def count_lines_in_file(file_path):
    """Count the number of lines in a text file."""
    try:
        with open(file_path, 'r') as file:
            line_count = sum(1 for _ in file)
        return line_count
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def run_prscs_single_chrom(chrom, root_dir, path_to_prscs, ref_dir, summary_stats, bfile, n_gwas, fold):
    """Run PRScs for a specific chromosome."""
    output_prs = os.path.join(root_dir, f"prs_fold_{fold}_chr{chrom}")
    b_file = bfile + ".bim"
    command = f"python {path_to_prscs} --ref_dir {ref_dir} --bim_file {b_file} --sst_file {summary_stats} --n_gwas {n_gwas} --chrom={chrom} --out {output_prs}"
    subprocess.run(command, shell=True, check=True) 
    print(f"Chromosome {chrom}: PRScs completed successfully.")

def run_prscs(fold, args):
    """Run PRScs in parallel for all chromosomes using mpire."""
    # Determine the number of GWAS samples by counting lines in the discovery sample list (excluding header)
    discovery_sample_id_lst = os.path.join(args.fold_data_dir, f"discovery_FID_IID_{fold}.txt")
    n_gwas = count_lines_in_file(discovery_sample_id_lst) - 1
    
    # Define the PRScs script and summary statistics file based on the model and phenotype type
    ss_prefix = os.path.join(args.gwas_dir, f"gwas_fold_{fold}.glm.{args.model}")
    summary_stats = f"{ss_prefix}_or_p.txt" if args.pheno_type == "binary" else f"{ss_prefix}_beta_p.txt"

    chrom_lst = list(range(1, 23))  # Chromosomes 1 to 22 for autosomes
    prscs_path = f"{args.src_path}/PRScs.py"
    bfile = os.path.join(args.fold_data_dir, f"target_{fold}")
    # Prepare arguments for each chromosome task
    cur_args = [(chrom, args.prscs_dir, prscs_path, args.ref_dir, summary_stats, bfile, n_gwas, fold) for chrom in chrom_lst]
    
    # Use mpire WorkerPool for parallel execution
    with WorkerPool(n_jobs=args.threads) as pool:
        # Execute the run_prscs_single_chrom function in parallel
        pool.map(run_prscs_single_chrom, cur_args, progress_bar=False)
    
    print(f"All chromosomes have been processed in fold {fold}.")
    #########################################################
    snp_files = sorted(glob.glob(os.path.join(args.prscs_dir, f"prs_fold_{fold}_chr*_pst_eff_*.txt")))
    all_snp_file = os.path.join(args.prscs_dir, f"prs_fold_{fold}_chr1_22.txt")
    # Open the output file for writing
    with open(all_snp_file, 'w') as outfile:
        for file in snp_files:
            with open(file, 'r') as infile:
                # Write the contents of each file to the output file
                outfile.writelines(infile.readlines())
    print(f"Concatenation complete. Combined file saved as {all_snp_file}.")

    plink_prs = os.path.join(args.prscs_dir, f"PGS_fold_{fold}_sum")
    command_plink = f"plink2 --bfile {bfile} --score {all_snp_file} 2 4 6 --out {plink_prs}"
    subprocess.run(command_plink, shell=True, check=True)
    print(f"PRS calculated for fold {fold}. Output: {plink_prs}")

def merge_files(args):
    # Read the main file
    df_main = pd.read_csv(args.pheno_file, sep='\t' )
    # Initialize a list to store all merged dataframes
    merged_dfs = []
    df_lst = []

    for fold in range(1, args.num_folds+1):
        file_path = os.path.join(args.prscs_dir, f"PGS_fold_{fold}_sum.sscore")
        cur_df = pd.read_csv(file_path, sep=r'\s+', header=0)
        cur_df.rename(columns={'#FID': 'FID'}, inplace=True)
        df_lst.append(cur_df)
    # Concatenate all DataFrames in the list into a single DataFrame
    combined_df = pd.concat(df_lst, ignore_index=True)
    df_merged = pd.DataFrame()
    if args.pheno_type == "alda":
        df_zpheno = pd.read_csv(os.path.join(args.fold_data_dir, "European_phenotype.txt"), sep='\t')
        print(df_zpheno)
        df_merged = pd.merge(df_main, df_zpheno, on=['FID', 'IID'], how='left')
        df_merged = pd.merge(df_merged, combined_df, on=['FID', 'IID'], how='left')
    else:
        df_merged = pd.merge(df_main, combined_df, on=['FID', 'IID'], how='left')
    df_merged.rename(columns={'SCORE1_AVG': f'PRScs_{args.pheno_type}'}, inplace=True)
    output_base_name = os.path.join(args.root_dir, f'ConLiGEN_phenotypes_{args.pheno_type}_with_PRSs')
    # Calculate the mean and standard deviation 
    mean_prs = df_merged[f'PRScs_{args.pheno_type}'].mean()
    std_prs = df_merged[f'PRScs_{args.pheno_type}'].std()
    # Step 4: Apply Z-standardization
    df_merged[f'Z_PRScs_{args.pheno_type}'] = (df_merged[f'PRScs_{args.pheno_type}'] - mean_prs) / std_prs

    # Save the final dataframe to a tab-separated text file
    df_merged.to_csv(f'{output_base_name}.txt', sep='\t', index=False)
    # Save the final dataframe to a comma-separated CSV file
    df_merged.to_csv(f'{output_base_name}.csv', sep=',', index=False)

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run the GWAS and PRS pipeline")
    parser.add_argument('--bfile', required=True, help="Path to PLINK bfile prefix")
    parser.add_argument('--pheno_file', required=True, help="Path to phenotype file")
    parser.add_argument('--pheno_name', required=True, help="Phenotype name in the file of --pheno_file")
    parser.add_argument('--pheno_type', choices=['binary', 'alda'], required=True, help="Type of phenotype")
    parser.add_argument('--covar_file', required=True, help="Path to covariates file")
    parser.add_argument('--ref_dir', required=True, help="Path to PRScs reference directory")
    parser.add_argument('--src_path', required=True, help="Path to PRScs.py script")
    # parser.add_argument('--path_to_prscs', required=True, help="Path to PRScs.py script")
    parser.add_argument('--num_folds', type=int, default=5, help="Number of folds for cross-validation")
    parser.add_argument('--root_dir', default='output', help="Root output directory")
    parser.add_argument('--fold_data_dir', help="Root output directory")
    parser.add_argument('--gwas_dir', help="Root output directory")
    parser.add_argument('--prscs_dir', help="Root output directory")
    parser.add_argument('--model', required=True, help="Model used for gwas.")
    parser.add_argument('--threads', type=int, default=1, help="Threads or processes number for PLINK2 or PRScs.")
    
    args = parser.parse_args()
    # Step 1: Create root directory
    args.gwas_dir = os.path.join(args.root_dir, "gwas")
    args.fold_data_dir = os.path.join(args.root_dir, "fold_data")
    args.prscs_dir = os.path.join(args.root_dir, "prscs")
    create_root_dir(args)
    # Step 2: extract and Split data
    split_data(args)
    if args.pheno_type != "binary":
        args.pheno_name = "zAldaTOTAL"
    # Step 4: Run GWAS for each fold
    for fold in range(1, args.num_folds+1):
        run_gwas(fold, args)
        extract_gwas_columns(fold, args)
        run_prscs(fold, args)

    merge_files(args)

if __name__ == "__main__":
    main()
