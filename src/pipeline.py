import os
import argparse
import subprocess
import pandas as pd
import glob
import multiprocessing

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

def extract_ref_allele(bim_file, output_file):
    # Check if the input .bim file exists
    if not os.path.exists(bim_file):
        print(f"Error: {bim_file} does not exist.")
        return
    
    # Construct the awk command
    awk_command = f"awk '{{print $2, $6}}' {bim_file} > {output_file}"
    
    try:
        # Run the awk command using subprocess
        subprocess.run(awk_command, shell=True, check=True)
        print(f"Reference alleles successfully extracted to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while extracting reference alleles: {e}")

def run_gwas(fold, root_dir, covar, pheno_name, threads_num):
    """Run GWAS for each fold."""
    # bed_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.bed")
    # bim_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.bim")
    # fam_f = os.path.join(root_dir, f"fold_data/discovery_{fold}.fam")
    b_f = os.path.join(root_dir, f"fold_data/discovery_{fold}")
    pheno = os.path.join(root_dir, f"fold_data/discovery_pheno_{fold}.txt")
    output_prefix = os.path.join(root_dir, f"gwas_fold_{fold}")

    # analysis_type = "--logistic --hide-covar" if pheno_type == 'binary' else "--linear --hide-covar"
    # command = f"plink --bfile {b_f} --pheno {pheno} --covar {covar} --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --keep-allele-order --real-ref-alleles {analysis_type} --adjust --all-pheno --ci 0.95 --make-bed --out {output_prefix}"
    
    # ref_allele = os.path.join(root_dir, f"fold_data/ref_allele_{fold}.txt")
    # extract_ref_allele(b_f+".bim", ref_allele)
    pgen_output = os.path.join(root_dir, f"fold_data/discovery_pgen_ref_fixed_{fold}") 
    command1 = f"plink2 --bfile {b_f} --ref-allele {b_f}.bim 6 2 --make-pgen --threads {threads_num} --out {pgen_output}"
    try:
        print(f"Running Step 1: Adjust REF/ALT Alleles with command: {command1}")
        subprocess.run(command1, shell=True, check=True)
        print("Step 1 completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during Step 1: {e}")
        return

    command2 = f"plink2 --pfile {pgen_output} --pheno {pheno} --pheno-name {pheno_name} --covar {covar} --covar-name Age,Sex,chip1,zPC1,zPC2,zPC3,zPC4 --covar-quantile-normalize --glm hide-covar no-firth omit-ref --ci 0.95 --threads {threads_num} --out {output_prefix}"#--pfile {pgen_output} --make-pgen 

    try:
        print(f"Running Step 2: Run GWAS with the Adjusted File using command: {command2}")
        subprocess.run(command2, shell=True, check=True)
        print("Step 2 completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during Step 2: {e}")
        return

    # analysis_type = "Response" if pheno_type == 'binary' else "--linear"
    # command = f"plink2 --bfile {b_f} --ref-allele {ref_allele} 2 1 --make-pgen --pheno {pheno} --pheno-name Response --covar {covar} --covar-name Age,Sex,zPC1,zPC2,zPC3,zPC4 --covar-variance-standardize hide-covar --glm firth-fallback --ci 0.95 --out {output_prefix}"#--covar-variance-standardize 
    # command = f"plink --bed {bed_f} --bim {bim_f} --fam {fam_f} --pheno {pheno} --covar {covar} --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --keep-allele-order --real-ref-alleles {analysis_type} --adjust --all-pheno --ci 0.95 --make-bed --out {output_prefix}"
    # subprocess.run(command, shell=True, check=True)
    # print(f"GWAS completed for fold {fold}. Output: {output_prefix}.assoc.{pheno_type}")

def extract_gwas_columns(gwas_file, output_prefix):
    """
    Extracts columns (SNP, A1, A2, BETA, SE), (SNP, A1, A2, OR, SE), 
    (SNP, A1, A2, BETA, P), or (SNP, A1, A2, OR, P) from a GWAS summary statistics file.
    
    The 'ID' column will be used as SNP, 'REF' as A2, and 'LOG(OR)_SE' as SE.
    """
    
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
    # or_se = gwas_data[['SNP', 'A1', 'A2', 'OR', 'SE']]
    or_se = gwas_data[['SNP', 'A1', 'A2', 'OR', 'SE']].dropna(subset=['OR', 'SE'])
    or_se_file = f"{output_prefix}_or_se.txt"
    or_se.to_csv(or_se_file, sep='\t', index=False)
    print(f"Generated or_se file: {or_se_file}")

    # # 3. beta_p file (SNP, A1, A2, BETA, P)
    # beta_p = gwas_data[['SNP', 'A1', 'A2', 'BETA', 'P']]
    # beta_p_file = f"{output_prefix}_beta_p.txt"
    # beta_p.to_csv(beta_p_file, sep='\t', index=False)
    # print(f"Generated beta_p file: {beta_p_file}")

    # 4. or_p file (SNP, A1, A2, OR, P)
    # or_p = gwas_data[['SNP', 'A1', 'A2', 'OR', 'P']]
    or_p = gwas_data[['SNP', 'A1', 'A2', 'OR', 'P']].dropna(subset=['OR', 'P'])
    or_p_file = f"{output_prefix}_or_p.txt"
    or_p.to_csv(or_p_file, sep='\t', index=False)
    print(f"Generated or_p file: {or_p_file}")

# Function to run PRScs for a specific chromosome
def run_prscs_single_chrom(chrom, root_dir, path_to_prscs, ref_dir, summary_stats, bfile, n_gwas, fold):
    """Run PRScs for a specific chromosome."""
    output_prs = os.path.join(root_dir, f"prs_fold_{fold}_chr{chrom}")
    b_file = bfile + ".bim"
    command = f"python {path_to_prscs} --ref_dir {ref_dir} --bim_file {b_file} --sst_file {summary_stats} --n_gwas {n_gwas} --chrom={chrom} --out {output_prs}"

    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Chromosome {chrom}: PRScs completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Chromosome {chrom}: PRScs failed. Error: {e}")

def run_prscs(root_dir, path_to_prscs, ref_dir, summary_stats, bfile, n_gwas, fold, threads_num):
    """Run PRScs in parallel for all chromosomes."""
    chrom_lst = list(range(1, 23))  # Chromosomes 1 to 22 for autosomes
    # Create a pool of workers, one for each chromosome
    pool = multiprocessing.Pool(processes=threads_num)
    args = [(chrom, root_dir, path_to_prscs, ref_dir, summary_stats, bfile, n_gwas, fold) for chrom in chrom_lst]
    # Run PRScs in parallel for each chromosome
    pool.starmap(run_prscs_single_chrom, args)
    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()

    #################################################################
    ## no parallel
    # """Run PRScs to calculate PRS."""
    # output_prs = os.path.join(root_dir, f"prs_fold_{fold}")
    # b_file = bfile + ".bim"
    # command = f"python {path_to_prscs} --ref_dir {ref_dir} --bim_file {b_file} --sst_file {summary_stats} --n_gwas {n_gwas} --chrom=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 --out {output_prs}"
    # subprocess.run(command, shell=True, check=True)
    #########################################################
    snp_files = sorted(glob.glob(os.path.join(root_dir, f"prs_fold_{fold}_chr*_pst_eff_*.txt")))
    all_snp_file = os.path.join(root_dir, f"prs_fold_{fold}_chr1_22.txt")
    # Open the output file for writing
    with open(all_snp_file, 'w') as outfile:
        for file in snp_files:
            with open(file, 'r') as infile:
                # Write the contents of each file to the output file
                outfile.writelines(infile.readlines())
    print(f"Concatenation complete. Combined file saved as {all_snp_file}.")

    plink_prs = os.path.join(root_dir, f"PGS_fold_{fold}_sum")
    command_plink = f"plink2 --bfile {bfile} --score {all_snp_file} 2 4 6 sum --out {plink_prs}"
    subprocess.run(command_plink, shell=True, check=True)
    print(f"PRS calculated for fold {fold}. Output: {plink_prs}")

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Run the GWAS and PRS pipeline")
    parser.add_argument('--bfile', required=True, help="Path to PLINK bfile prefix")
    parser.add_argument('--pheno_file', required=True, help="Path to phenotype file")
    parser.add_argument('--pheno_name', required=True, help="Phenotype name in the file of --pheno_file")
    parser.add_argument('--covar_file', required=True, help="Path to covariates file")
    parser.add_argument('--pheno_type', choices=['binary', 'continuous'], required=True, help="Type of phenotype")
    parser.add_argument('--ref_dir', required=True, help="Path to PRScs reference directory")
    parser.add_argument('--path_to_prscs', required=True, help="Path to PRScs.py script")
    parser.add_argument('--num_folds', type=int, default=5, help="Number of folds for cross-validation")
    parser.add_argument('--root_dir', default='output', help="Root output directory")
    parser.add_argument('--threads', type=int, default=1, help="Threads or processes number for PLINK2 or PRScs.")

    args = parser.parse_args()

    # Step 1: Create root directory
    create_root_dir(args.root_dir)
    # Step 2: Extract phenotype
    phenotype_file = extract_phenotype(args.pheno_file, args.root_dir, args.pheno_type)
    # Step 3: Split data
    split_data(phenotype_file, args.bfile, os.path.join(args.root_dir, "fold_data"), args.num_folds)

    # Step 4: Run GWAS for each fold
    for fold in range(1, args.num_folds+1):
    # for fold in range(2, args.num_folds+1):
        run_gwas(fold, args.root_dir, args.covar_file, args.pheno_name, args.threads)
        # Step 5: Calculate PRS for each fold
        summary_stats = os.path.join(args.root_dir, f"gwas_fold_{fold}.Response.glm.logistic")
        summary_stats_prefix = os.path.join(args.root_dir, f"gwas_fold_{fold}.glm.logistic")
        # filter_gwas(summary_stats, args.bfile+".bim", summary_stats_prefix)
        extract_gwas_columns(summary_stats, summary_stats_prefix)
        target_bfile = os.path.join(args.root_dir, f"fold_data/target_{fold}")
        run_prscs(args.root_dir, args.path_to_prscs, args.ref_dir, summary_stats_prefix+"_or_p.txt", target_bfile, 1893, fold, args.threads)

if __name__ == "__main__":
    main()
