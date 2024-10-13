import os
import subprocess

def merge_bfiles(merge_list_file, output_prefix):
    """
    Merge multiple PLINK binary files into one using PLINK2.

    Parameters:
        merge_list_file (str): Path to the file listing binary files to merge.
        output_prefix (str): Prefix for the output merged binary files.
    """
    # Command to merge binary files
    merge_cmd = [
        "plink2",
        "--pmerge-list", merge_list_file,  # Use --pmerge-list instead of --merge-list
        "--make-bed",
        "--out", output_prefix
    ]
    
    print("Merging PLINK files...")
    subprocess.run(merge_cmd, check=True)
    print(f"Merging completed. Output files prefixed with: {output_prefix}")

def run_plink_qc(bfile_prefix):
    """
    Perform quality control on genotype data using PLINK2.
    
    Parameters:
        bfile_prefix (str): Prefix of the binary files (without .bed/.bim/.fam).
    """
    # Step 1: Initial Quality Control
    qc_step1_cmd = [
        "plink2", 
        "--bfile", bfile_prefix,
        "--maf", "0.01",  # Remove SNPs with MAF < 1%
        "--hwe", "1e-6",  # Remove SNPs deviating from HWE (P < 10^-6)
        "--geno", "0.05", # Remove SNPs with poor genotyping rate (< 95%)
        "--mind", "0.05", # Remove individuals with low genotype rates (< 95%)
        "--make-bed", 
        "--out", f"{bfile_prefix}_qc_step1"
    ]
    
    print("Running PLINK QC Step 1...")
    subprocess.run(qc_step1_cmd, check=True)
    
    # Step 2: Exclude individuals based on sex inconsistencies
    qc_step2_cmd = [
        "plink2", 
        "--bfile", f"{bfile_prefix}_qc_step1",
        "--check-sex", # Check for sex discrepancies
        "--out", f"{bfile_prefix}_qc_step2"
    ]
    
    print("Running PLINK QC Step 2: Check Sex...")
    subprocess.run(qc_step2_cmd, check=True)

    # Read the sex check results and create a list of individuals to exclude
    with open(f"{bfile_prefix}_qc_step2.sexcheck", 'r') as f:
        lines = f.readlines()
        to_exclude = [line.split()[0] for line in lines[1:] if line.strip() and (line.split()[1] != line.split()[2])]

    # Step 3: Create a list of individuals to exclude from the QC
    exclude_file = f"{bfile_prefix}_exclude.txt"
    with open(exclude_file, 'w') as f:
        for individual in to_exclude:
            f.write(f"{individual}\n")

    # Step 4: Perform final QC
    final_qc_cmd = [
        "plink2", 
        "--bfile", f"{bfile_prefix}_qc_step1",
        "--remove", exclude_file,  # Exclude individuals
        "--make-bed",
        "--out", f"{bfile_prefix}_final_qc"
    ]
    
    print("Running PLINK Final QC...")
    subprocess.run(final_qc_cmd, check=True)
    
    print(f"Quality control completed. Final dataset: {bfile_prefix}_final_qc")


# Example usage
if __name__ == "__main__":
    bfile_lst = "/hpcfs/users/a1236780/Repos/LiPGS/data/bfiles_list.txt"
    # Output prefix for the merged PLINK files
    output_prefix = "/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen_combined_european_bfiles"
    
    # Merge the bfiles
    merge_bfiles(bfile_lst, output_prefix)

    run_plink_qc(output_prefix)
