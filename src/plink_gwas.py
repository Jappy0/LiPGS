import subprocess
import argparse

def run_plink(bfile, pheno, covar, covar_name, output, allow_no_sex=False):
    # Build the command
    command = [
        "plink",
        "--bfile", bfile,
        "--pheno", pheno,
        "--covar", covar,
        "--covar-name", covar_name,
        "--logistic",
        "--out", output
    ]
    
    if allow_no_sex:
        command.append("--allow-no-sex")
    
    # Run the command
    try:
        subprocess.run(command, check=True)
        print(f"PLINK analysis completed successfully. Output saved to: {output}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running PLINK: {e}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run GWAS analysis using PLINK")
    parser.add_argument("--bfile", required=True, help="Base filename for input files (without .bed, .bim, .fam)")
    parser.add_argument("--pheno", required=True, help="Phenotype file")
    parser.add_argument("--covar", required=True, help="Covariate file")
    parser.add_argument("--covar-name", required=True, help="Comma-separated list of covariate names")
    parser.add_argument("--output", required=True, help="Output file prefix")
    parser.add_argument("--allow-no-sex", action='store_true', help="Allow analysis without sex information")

    # Parse the arguments
    args = parser.parse_args()

    # Run PLINK with the provided arguments
    run_plink(args.bfile, args.pheno, args.covar, args.covar_name, args.output, args.allow_no_sex)

if __name__ == "__main__":
    main()
