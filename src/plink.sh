#!/bin/bash
#SBATCH --job-name=plink_qc         # Job name
#SBATCH --output=job_%j.out        # Output file (%j expands to job ID)
#SBATCH --error=job_%j.err         # Error file (%j expands to job ID)
#SBATCH -p batch                    # partition 
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        
#SBATCH --cpus-per-task=12          # Number of CPU cores per task
#SBATCH --mem=64G                   # Memory required per node (adjust based on your data)
#SBATCH --time=5:00:00             # Time limit (hours:minutes:seconds)

module load Anaconda3/2024.06-1
source activate PGS

# Set variables for file names
INPUT_GENOTYPE="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/Conligen_filtered_withrsid"  # Replace with your actual input genotype file (without extension)
OUTPUT_PREFIX="/hpcfs/users/a1236780/Repos/LiPGS/output/Conligen_filtered_genotype"  # Prefix for filtered output files

# Set the number of threads based on --cpus-per-task above
# THREADS=$SLURM_CPUS_PER_TASK

plink --version
echo "Input file: $INPUT_GENOTYPE"
# echo "Output prefix: $OUTPUT_PREFIX"

# Run PLINK with the specified quality control filters
# plink --bfile $INPUT_GENOTYPE --geno 0.05 --maf 0.01 --hwe 1e-6 --mind 0.05 --allow-no-sex --make-bed --out ${OUTPUT_PREFIX}_qc         
# plink \
#   --bfile $INPUT_GENOTYPE \
#   --geno 0.05 \                    # Exclude SNPs with >5% missing data (i.e., <95% genotyping rate)
#   --maf 0.01 \                     # Exclude SNPs with MAF < 1%
#   --hwe 1e-6 \                      # Exclude SNPs deviating from Hardy-Weinberg Equilibrium
#   --mind 0.05 \                    # Exclude individuals with >5% missing data (i.e., <95% genotype rate)
#   --check-sex \                     # Check for sex inconsistencies
#   --make-bed \                      # Create a new binary PLINK file
#   --out ${OUTPUT_PREFIX}_qc         # Output file prefix

# echo "Quality control complete. Output saved to ${OUTPUT_PREFIX}_qc."
##############################################################
# Run PLINK for GWAS
OUTPUT_DIR="/hpcfs/users/a1236780/Repos/LiPGS/output"  # Output directory for results
Original_pheno="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/ConLiGen_phenotypes_2021_Niguse_withPCs_Country_updated.txt"

European_pheno_binary="/hpcfs/users/a1236780/Repos/LiPGS/output/ConLiGen_european_phenotypes_binary.txt"
# European_pheno_binary="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/binary_outcome.txt"
covariates="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/covariates_for_adjustment.txt"
# Filter for records with Origin equal to European
awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 || $8 == "European" {print $1, $2, $3+1}' "$Original_pheno" > "$European_pheno_binary"

# European_pheno_AldaTOTAL="/hpcfs/users/a1236780/Repos/LiPGS/output/ConLiGen_european_phenotypes_AldaTOTAL.txt"
# awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 || $8 == "European" {print $1, $2, $4}' "$Original_pheno" > "$European_pheno_AldaTOTAL"

# Step 4: Conduct GWAS for binary response (e.g., lithium treatment response)
plink --bfile $INPUT_GENOTYPE --pheno $European_pheno_binary --covar $covariates --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --logistic --out ${OUTPUT_PREFIX}_gwas_binary

#--adjust --allow-no-sex 
# Step 5: Conduct GWAS for continuous response (e.g., ALDA total score)
# plink --bfile $INPUT_GENOTYPE  --pheno $European_pheno_AldaTOTAL --covar $covariates --covar-name Age,Sex,Chip,zPC1,zPC2,zPC3,zPC4 --allow-no-sex --linear --adjust --all-pheno --out ${OUTPUT_PREFIX}_gwas_continuous

# Deactivate the conda environment
conda deactivate