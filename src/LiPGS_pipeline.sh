#!/bin/bash
#SBATCH --job-name=LiPGS_pipeline         # Job name
#SBATCH --output=job_%j.out        # Output file (%j expands to job ID)
#SBATCH --error=job_%j.err         # Error file (%j expands to job ID)
#SBATCH -p batch                    # partition 
#SBATCH -N 1                        # number of nodes
#SBATCH -n 22                    
#SBATCH --mem=128G                   # Memory required per node (adjust based on your data)
#SBATCH --time=50:00:00             # Time limit (hours:minutes:seconds)

module load Anaconda3/2024.06-1
source activate PGS
N_THREADS=22
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Define paths directly in the script
BFILE="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/Conligen_filtered_withrsid"
PHENO_FILE="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/ConLiGen_phenotypes_2021_Niguse_withPCs_Country_updated.txt"
COVAR_FILE="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/covariates_for_adjustment.txt"
PHENO_TYPE="binary"  # Change to "continuous" if using a continuous trait
REF_DIR="/hpcfs/users/a1236780/Repos/PGS/PRScs/LD_reference/1KG/ldblk_1kg_eur"
PATH_TO_PRScs="/hpcfs/users/a1236780/Repos/LiPGS/src/PRScs.py"
NUM_FOLDS=5  # Number of folds for cross-validation
ROOT_DIR="output_binary"  # Change to "output_alda" for continuous traits

# Create root directory
mkdir -p "$ROOT_DIR"

# Run the Python pipeline
python pipeline.py \
    --bfile "$BFILE" \
    --pheno_file "$PHENO_FILE" \
    --covar_file "$COVAR_FILE" \
    --pheno_type "$PHENO_TYPE" \
    --ref_dir "$REF_DIR" \
    --path_to_prscs "$PATH_TO_PRScs" \
    --num_folds "$NUM_FOLDS" \
    --root_dir "$ROOT_DIR" \
    --threads "$N_THREADS"

echo "Pipeline execution completed."


# Run Snakemake, assuming Snakefile is in the same directory
# snakemake extract_phenotype split_data
# snakemake -F
# rm -rf .snakemake/ extract_phenotype split_data --configfile config.yaml 
# nextflow run main.nf -c nextflow.config
