#!/bin/bash
#SBATCH --job-name=PRScs         # Job name
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

# python PRScs.py --ref_dir=/hpcfs/users/a1236780/Repos/PGS/PRScs/LD_reference/1KG/ldblk_1kg_eur --bim_prefix=/hpcfs/users/a1236780/Repos/LiPGS/output/Conligen_filtered_genotype_gwas_binary --sst_file=/hpcfs/users/a1236780/Repos/LiPGS/output/Conligen_filtered_genotype_gwas_binary.assoc.logistic --n_gwas=4652947 --phi=1e-2 --out_dir=./eur

python PRScs.py --ref_dir=/hpcfs/users/a1236780/Repos/PGS/PRScs/LD_reference/1KG/ldblk_1kg_eur --bim_prefix=/hpcfs/users/a1236780/Repos/PGS/PRScs/test_data/test --sst_file=/hpcfs/users/a1236780/Repos/PGS/PRScs/test_data/sumstats_se.txt --n_gwas=200000 --chrom=22 --phi=1e-2 --out_dir=./test

INPUT_GENOTYPE="/hpcfs/users/a1236780/Repos/LiPGS/ConLiGen/Conligen_filtered_withrsid"
plink --bfile $INPUT_GENOTYPE --score test_pst_eff_a1_b0.5_phi1e-02_chr22.txt 2 4 6 --out prs_output