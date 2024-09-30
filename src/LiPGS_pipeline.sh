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

# Run Snakemake, assuming Snakefile is in the same directory
snakemake --jobs 10 --latency-wait 60 \
    --cluster "sbatch --partition=batch --cpus-per-task=12 --mem=64G --time=5:00:00" \
    --cluster-status "scontrol show job {jobid}" \
    --keep-going

# Explanation:
# --jobs 10: Allow up to 10 parallel jobs
# --latency-wait 60: Wait 60 seconds before considering a job failure due to filesystem latency
# --cluster "sbatch ...": Specifies how to submit jobs to SLURM for each Snakemake rule
# --cluster-status "scontrol show job {jobid}": Command to check job status in SLURM
# --keep-going: Keep running the workflow despite errors in individual jobs