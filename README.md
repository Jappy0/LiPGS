# LiPGS

### A pipeline to generate Polygenic scores for predicting Lithium treatment response (LiPGS) in patients with Bipolar Disorder (BD)

**Important Note:**  Please note that this is not the official code from the original research. This pipeline is based on quality-controlled and imputed genotype data. For details on quality control, imputation, and original results, refer to the original study.

Amare, A.T., Thalamuthu, A., Schubert, K.O. et al. Association of polygenic score and the involvement of cholinergic and glutamatergic pathways with lithium treatment response in patients with bipolar disorder. Mol Psychiatry 28, 5251–5261 (2023). https://doi.org/10.1038/s41380-023-02149-1

## Preparation

### Clone this repo
```
git clone https://github.com/JappyPing/LiPGS.git
```

### Create and activate virtual environment for easy deployment
```
conda create -n LiPGS python
```
```
conda activate LiPGS 
```
### Install dependencies
```
conda install bioconda::plink2
```
```
pip install mpire
```

## Usage
```
usage: pipeline.py [-h] --bfile BFILE --pheno_file PHENO_FILE --pheno_name PHENO_NAME --pheno_type
                   {binary,alda} --covar_file COVAR_FILE --ref_dir REF_DIR --src_path SRC_PATH
                   [--num_folds NUM_FOLDS] [--root_dir ROOT_DIR] [--fold_data_dir FOLD_DATA_DIR]
                   [--gwas_dir GWAS_DIR] [--prscs_dir PRSCS_DIR] --model MODEL [--threads THREADS]

Run the GWAS and PRS pipeline

options:
  -h, --help            show this help message and exit
  --bfile BFILE         Path to PLINK bfile prefix
  --pheno_file PHENO_FILE
                        Path to phenotype file
  --pheno_name PHENO_NAME
                        Phenotype name in the file of --pheno_file
  --pheno_type {binary,alda}
                        Type of phenotype
  --covar_file COVAR_FILE
                        Path to covariates file
  --ref_dir REF_DIR     Path to PRScs reference directory
  --src_path SRC_PATH   Path to PRScs.py script
  --num_folds NUM_FOLDS
                        Number of folds for cross-validation
  --root_dir ROOT_DIR   Root output directory
  --fold_data_dir FOLD_DATA_DIR
                        Splited data output directory
  --gwas_dir GWAS_DIR   GWAS output directory
  --prscs_dir PRSCS_DIR
                        PRScs output directory
  --model MODEL         Model used for gwas.
  --threads THREADS     Threads or processes number for PLINK2 or PRScs
```

## Note
As this pipeline utilize [PRScs](https://github.com/getian107/PRScs?tab=readme-ov-file) to generate PRSs, please include the following configuration to your job submission script to specify the threads (e.g., ```N_THREADS=22```) used in scipy if you are running this pipeline on a shared computer cluster. Refer [PRScs](https://github.com/getian107/PRScs?tab=readme-ov-file) for details.

```
N_THREADS=22
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS
```

## Example
```
python ./src/pipeline.py --bfile /path/to/your_bfile_prefix 
                    --pheno_file /path/to/your_phenotype_file 
                    --pheno_name Response 
                    --covar_file /path/to/your_covariates_file 
                    --pheno_type binary 
                    --ref_dir /path/to/your_LD_reference 
                    --src_path /path/to/src 
                    --num_folds 5 
                    --root_dir /results 
                    --threads 22 
                    --model logistic.hybrid
```

## Support
Feel free to contact me if you have any questions and found any bugs.