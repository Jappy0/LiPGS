#!/usr/bin/env nextflow

// Set the root directory based on the phenotype type (binary or continuous)
root_dir = params.plink_pheno_type == 'binary' ? 'output_binary' : 'output_alda'

// Step 1: Create the root directory if it does not exist
process createRootDir {
    output:
    dir root_dir

    script:
    """
    mkdir -p ${root_dir}
    """
}

// Step 2: Extract the phenotype according to the phenotype type (binary or continuous)
process extractPhenotype {
    input:
    file pheno_file from params.plink_pheno_file

    output:
    file("${root_dir}/European_phenotype.txt")

    script:
    if (params.plink_pheno_type == 'binary') {
        """
        awk 'BEGIN {FS="\\t"; OFS="\\t"} NR==1 || \$8 == "European" {print \$1, \$2, \$3+1}' ${pheno_file} > ${root_dir}/European_phenotype.txt
        """
    } else {
        """
        awk 'BEGIN {FS="\\t"; OFS="\\t"} NR==1 || \$8 == "European" {print \$1, \$2, \$4}' ${pheno_file} > ${root_dir}/European_phenotype.txt
        """
    }
}

// Step 3: Use `split_data.py` to generate datasets for the 5-fold cross-validation
process splitData {
    input:
    file phenotype_file from extractPhenotype.out

    output:
    file("${root_dir}/fold_data/discovery_*")

    script:
    """
    python split_data.py --phenotype_file ${phenotype_file} --plink_prefix ${params.plink_bfile} --output_dir ${root_dir}/fold_data --num_groups ${params.n_folds}
    """
}

// Step 4: Iteratively run GWAS for each fold using discovery data
process runGwas {
    input:
    file bed_f from "${root_dir}/fold_data/discovery_${fold}.bed"
    file bim_f from "${root_dir}/fold_data/discovery_${fold}.bim"
    file fam_f from "${root_dir}/fold_data/discovery_${fold}.fam"
    file pheno from "${root_dir}/fold_data/discovery_${fold}.txt"
    file covar from params.plink_covar_file

    output:
    file("${root_dir}/gwas_fold_${fold}.assoc.${params.plink_pheno_type}")
    file("${root_dir}/gwas_fold_${fold}.log")

    script:
    def analysis_type = params.plink_pheno_type == 'binary' ? '--logistic' : '--linear'

    """
    plink --bed ${bed_f} \
          --bim ${bim_f} \
          --fam ${fam_f} \
          --pheno ${pheno} \
          --covar ${covar} \
          --covar-name ${params.plink_covar_name} \
          --allow-no-sex \
          --keep-allele-order \
          --real-ref-alleles \
          ${analysis_type} \
          --adjust \
          --all-pheno \
          --ci 0.95 \
          --make-bed \
          --out ${root_dir}/gwas_fold_${fold}
    """
}

// Step 5: Iteratively calculate PRS using PRScs with GWAS results
process runPrscs {
    input:
    file ref_dir from params.prscs_ref_dir
    file summary_stats from runGwas.out.collect { f -> f.name.contains(".assoc") }
    file bed_f from "${root_dir}/fold_data/target_${fold}.bed"
    file bim_f from "${root_dir}/fold_data/target_${fold}.bim"
    file fam_f from "${root_dir}/fold_data/target_${fold}.fam"

    output:
    file("${root_dir}/prs_fold_${fold}")
    file("${root_dir}/plink_prs_fold_${fold}.txt")

    script:
    """
    # Step 1: Add A2 to summary stats
    awk 'BEGIN {FS="\\t"; OFS="\\t"} {if (NR==1) print \$0, "A2"; else print \$0, \$5;}' ${summary_stats} > ${root_dir}/prs_fold_${fold}.a2

    # Step 2: Run PRScs to calculate PRS
    python ${params.prscs_path} \
        --ref_dir ${ref_dir} \
        --bim_file ${bim_f} \
        --sst_file ${root_dir}/prs_fold_${fold}.a2 \
        --n_gwas ${params.n_gwas} \
        --out ${root_dir}/prs_fold_${fold}

    # Step 3: Run PLINK to score the PRS
    plink --bed ${bed_f} \
          --bim ${bim_f} \
          --fam ${fam_f} \
          --score ${root_dir}/prs_fold_${fold} 2 4 6 \
          --out ${root_dir}/plink_prs_fold_${fold}
    """
}

// Step 6: Combine all PRS scores into a single output
process combinePrs {
    input:
    file prs_files from runPrscs.out.collect { f -> f.name.contains("plink_prs") }

    output:
    file("${root_dir}/combined_prs.txt")

    script:
    """
    cat ${prs_files.join(' ')} > ${root_dir}/combined_prs.txt
    """
}

// Workflow execution order
workflow {
    createRootDir()

    extractPhenotype()

    splitData()

    runGwas()

    runPrscs()

    combinePrs()
}
