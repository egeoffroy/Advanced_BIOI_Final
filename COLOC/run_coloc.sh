#!/bin/bash
time python3 summary-gwas-imputation/src/run_coloc.py \
-gwas_mode bse \
-gwas coloc/GWAS_YRI_WBC.txt.gz \
-eqtl_mode bse \
-eqtl coloc/eQTL_YRI_WBC.txt.gz \
-gwas_sample_size FROM_GWAS \
-eqtl_sample_size 107 \
-parsimony 8 \
-output coloc/output_YRI_WBC_coloc.txt.gz  > /dev/null
