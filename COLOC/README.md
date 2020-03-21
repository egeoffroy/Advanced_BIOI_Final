## Simulated Dataset
To run COLOC with a simulated dataset of 50 SNPs and 200 samples:

```
Rscript coloc_sim_data.R
```

## YRI and Wojcik GWAS Summary Statistics
To run COLOC with the International HapMap Project Yoruba cohort and the Wojcik White Blood Cell Count GWAS Summary Statistics, run the appropriate 01 scripts necessary to format your input data. Then run:

```
Rscript 02_make_coloc_YRI.R
time bash 03_run_coloc.sh
```
