#!/bin/bash
#makes YRI_plink.frq file
plink --bfile /home/wheelerlab1/Data/Stranger_et_al_pop_eQTLs/HapMap3-genotypes/hapmap3_r2_b36_fwd.consensus.qc.poly --freq --keep YRI_individuals.txt --out YRI_plink
