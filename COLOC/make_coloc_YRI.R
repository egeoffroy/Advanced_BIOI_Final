library(data.table)
library(dplyr)
library(R.utils)

"%&%" = function(a,b) paste(a,b,sep="")
phenos <- c("WBC")
chrs <- c(1:22)
pops <- c("YRI") #do combined pops later
pops_sample_size <- c(107) #R doesn't have dicts so we're doing it a slgihtly more ratchet way
sig_gene_SNPs <- fread("/home/elyse/snps.txt", header = F) #so we don't run all the SNPs b/c it takes forever
sig_gene_SNPs <- sig_gene_SNPs$V1

for(pop in 1:length(pops)){ #read in pop's .frq file for MAF
  frq <- fread("/home/elyse/YRI_plink.frq")
  frq <- frq %>% dplyr::select(SNP, MAF)

  for(pheno in phenos){ #read in GEMMA output file
    GEMMA_result <- fread("/home/elyse/Wojcik_build37/31217584-GCST008039-EFO_0004309-build37.f.tsv.gz", header = T)
    GEMMA_result$chr_pos <- paste(gsub("chr", "", GEMMA_result$chromosome), GEMMA_result$base_pair_location, sep = ":")
    GEMMA_for_COLOC <- GEMMA_result %>% dplyr::select(variant_id, beta, standard_error, effect_allele_frequency) #subset to COLOC input
    GEMMA_for_COLOC$sample_size <- 28608
    colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
    GEMMA_for_COLOC <- GEMMA_for_COLOC[complete.cases(GEMMA_for_COLOC),] #COLOC does not like missing values
    #GWAS_write <- data.frame(panel_variant_id = character(), effect_size = numeric(), standard_error = numeric(), frequency = numeric(), sample_size = numeric(), stringsAsFactors = F)
    GWAS_write <- GEMMA_for_COLOC
    eQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F)

    for(chr in chrs){ #yes triple loops are ratchet
      #system("zcat -f /home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/" %&% pops[pop] %&% "_Nk_10_PFs_chr" %&% chr %&% "pcs_3.meqtl.cis.* > /home/elyse/COLOC/meQTL_input.txt") #fread doesn't seem to like wildcards so we're gonna do this the ugly way
      meqtl <- fread("/home/elyse/meqtl_YRI_input.txt", nThread = 40) #read in matrix eQTL results
      meqtl$se <- meqtl$beta / meqtl$statistic #make your own standard error since it's not in the meQTL output
      meqtl$n_samples <- pops_sample_size[pop]
      meQTL_for_COLOC <- left_join(meqtl, frq, by = c("snps" = "SNP")) #add freq to COLOC input
      meQTL_for_COLOC <- meQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) #subset to COLOC input
      colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
      meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]

      #GWAS_write <- rbind(GWAS_write, GEMMA_for_COLOC_chr)
      eQTL_write <- rbind(eQTL_write, meQTL_for_COLOC)
    }

    snps_in_both <- intersect(GWAS_write$panel_variant_id, eQTL_write$variant_id) #is there a better way to do this? Probably. Do I feel like figuring it out? Nah.
    snps_in_all <- intersect(snps_in_both, sig_gene_SNPs)
    GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_all)
    eQTL_write <- subset(eQTL_write, variant_id %in% snps_in_all)
    #GWAS_write <- GWAS_write[order(GWAS_write$gene_id),] #don't order a column that doesn't exist
    eQTL_write <- eQTL_write[order(eQTL_write$gene_id),]

    fwrite(eQTL_write, "/home/elyse/coloc/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    gzip("/home/elyse/coloc/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", destname = "/home/elyse/coloc/eQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
    fwrite(GWAS_write, "/home/elyse/coloc/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    gzip("/home/elyse/coloc/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", "/home/elyse/coloc/GWAS_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz")
    print("Completed with " %&% pops[pop] %&% ", for " %&% pheno %&% ".")
  }
}
