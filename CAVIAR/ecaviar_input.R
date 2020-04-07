library(data.table)
library(dplyr)
chr <- c(1:22)

gwas <- fread("Wojcik_build37/31217584-GCST008053-EFO_0004339-build37.f.tsv.gz", header=T, stringsAsFactors=F)
eqtl <- fread("eQTL_ecaviar_Z.txt", header=F, stringsAsFactors=F)
for(i in chr){
  gwas1 <- gwas %>% filter(chromosome == i)
  snps <- gwas1$variant_id
  gwas1 <- data.frame(gwas1$panel_variant_id, gwas1$standard_error, gwas1$effect_size)
  colnames(gwas1) <- c("variant_id", "se", "effect_size")
  gwas1$Z <- as.numeric(gwas1$effect_size)/as.numeric(gwas1$se)
  gwas1 <- data.frame(gwas1$variant_id, gwas1$Z)

  eqtl1 <- subset(eqtl, eqtl$V1 == snps)
  write.table(gwas1, paste('gwas_', i, '_height.txt', sep= ''), quote=F, row.names=F, col.names=F)
  write.table(eqtl1, paste('eqtl_YRI_', i, '_height.txt', sep= ''), quote=F, row.names=F, col.names=F)
}
