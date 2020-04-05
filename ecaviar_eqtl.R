library(dplyr)
library(data.table)

eqtl <- fread('coloc/Sig_Height_eQTL_YRI_Height.txt.gz', header = T, stringsAsFactors=F)
eqtl <- unique(eqtl)
for(i in c(1:22)){
  snps_1 <- fread(paste('snps_', i, '.txt', sep = ''), header = F, stringsAsFactors=F)
  eqtl1 <- eqtl %>% filter(eqtl$variant_id %in% snps_1$V1)
  eqtl1$Z <- as.numeric(eqtl1$slope)/as.numeric(eqtl1$slope_se)
  eqtl2 <- data.frame(eqtl1$variant_id, eqtl1$Z)
  eqtl2 <- unique(eqtl2)
  eqtl2 <- eqtl2[!duplicated(eqtl2$eqtl1.variant_id),]
  write.table(eqtl2, paste('caviar/eqtl_chr', i '_one_val.txt', sep =''), quote = F, row.names=F, col.names=F)
}
