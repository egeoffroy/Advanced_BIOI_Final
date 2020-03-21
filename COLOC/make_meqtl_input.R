library(data.table)
data <- fread("/home/wheelerlab3/files_for_revisions_plosgen/meqtl_results/hapmap/next", header = F, stringsAsFactors=F)
colnames(data) <- c("snps", "gene", "statistic", "pvalue", "FDR", "beta")
write.table(data, "/home/elyse/meqtl_YRI_input.txt", col.names=T, row.names=F, quote=F)
