library(dplyr)
library(data.table)
data <- fread("Wojcik_build37/31217584-GCST008049-EFO_0004308-build37.f.tsv.gz", header = T, stringsAsFactors=F)
data <- data.frame(data$variant_id, data$beta, data$standard_error)
colnames(data) <- c("variant_id", "beta", "standard_error")
data <- transform(data, zscore = beta/standard_error)
#data$zscore <- data$beta/data$standard_error
print(head(data))
data <- data.frame(data$variant_id, data$zscore)
colnames(data) <- c("rsid", "zscore")
write.table(data, "Wojcik_WBC_rs_z.txt",  quote = F, row.names=F)
