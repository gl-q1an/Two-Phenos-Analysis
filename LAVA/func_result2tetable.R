args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript lava_resulttable.R univfile bivarfile output")
} else {
  univfile <- args[1]
  bivarfile <- args[2]
  output <- args[3]
}

univ <- read.table(univfile, header = TRUE)
bivar <- read.table(bivarfile, header = TRUE)

if (nrow(bivar) == 0) {
  warning1 <- paste0(bivarfile," with 0 line, there may not be significant results, Please check!")
  cat(warning1, file = output)
  stop(warning1)
}

bivar$p_fdr <- p.adjust(bivar$p, method='BH')
bivar$p_bfn <- bivar$p*nrow(bivar)
bivar$p_bfn[bivar$p_bfn>1] <- 1

bivar_sig <- subset(bivar, p < 0.05)

if (nrow(bivar_sig) == 0) {
  warning2 <- paste0(bivarfile," p<0.05 with 0 line, there may not be significant results, Please check!")
  cat(warning2, file = output)
  stop(warning2)
}

library(tidyr)

univ <- univ[c(1,7:9)]
univ_w <- pivot_wider(univ, names_from = phen, values_from = c(h2.obs, p))
phe1 <- bivar_sig[1,7]
phe2 <- bivar_sig[1,8]
univ_colname <- names(univ_w)
univ_colname <- gsub(phe1, "trait1", univ_colname)
univ_colname <- gsub(phe2, "trait2", univ_colname)
colnames(univ_w) <- univ_colname

result <- merge(univ_w,bivar_sig,by="locus")
colnames <- c("phen1","phen2","locus","chr","start","stop","n.snps","n.pcs",
              "h2.obs_trait1","p_trait1","h2.obs_trait2","p_trait2",
              "rho","rho.lower","rho.upper","r2","r2.lower","r2.upper","p","p_fdr","p_bfn")
result <- result[,colnames]
colnames(result)[1:2] <- c("trait1","trait2")

write.table(result,output,row.names = FALSE,quote = FALSE)

cat(paste0(phe1," & ",phe2, " is finished!\n"))
