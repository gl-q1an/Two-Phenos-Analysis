suppressMessages(library(optparse))
help_x <- "The exposure phenotype."
help_y <- "The outcome phenotype."
help_xf <- "The exposure file path."
help_yf <- "The outcome file path."
help_pk <- "The plink path"
help_r <- "Clump and the LD reference bfile name."
help_o <- "The result and output dir"
help_p <- "The p threshold, default 5e-8"
help_t <- "LDlinkR token number"
help_description <- "The MR analysis"

option_list <- list(
  make_option(c("-x", "--exposure"), type = "character", help = help_x),
  make_option(c("-y", "--outcome"), type = "character", help = help_y),
  make_option(c("-f", "--exposurefile"), type = "character", help = help_xf),
  make_option(c("-F", "--outcomefile"), type = "character", help = help_yf),
  make_option(c("-P", "--plink"), type = "character",help = help_pk),
  make_option(c("-R", "--ref"), type = "character", help = help_r),
  make_option(c("-o", "--outputdir"), type = "character", help = help_o),
  make_option(c("-p", "--pthresh"),  type = "numeric", default = 5e-8,help = help_p),
  make_option(c("-t", "--token"),  type = "character",help = help_t))

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", 
                               option_list = option_list,
                               description=help_description))
 
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ieugwasr))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(LDlinkR))
suppressMessages(library(RadialMR))

exposure <- as.character(opt$exposure)
exposure_file <- as.character(opt$exposurefile)
outcome <- as.character(opt$outcome)
outcome_file <- as.character(opt$outcomefile)

result_dir <- as.character(opt$outputdir)
sigstate <- 0

p_thresh <- as.numeric(opt$pthresh)

plink_path <- as.character(opt$plink)
ld_ref_path <- as.character(opt$ref)

token_num <- as.character(opt$token)

cat("=================================================\n")
cat("                   MR analysis                   \n")
cat("MR analysis config:\n")
cat(paste0("[1]exposure: ",exposure,"\n"))
cat(paste0("[2]exposurefile: ",exposure_file,"\n"))
cat(paste0("[3]outcome: ",outcome,"\n"))
cat(paste0("[4]outcome_file: ",outcome_file,"\n"))
cat(paste0("[5]output_dir: ",result_dir,"\n"))
cat(paste0("[6]p_thresh: ",p_thresh,"\n"))
cat(paste0("[7]plink_path: ",plink_path,"\n"))
cat(paste0("[8]ld_ref_path: ",ld_ref_path,"\n"))
cat("=================================================\n")


Result <- data.frame(
  Trait1 = rep(exposure,6),
  Trait2 = rep(outcome,6),
  p_thresh = rep(p_thresh,6)
)

# 1. Data Prepare (Clump)
exp_dat <- read_exposure_data(
  filename = exposure_file,
  clump = FALSE,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Frq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  min_pval = 1e-200, 
  log_pval = FALSE, 
  chr_col = "CHR",
  pos_col = "BP"
)

# Clumping instruments
exp_dat <- exp_dat %>% 
  rename(rsid = SNP,pval = pval.exposure)

exp_dat_clumped <- ld_clump(
  dat = exp_dat,
  clump_kb = 10000, 
  clump_r2 = 0.001, 
  clump_p = p_thresh,
  plink_bin = plink_path,
  bfile = ld_ref_path
)

exp_dat_clumped <- exp_dat_clumped %>% 
  rename(SNP = rsid,pval.exposure = pval) %>% 
  mutate(exposure=!!exposure)

Result <- Result %>%
  mutate(N_sig_SNPs=length(exp_dat_clumped$SNP))

print(paste0("Number of IVs: ", as.character(length(exp_dat_clumped$SNP))))

# Outcome Data
out_dat <- read_outcome_data(
  filename = outcome_file, 
  snps = exp_dat_clumped$SNP, 
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Frq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  min_pval = 1e-200, 
  log_pval = FALSE, 
  chr_col = "CHR",
  pos_col = "BP"
)

out_dat <- out_dat %>%
  mutate(outcome=!!outcome)

# Identifying & printing exposure instruments missing from outcome GWAS
missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
print("List of IVs missing from outcome GWAS:")
for (i in 1:length(missing_IVs)) {
  print(paste0(missing_IVs[i]))
}

# Replacing missing instruments from outcome GWAS with proxies
if(length(missing_IVs) == 0) {
  print("All exposure IVs found in outcome GWAS.")
} else {
  print("Some exposure IVs missing from outcome GWAS.")
  out_full <- fread(outcome_file)
  for (i in 1:length(missing_IVs)) {
    proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = token_num, file = FALSE)
    proxies <- proxies[proxies$R2 > 0.8, ]
    proxy_present = FALSE
    if(length(proxies$RS_Number) == 0){
      print(paste0("No proxy SNP available for ", missing_IVs[i]))
    } else {
      for (j in 1:length(proxies$RS_Number)) {
        proxy_present <- proxies$RS_Number[j] %in% out_full$SNP
        if (proxy_present) {
          proxy_SNP = proxies$RS_Number[j]
          proxy_SNP_allele_1 = str_sub(proxies$Alleles[j], 2, 2)
          proxy_SNP_allele_2 = str_sub(proxies$Alleles[j], 4, 4)
          original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
          original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
          break
        }
      }
    }
    
    if(proxy_present == TRUE) {
      print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
      proxy_row <- out_dat[1, ]
      proxy_row$SNP = missing_IVs[i]
      proxy_row$beta.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
      proxy_row$se.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
      if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
      if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
      proxy_row$pval.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "P"])
      proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N"])
      proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "chr.exposure"])
      proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "pos.exposure"])
      if("Frq" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "Frq"])
      out_dat <- rbind(out_dat, proxy_row)
    }
    
    if(proxy_present == FALSE) {
      print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
    }
  }
}

cat("Harmonise DATA!\n")

dat <- harmonise_data(
  exposure_dat = exp_dat_clumped, 
  outcome_dat = out_dat, 
  action = 2
)

rm(exp_dat,out_full)

# Calculate the F-value
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("Frq is not found, F-value can't be calculated!")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("BETA is not found, F-value can't be calculated!")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("SE is not found, F-value can't be calculated!")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("SampleSize is not found, F-value can't be calculated!")
    return(dat)}
  
  
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}
dat1 <- get_f(dat)

cat(paste0("SNPs after f-stastic: ",as.character(nrow(dat1)),"\n"))

# Remove outliers
# cat("Remove outliers!\n")

# dat_radial <- format_radial(dat1$`beta.exposure`,dat1$`se.exposure`,dat1$`beta.outcome`,dat1$`se.outcome`,dat1$SNP)
# radial_ivw_obj <- ivw_radial(dat_radial,0.05,1,0.0001)  # Cochran Q
# radial_egger_obj <- egger_radial(dat_radial,0.05,1)   # Rucker's Q

# radial_ivw_dat <- radial_ivw_obj$data %>%
#   filter(Outliers=="Variant")
# radial_egger_dat <- radial_egger_obj$data %>%
#   filter(Outliers=="Variant")

# cat("SNPs qc is done!\n")
dat_qc <- dat1

# dat_qc <- subset(dat1, SNP %in% radial_ivw_dat$SNP & SNP %in% radial_egger_dat$SNP)

# cat(paste0("SNPs after radial qc: ",as.character(nrow(dat_qc)),"\n"))

mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj) 
{
  library(TwoSampleMR)
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
  method_list_modified <- stringr::str_replace_all(method_list, "mr_raps","mr_raps_modified")
  
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),function(x1)
  {
    x <- subset(x1, mr_keep)
    
    if (nrow(x) == 0) {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list_modified, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    }
    )
    
    methl <- mr_method_list()
    mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1], 
                         method = methl$name[match(method_list, methl$obj)], 
                         nsnp = sapply(res, function(x) x$nsnp), 
                         b = sapply(res, function(x) x$b), 
                         se = sapply(res, function(x) x$se), 
                         pval = sapply(res, function(x) x$pval))
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    
    return(mr_tab)
  }
  )
  return(mr_tab)
}

tsmr1<-mr(dat_qc, method_list=c("mr_weighted_median"))
tsmr2<-mr(dat_qc, method_list=c("mr_ivw"))
tsmr3<-mr(dat_qc, method_list=c("mr_egger_regression"))
tsmr4<-mr(dat_qc, method_list=c("mr_weighted_mode"))
tsmr5<-mr_modified(dat_qc, method_list=c("mr_raps"))
tsmr6<-mr(dat_qc, method_list=c("mr_wald_ratio"))

method <- c("mr_weighted_median","mr_ivw","mr_egger_regression","mr_weighted_mode","mr_raps","mr_wald_ratio")

tsmr_list <- list(tsmr1,tsmr2,tsmr3,tsmr4,tsmr5,tsmr6)

lowerCI <- function(beta,df,SE){
  return(beta - (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
}

upperCI <- function(beta,df,SE){
  return(beta + (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
}

CI <- 0.95

for (i in 1:length(tsmr_list)){
  df_i <- tsmr_list[[i]]
  
  n_snp_i <- df_i$nsnp
  beta_i <- df_i$b
  se_i <- df_i$se

  pval_i <- df_i$pval

  if (is.numeric(pval_i)){
    if (pval_i < 0.05){
      sigstate <- 1
  }}

  if (i %in% c(1, 2, 5, 6)) {
    beta_low95_i <- df_i$b-1.96*df_i$se
    beta_up95_i <- df_i$b+1.96*df_i$se
  } else if (i==3) {
    beta_low95_i <- mapply(lowerCI, df_i$b, df_i$nsnp - 2, df_i$se)
    beta_up95_i <- mapply(upperCI, df_i$b, df_i$nsnp - 2, df_i$se)
  } else if (i==4) {
    beta_low95_i <- mapply(lowerCI, df_i$b, df_i$nsnp - 1, df_i$se)
    beta_up95_i <- mapply(upperCI, df_i$b, df_i$nsnp - 1, df_i$se)
  }
  
  Result$method[i] <- method[i]
  Result$nsnp[i] <- ifelse(length(n_snp_i) == 0, NA, n_snp_i)
  Result$beta[i] <- ifelse(length(beta_i) == 0, NA, beta_i)
  Result$se[i] <- ifelse(length(se_i) == 0, NA, se_i)
  Result$beta_low95[i] <- ifelse(length(beta_low95_i) == 0, NA, beta_low95_i)
  Result$beta_up95[i] <- ifelse(length(beta_up95_i) == 0, NA, beta_up95_i)
  Result$pval[i] <- ifelse(length(pval_i) == 0, NA, pval_i)
}  

cat("Result table is built!\n")

if (sigstate==1) {
  result_prefix <- paste("tsmr",as.character(p_thresh),"sig",exposure,outcome,sep="_")
} else {
  result_prefix <- paste("tsmr",as.character(p_thresh),"no",exposure,outcome,sep="_")
}

outprefix <- file.path(result_dir,result_prefix)

save(dat_qc,file = paste0(outprefix,".rdata"))
write.table(Result,paste0(outprefix,".csv"),row.names = F,quote = F,sep = ",")
cat(paste0("MR analysis between ",exposure," and ",outcome, " is done!!\n"))