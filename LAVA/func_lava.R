library(optparse)
option_list <- list(
  make_option(c("-c", "--lavaconfig"), type = "character", default = NULL,
              help = "The LAVA info file"),
  make_option(c("-p", "--phe1"), type = "character", default = NULL,
              help = "The Pheno 1"),
  make_option(c("-q", "--phe2"), type = "character", default = NULL,
              help = "The Pheno 2"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "The outfile prefix"),
  make_option(c("-r", "--refdir"), type = "character", default = NULL,
              help = "The reference dir"),
  make_option(c("-m", "--max"), type = "integer", default = NULL,
              help = "max/one task"))

opt <- parse_args(OptionParser(usage = "Rscript func_lava.R [options]", option_list = option_list))

inputinfo <- opt$lavaconfig
phe1 <- opt$phe1
phe2 <- opt$phe2
outfile <- opt$outfile
refdir <- opt$refdir
core_num <- as.numeric(opt$max)

library(LAVA)
library(parallel)

input <- process.input(input.info.file=inputinfo,
					  sample.overlap.file=NULL,
					  ref.prefix=paste0(refdir,"/g1000_eur"),
					  phenos=c(phe1,phe2))

ls(input)
ls(input$sum.stats)
head(input$sum.stats[[phe1]])
head(input$sum.stats[[phe2]])

loci <- read.loci(paste0(refdir,"/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
n.loc <- nrow(loci)

#### ANALYSE ####
print(paste("Starting LAVA analysis for",n.loc,"loci"))
progress <- ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (used for printing the progress)
u=b=NULL

# Using mclapply to analyse the loci in parallel
# adjust the mc.cores argument according to the number of available cores (e.g. N avail cores - 1)
# if you dont want to parallelise you can use lapply or a for loop instead
out <- mclapply(1:n.loc, mc.cores=core_num, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)]))) 
	locus <- process.locus(loci[i,], input)

	if (!is.null(locus)) {
		loc.info <- data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
		ub <- run.univ.bivar(locus, univ.thresh= 2e-5)
		u <- cbind(loc.info, ub$univ)
		if (!is.null(ub$bivar)) {
			b <- cbind(loc.info, ub$bivar)
		}
	}
	return(list(univ=u, bivar=b))
})

#save the env
save.image(file=paste0(outfile,".RData"))

### WRITE OUTPUT ###
write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(outfile,".univ"), row.names=F, quote=F)
write.table(do.call(rbind, lapply(out, "[[", "bivar")), paste0(outfile,".bivar"), row.names=F, quote=F)