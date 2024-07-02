suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

plot_lava <- function(lavaresult, lociref, plotout, df=NULL,
   color_combsig="red", color_sig="blue", color_nsig="gray"){
  
  # Make x-axis
  loci <- fread(lociref, header=TRUE) 
  
  loci_pro <- loci %>%
    group_by(CHR) %>%
    summarise(chr_len=max(STOP), chr_pos=(min(START)+max(STOP))/2) %>%
    mutate(chr_begin=cumsum(as.numeric(chr_len))-chr_len) %>%
    mutate(br=chr_begin+chr_pos) 
  
  loci_ref <- loci_pro %>%
    left_join(loci, ., by=c("CHR"="CHR")) %>%
    mutate(loci_pos=(START+STOP)/2+chr_begin) %>%
    select(c("LOC","loci_pos"))
  
  x_lab <- loci_pro %>%
    select(c("CHR","br"))
  
  # Read the result
  data <- read.table(lavaresult, header = TRUE) %>% 
    left_join(., loci_ref, by=c("locus"="LOC")) %>%
    mutate(Log10Psign=-log10(p)*sign(rho)) %>%
    mutate(color = ifelse(p < 0.05, color_sig, color_nsig))
  
  if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
    for (i in 1:nrow(data)) {
      snps_in_range <- df %>% 
        filter(CHR==data$chr[i], BP>=data$start[i], BP<=data$stop[i])
      
      if (nrow(snps_in_range) > 0 & data$color[i] == color_sig) {
        data$color[i] <- color_combsig
      }
    }
  }
  
  color_mapping <- setNames(data$color, data$color)
  
  # Plot
  p <- ggplot(data, aes(x=loci_pos, y=rho ,color = color)) +
    geom_point(size=2, alpha=0.8, stroke=0) +
    geom_segment(aes(x = loci_pos, y = rho, xend = loci_pos, yend = 0)) +
    
    scale_x_continuous(label = x_lab$CHR, breaks= x_lab$br ) +
    scale_color_manual(values = color_mapping) +
    labs(x="Chromosome", y="rho") +
    ggtitle("Local genetic correlation") +
    
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth  = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      # axis.line = element_line(linetype = "solid", color = "black")
    ) 
  
  loci_pro <- arrange(loci_pro,"CHR")
  
  for (i in seq(1,22,2)){
    p <- p +
      annotate("rect",xmin = loci_pro$chr_begin[i], xmax = loci_pro$chr_begin[i]+loci_pro$chr_len[i],
               ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.2)
  }
  
  ggsave(plotout, plot = p, width = 10, height = 4)
  return(p)
}


if (!exists("source_sign")) {
  suppressMessages(library(optparse))
  option_list <- list(
    make_option(c("-l", "--lavaresult"), type = "character", help = "The LAVA result file"),
    make_option(c("-r", "--lociref"), type = "character", help = "The loci file"),
    make_option(c("-c", "--color"), type = "character", help = "2 colors"),
    make_option(c("-o", "--outplot"), type = "character", help = "The output plot"))

  opt <- parse_args(OptionParser(usage = "Rscript plot_lava.R [options]", option_list = option_list))

  lavaresult <- opt$lavaresult
  lociref <- opt$lociref
  color <- opt$color
  outplot <- opt$outplot

  color <- unlist(strsplit(color, ","))

  plot_lava(lavaresult, lociref, outplot, color_sig=color[1],  color_nsig=color[2])
}