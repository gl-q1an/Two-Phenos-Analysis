suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

plot_cfdr <- function(cfdrresult, plotout, pheno1, pheno2){

  fdrresult <- fread(cfdrresult, header=TRUE) %>%
    na.omit(.) %>%
    mutate(logFDR=-log10(FDR))

  chr_tot <- fdrresult %>%
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len)
  
  don <- chr_tot %>% 
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(fdrresult, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot)

  axisdf = don %>% 
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


  p <- ggplot(don, aes(x=BPcum, y=logFDR)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("#004e74","grey"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  labs(x="Chromosome", y="-log10(FDR)") +
  ggtitle(paste0(pheno1," & ", pheno2, " Manhattan Plot")) +

  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +

  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
  )

  for (i in seq(1,22,2)){
    p <- p +
      annotate("rect",xmin = chr_tot$tot[i], xmax = chr_tot$tot[i+1],
               ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.2)
  }
  
  ggsave(plotout, plot = p, width = 10, height = 4)
  
  return(p)
}

