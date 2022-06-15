setwd("~/countsplit_paper/Fig1_globalnull/")
load("global_null_res_may_23.RData")
load("pvals_DE_pt.RData")
library(tidyverse)

pvals_pseudotime$ID <- 1:NROW(pvals_pseudotime)
pvals_DE_pt$ID <- 1:NROW(pvals_DE_pt)

test <- full_join(pvals_pseudotime, pvals_DE_pt)



naivecol <- "#E763F3"
cellsplitcol <- "#00BF7D"
countsplitcol <- "#E7861B"
jackstrawcol <- "darkblue"
DEcol <- "cyan"

pall <- ggplot(data=test)+
  geom_qq(aes(sample=Naive, col="Naive"), distribution = "qunif")+
  geom_qq(aes(sample=CellSplitting, col="Cell Splitting"), distribution = "qunif")+
  geom_qq(aes(sample=Jackstraw, col="Jackstraw"), distribution = "qunif")+
  geom_qq(aes(sample=DE, col="PseudotimeDE"), distribution="qunif")+
  geom_qq(aes(sample=CountSplitting, col="Count Splitting, \u03B5 =0.5"), distribution = "qunif")+
  ggtitle("All")+
  coord_fixed()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, col="red")+
  #xlab("Unif(0,1) Quantiles")+
  ylab("Sample Quantiles")+
  scale_color_manual(
    limits = c("Naive", "Cell Splitting", "Jackstraw", "PseudotimeDE", "Count Splitting, \u03B5 =0.5"),
    values=c(naivecol, cellsplitcol, jackstrawcol,DEcol, countsplitcol))+
  #ggtitle("P-values for each gene across 2000 realizations")+
  labs(col="Method")+xlab(" ")
  
plow <- ggplot(data=test %>% filter(Intercept=="1"))+
    geom_qq(aes(sample=Naive, col="Naive"), distribution = "qunif")+
    geom_qq(aes(sample=CellSplitting, col="Cell Splitting"), distribution = "qunif")+
    geom_qq(aes(sample=Jackstraw, col="Jackstraw"), distribution = "qunif")+
    geom_qq(aes(sample=DE, col="PseudotimeDE"), distribution="qunif")+
    geom_qq(aes(sample=CountSplitting, col="Count Splitting, \u03B5 =0.5"), distribution = "qunif")+
    ggtitle( expression(Lambda[ij]~ "="~1))+
  coord_fixed()+
    theme_bw()+
    geom_abline(slope=1, intercept=0, col="red")+
    xlab("Unif(0,1) Quantiles")+#ylab("Sample Quantiles")+
    scale_color_manual(
      limits = c("Naive", "Cell Splitting", "Jackstraw", "PseudotimeDE", "Count Splitting, \u03B5 =0.5"),
      values=c(naivecol, cellsplitcol, jackstrawcol,DEcol, countsplitcol))+
    #ggtitle("P-values for each gene across 2000 realizations")+
    labs(col="Method")
  
phigh <- ggplot(data=test %>% filter(Intercept=="10"))+
    geom_qq(aes(sample=Naive, col="Naive"), distribution = "qunif")+
    geom_qq(aes(sample=CellSplitting, col="Cell Splitting"), distribution = "qunif")+
    geom_qq(aes(sample=Jackstraw, col="Jackstraw"), distribution = "qunif")+
    geom_qq(aes(sample=DE, col="PseudotimeDE"), distribution="qunif")+
    geom_qq(aes(sample=CountSplitting, col="Count Splitting, \u03B5 =0.5"), distribution = "qunif")+
  ggtitle( expression(Lambda[ij]~ "="~10))+
  coord_fixed()+
    theme_bw()+
    geom_abline(slope=1, intercept=0, col="red")+xlab(" ")+
   # xlab("Unif(0,1) Quantiles")+#ylab("Sample Quantiles")+
    scale_color_manual(
      limits = c("Naive", "Cell Splitting", "Jackstraw", "PseudotimeDE", "Count Splitting, \u03B5 =0.5"),
      values=c(naivecol, cellsplitcol, jackstrawcol,DEcol, countsplitcol))+
    #ggtitle("P-values for each gene across 2000 realizations")+
    labs(col="Method")

library(patchwork)
pall+plow+phigh+plot_layout(guides="collect")
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/paper_v15/Figures/Fig1.png",
       width=12, height=3.5, units="in") 
