setwd("~/countsplit_paper/Fig2_overdisp")
setwd("./overdispres_pt")
file_names <- dir(".", pattern="Feb2*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(res) <- c("j", "trueCoeff", "intercept", "pval", "fitCoeff", "cor", "eps",
                "n", "p", "b")
library(tidyverse)
library(patchwork)
res <- res %>% mutate(intercept = exp(intercept)) %>%
  mutate(overdisp1 = round(intercept/b,2))

naivecol <- "#E763F3"
countsplitcol <- "#E7861B"

res$method <- ifelse(res$eps==0.5, "Count splitting, \u03B5 =0.5", "Naive")
res$overdisp2 <- paste("\u039B/b =", res$overdisp1)


### New plot- only the Lambda = 5 panels. 
ggplot(data=res %>% 
         filter(pval <= 1, overdisp1 != 0.25, overdisp1 != 5, intercept > 4, intercept < 5.1),
       aes(sample=as.numeric(pval), col=method, group=method))+
  geom_qq(distribution="qunif")+
  facet_grid(cols=vars(overdisp2))+
  geom_abline(slope=1, col="red")+ggtitle("P-values under a negative binomial model")+
  labs(col="Method")+
  coord_fixed()+
theme_bw()+scale_color_manual(
  limits=c("Count splitting, \u03B5 =0.5", "Naive"),
  values=c(countsplitcol, naivecol))+
  xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")
#ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/paper_v22/Figures/overdisp_new.png",
#       width=13, height=4)
