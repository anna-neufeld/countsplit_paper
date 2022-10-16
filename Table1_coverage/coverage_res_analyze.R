library(tidyverse)
library(patchwork)
setwd("~/countsplit_paper/Table1_coverage/res")

setwd("~/countsplit_paper/Table1_coverage/res")
file_names <- dir("~/countsplit_paper/Table1_coverage/res", pattern="april_13.*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(res) <- c("j", "trueCoeff", "intercept", "pval", "coverage1" , "coverage2", "coverage3", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

setwd("~/Dropbox/Pseudotime : PCA NEW/Paper/FinalCode/Coverage/clusterres")
file_names <- dir("~/Dropbox/Pseudotime : PCA NEW/Paper/FinalCode/Coverage/clusterres", pattern="april_13.*")
cluster_res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(cluster_res) <- c("j", "trueCoeff", "intercept", "pval","coverage1", "coverage2", "coverage3", "fitCoeff", "cor", "eps", "type", "propImp",
                        "n", "p", "prop1")


res <- res %>% mutate(type = "pt")
cluster_res <- res %>% mutate(type = "cluster")
res <- rbind(res, cluster_res)

as.data.frame(res %>% filter(n==200) %>% 
  group_by(type, exp(intercept), j<= 10) %>% summarize(mean(coverage2)))



cons_res <- res %>%
  group_by(j, trueCoeff, intercept, eps, type, n) %>%
  summarize(cov1 = mean(coverage1), 
            cov2 = mean(coverage2), 
            cov3 = mean(coverage3),
            num = n())

