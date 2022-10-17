library(tidyverse)
library(patchwork)
setwd("~/countsplit_paper/Table1_coverage/res")

setwd("~/countsplit_paper/Table1_coverage/res")
file_names <- dir("~/countsplit_paper/Table1_coverage/res", pattern="oct16.*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(res) <- c("j", "trueCoeff", "intercept", "pval", "coverage2", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

setwd("~/countsplit_paper/Table1_coverage/clusterres")
file_names <- dir("~/countsplit_paper/Table1_coverage/clusterres", pattern="oct16.*")
cluster_res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(cluster_res) <- c("j", "trueCoeff", "intercept", "pval", "coverage2","fitCoeff", "cor", "eps", "type", "propImp",
                        "n", "p", "prop1")


res <- res %>% mutate(type = "pt")
cluster_res <- cluster_res %>% mutate(type = "cluster")
res <- rbind(res, cluster_res)

as.data.frame(res %>%
  group_by(type, exp(intercept), trueCoeff==0) %>% summarize(mean(coverage2)))
