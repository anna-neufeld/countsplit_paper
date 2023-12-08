## define packages to install
packages <- c("tidyverse", "patchwork")

## install all packages that are not already installed
user_name = "wenhaop"
lib_dir = paste("/home/users/", user_name, "/R_lib", sep="")
install.packages(
    setdiff(packages, rownames(installed.packages(lib_dir))),
    repos = "http://cran.us.r-project.org",
    lib=lib_dir
) # install packages only once to avoid errors when parallel computing

library(tidyverse, lib="/home/users/wenhaop/R_lib")
library(patchwork, lib="/home/users/wenhaop/R_lib")
setwd("~/countsplit_paper/Fig3_mainsims_resubmit/")

## simulation name
simname <- "big_test"

### Make sure that countsplit_paper repository is saved in your home directory. 

#### DATA FOR FIGS 2-4
setwd("~/countsplit_paper/Fig3_mainsims_resubmit/res")
file_names <- dir("~/countsplit_paper/Fig3_mainsims_resubmit/res/", pattern=paste(simname, "*", sep="")) 
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(res) <- c("j", "trueCoeff", "intercept", "pval", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

setwd("~/countsplit_paper/Fig3_mainsims_resubmit/clusterres")
file_names <- dir("~/countsplit_paper/Fig3_mainsims_resubmit/clusterres", pattern=paste(simname, "*", sep="")) 
cluster_res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(cluster_res) <- c("j", "trueCoeff", "intercept", "pval", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

res <- res %>% mutate(intercept2=round(exp(intercept)))
cluster_res <- cluster_res %>% mutate(intercept2=round(exp(intercept)))


null_res <- res %>% filter(trueCoeff==0, prop1==0.5)
power_res <- res %>% filter(j <= 250)
null_cluster_res <- cluster_res %>% filter(trueCoeff==0, prop1==0.5)
power_cluster_res <- cluster_res %>% filter(j <= 250)

set.seed(1)
subsamp <- sample(1:NROW(null_res), size=10^5)
null_res_subset <- null_res[subsamp,]
subsamp <- sample(1:NROW(null_cluster_res), size=10^5)
null_clusterres_subset <- null_cluster_res[subsamp,]

#null_res_subset$intercept3 = "exp(\u03B20)=25"
null_res_subset$intercept3 = "High Intercept"
null_res_subset$intercept3[null_res_subset$intercept2==3] =  "Low Intercept"
null_clusterres_subset$intercept3 = "High Intercept"
null_clusterres_subset$intercept3[null_clusterres_subset$intercept2==3] = "Low Intercept"

null_res_subset$intercept3 = ordered(null_res_subset$intercept3, levels=c("Low Intercept",
                                                                          "High Intercept"))
null_clusterres_subset$intercept3 = ordered(null_clusterres_subset$intercept3, levels=c("Low Intercept",
                                                                          "High Intercept"))

setwd("~/countsplit_paper/Fig3_mainsims_resubmit/")
###### TYPE 1 ERROR FIGURE
p1null <- ggplot(data=null_res_subset, aes(sample=as.numeric(pval), col=as.factor(eps), group=eps))+
  geom_qq(distribution=stats::qunif)+
  facet_grid(col=vars(intercept3))+
  geom_abline(slope=1, col="red")+
  ggtitle("Type 1 error")+
  labs(col=expression(epsilon))+
  xlab("Unif(0,1) Quantiles")+
  ylab("Sample Quantiles")+
  coord_fixed()+
  theme_bw()+
  guides(col="none")
p2null <- ggplot(data=null_clusterres_subset, aes(sample=as.numeric(pval), col=as.factor(eps), group=eps))+
  geom_qq(distribution=stats::qunif)+
  facet_grid(col=vars(intercept3))+
  geom_abline(slope=1, col="red")+
  ggtitle("Type 1 error")+
  xlab("Unif(0,1) Quantiles")+
  ylab("Sample Quantiles")+
  labs(col=expression(epsilon))+
  coord_fixed()+
  theme_bw()+
  guides(col="none")
### Takeaways: Type 1 error is controlled! Epsilon does not matter at ALL for this. 
# p1+p2+plot_layout(nrow=1, guides="collect")
p1null+p2null+plot_layout(nrow=1, guides="collect")
ggsave(paste("figures/",simname, "_Fig2.png", sep=""))

##### DETECTION FIGURE
detection_res <- power_res %>% filter(j==1) %>% group_by(trueCoeff, eps,prop1) %>%
  summarize(avcor = mean(abs(cor)))

detection_res$intercept_dist = "0% Low Intercept"
detection_res$intercept_dist[detection_res$prop1 == 0.5] = "50% Low Intercept"
detection_res$intercept_dist[detection_res$prop1 == 1] = "100% Low Intercept"


detection_res$intercept_dist <- ordered(detection_res$intercept_dist ,
                                        levels=c("0% Low Intercept",
                                                 "50% Low Intercept",
                                              "100% Low Intercept"))


detection_res_cluster <- power_cluster_res %>% group_by(trueCoeff, eps, 
                                                        prop1) %>% summarize(avcor = mean((cor)))
detection_res_cluster$intercept_dist = "0% Low Intercept"
detection_res_cluster$intercept_dist[detection_res$prop1 == 0.5] = "50% Low Intercept"
detection_res_cluster$intercept_dist[detection_res$prop1 == 1] = "100% Low Intercept"

detection_res_cluster$intercept_dist <- ordered(detection_res_cluster$intercept_dist ,
                                        levels=c("0% Low Intercept",
                                                 "50% Low Intercept",
                                                 "100% Low Intercept"))



### DETECTION FIGURE
p1detect <- ggplot(data=detection_res, aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Detection")+
  ylab(expression('Correlation between '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
detection_res_cluster$avcor[detection_res_cluster$avcor < 0] <- 0
p2detect <- ggplot(data=detection_res_cluster, aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Detection")+
  ylab(expression('Adjusted Rand Index between  '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
p1detect+p2detect+plot_layout(guides="collect", nrow=2)
ggsave(paste("figures/",simname, "_Fig3.png", sep=""))



###### POWER
power_res <- power_res %>% filter(prop1==0.5)
power_cluster_res <- power_cluster_res %>% filter(prop1==0.5)
power_res$intercept3 = "High Intercept"
power_res$intercept3[power_res$intercept2==3] =  "Low Intercept"
power_cluster_res$intercept3 = "High Intercept"
power_cluster_res$intercept3[power_cluster_res$intercept2==3] = "Low Intercept"

power_res$intercept3 = ordered(power_res$intercept3, levels=c("Low Intercept",
                                                                          "High Intercept"))
power_cluster_res$intercept3 = ordered(power_cluster_res$intercept3, levels=c("Low Intercept",
                                                             "High Intercept"))

                                                                                        
p1power <- ggplot(data=power_res, aes(x=abs(fitCoeff), y=as.numeric(pval<0.05), col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm", method.args=list(family="binomial"))+
  facet_grid(cols=vars(intercept3))+
  xlab(expression(beta(widehat(L)(X^{train}), bold(X)^{test}))) + ylab("Power")+ggtitle("Power")+
  labs(col=expression(epsilon))+
  theme_bw()

p2power <- ggplot(data=power_cluster_res, aes(x=abs(fitCoeff), y=as.numeric(pval<0.05), col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm", method.args=list(family="binomial"))+
  facet_grid(cols=vars(intercept3))+
  xlab(expression(beta(widehat(L)(X^{train}), bold(X)^{test}))) + ylab("Power")+ggtitle("Power")+
  labs(col=expression(epsilon))+
  theme_bw()

p1power+p2power+plot_layout(guides="collect", nrow=1)
ggsave(paste("figures/",simname, "_Fig4.png", sep=""))


p1null+ p1detect+ p1power+
  p2null+p2detect+p2power+
  plot_layout(guides="collect", nrow=2)
ggsave(paste("figures/",simname, "_megaFig2.eps", sep=""), width=15, height=6.5)
