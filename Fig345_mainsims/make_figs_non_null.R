library(tidyverse)
library(patchwork)
setwd("~/countsplit_paper/Fig345_mainsims")

### Make sure that countsplit_paper repository is saved in your home directory. 

#### DATA FOR FIGS 2-4
setwd("~/countsplit_paper/Fig345_mainsims/res")
file_names <- dir("~/countsplit_paper/Fig345_mainsims/res", pattern="jan11sim*") 
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(res) <- c("j", "trueCoeff", "intercept", "pval", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

setwd("~/countsplit_paper/Fig345_mainsims/clusterres")
file_names <- dir("~/countsplit_paper/Fig345_mainsims/clusterres", pattern="jan11sim*") 
cluster_res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))
names(cluster_res) <- c("j", "trueCoeff", "intercept", "pval", "fitCoeff", "cor", "eps", "type", "propImp",
                "n", "p", "prop1")

res <- res %>% mutate(intercept2=round(exp(intercept)))
cluster_res <- cluster_res %>% mutate(intercept2=round(exp(intercept)))


null_res <- res %>% filter(trueCoeff==0, prop1==0.5)
power_res <- res %>% filter(trueCoeff !=0)
null_cluster_res <- cluster_res %>% filter(trueCoeff==0, prop1==0.5)
power_cluster_res <- cluster_res %>% filter(trueCoeff !=0)

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

###### TYPE 1 ERROR FIGURE
p1 <- ggplot(data=null_res_subset, aes(sample=as.numeric(pval), col=as.factor(eps), group=eps))+geom_qq(distribution="qunif")+
  facet_grid(col=vars(intercept3))+
  geom_abline(slope=1, col="red")+ggtitle("Trajectory Estimation")+
  labs(col=expression(epsilon))+xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  coord_fixed()+
  theme_bw()
p2 <- ggplot(data=null_clusterres_subset, aes(sample=as.numeric(pval), col=as.factor(eps), group=eps))+geom_qq(distribution="qunif")+
  facet_grid(col=vars(intercept3))+
  geom_abline(slope=1, col="red")+ggtitle("Clustering")+
  xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  labs(col=expression(epsilon))+coord_fixed()+
  theme_bw()
### Takeaways: Type 1 error is controlled! Epsilon does not matter at ALL for this. 
p1+p2+plot_layout(nrow=1, guides="collect")
#ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/paper_v15/Figures/Fig2.png")

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



### POWER FIGURE
p1 <- ggplot(data=detection_res %>% filter(intercept_dist=="50% Low Intercept"), 
             aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  #facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Trajectory Estimation")+
  ylab(expression('Correlation between '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
detection_res_cluster$avcor[detection_res_cluster$avcor < 0] <- 0
p2 <- ggplot(data=detection_res_cluster %>% filter(intercept_dist=="50% Low Intercept"), aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  #facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Clustering")+
  ylab(expression('Adjusted Rand Index between  '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
p1+p2+plot_layout(guides="collect", nrow=1)


p1 <- ggplot(data=detection_res, aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Trajectory Estimation")+
  ylab(expression('Correlation between '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
detection_res_cluster$avcor[detection_res_cluster$avcor < 0] <- 0
p2 <- ggplot(data=detection_res_cluster, aes(x=abs(trueCoeff), y=avcor, col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm",
              method.args=list(family="binomial"))+
  facet_grid(col=vars(intercept_dist))+
  labs(col=expression(epsilon))+ggtitle("Clustering")+
  ylab(expression('Adjusted Rand Index between  '*widehat(L)(X^{train})*' and '*L))+
  xlab(expression(beta['1j']))+
  theme_bw()
p1+p2+plot_layout(guides="collect", nrow=2)



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

                                                                                        
p1 <- ggplot(data=power_res, aes(x=abs(fitCoeff), y=as.numeric(pval<0.05), col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm", method.args=list(family="binomial"))+
  facet_grid(cols=vars(intercept3))+
  xlim(0,15)+
  xlab(expression(beta(widehat(L)(X^{train}), bold(X)^{test}))) + ylab("Power")+ggtitle("Trajectory Estimation")+
  labs(col=expression(epsilon))+
  theme_bw()

p2 <- ggplot(data=power_cluster_res, aes(x=abs(fitCoeff), y=as.numeric(pval<0.05), col=as.factor(eps)))+
  geom_smooth(se=FALSE, method="glm", method.args=list(family="binomial"))+
  facet_grid(cols=vars(intercept3))+
  xlab(expression(beta(widehat(L)(X^{train}), bold(X)^{test}))) + ylab("Power")+ggtitle("Clustering")+
  labs(col=expression(epsilon))+
  theme_bw()+
  xlim(0,2)

p1+p2+plot_layout(guides="collect", nrow=1)
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/paper_v17/Figures/Fig4.png", width=8, height=4)
