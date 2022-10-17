set.seed(1)
n <- 500
p <- 5
X <- data.frame(matrix(rpois(n*p, 5), nrow=n, ncol=p))
X5temp <- X[,5]
X1temp <- X[,1]
X[,1] <- X5temp
X[,5] <- X1temp
clusters <- kmeans(log(X+1), centers=2)$cluster
border <- sum(clusters==1)
names(X) <- paste("gene", 1:p, sep="")


X$cell <- 1:n
Xlong <- reshape(X, direction="long", varying=list(1:p))
names(Xlong) <- c("cell", "gene", "count", "id")

#p1 <- ggplot(data=Xlong, aes(x=cell, y=gene, fill=count))+
#  geom_tile(col="black")+xlab("cell")+theme_bw()+
#  xlab("Cells ordered by index")+
#  theme(axis.ticks.x = element_blank(),
#        axis.text.x = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
#        panel.background = element_blank())+
#  scale_y_reverse(breaks=0:5, labels=c("", paste("Gene", 1:5)))+
#  ylab("")+
#  geom_text(aes(label="Cluster 2", x=border+30, y=0.25), col="white")+
#  geom_text(aes(label="Cluster 1", x=30, y=0.25), col="white")+
#  xlim(0,n+100)
#clusters

X2 <- rbind(X[clusters==1,], X[clusters==2,])
X2$cell <- 1:n
X2$cluster <- c(rep(1, sum(clusters==1)), rep(2, sum(clusters==2)))
Xlong2 <- reshape(X2, direction="long", varying=list(1:p))
names(Xlong2) <- c("cell", "cluster", "gene", "count", "id")


#p2 <- ggplot(data=Xlong2, aes(x=cell, y=gene, fill=count, col=count))+
#  geom_tile(col="black")+xlab("cell (ordered by cluster)")+
 # geom_vline(xintercept=border, col="black", lwd=2)+
 # theme_bw()+
 # theme(axis.ticks.x = element_blank(),
 #       axis.text.x = element_blank(),
  #      panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
   #     panel.border = element_blank(),
  #      panel.background = element_blank())+
  #xlab("Cells ordered by cluster")+
  #geom_text(aes(label="Cluster 2", x=border+30, y=0.25), size=4)+
  #geom_text(aes(label="Cluster 1", x=30, y=0.25), size=4)+
  #scale_y_reverse(breaks=0:5, labels=c("", paste("Gene", 1:5)))+
  #ylab("")+
  #xlim(0,n+100)+
  #geom_text(aes(label=paste("p =", round(summary(glm(X$gene1~clusters, family="poisson"))$coefficients[2,4],3)), x=550, y=1), size=2)+
  #geom_text(aes(label=paste("p =", round(summary(glm(X$gene2~clusters, family="poisson"))$coefficients[2,4],3)), x=550, y=2), size=2)+
  #geom_text(aes(label=paste("p =", round(summary(glm(X$gene3~clusters, family="poisson"))$coefficients[2,4],3)), x=550, y=3), size=2)+
  #geom_text(aes(label=paste("p =", round(summary(glm(X$gene4~clusters, family="poisson"))$coefficients[2,4],3)), x=600, y=4),size=2)+
  #geom_text(aes(label=paste("p =", round(summary(glm(X$gene5~clusters, family="poisson"))$coefficients[2,4],3)), x=600, y=5), size=2)



#p1+p2 + plot_layout(guides="collect") & scale_fill_gradient(low="white", high="darkred")
#ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat Revision August 2022/Figures/intro_heatmap.png")


Xlong$gene2 <- paste("Gene", Xlong$gene)
Xlong2$gene2 <- paste("Gene", Xlong2$gene)
Xlong2$pval <- NA
Xlong2$pval[Xlong2$gene==1] <- round(summary(glm(X$gene1~clusters, family="poisson"))$coefficients[2,4],3)
Xlong2$pval[Xlong2$gene==2] <- round(summary(glm(X$gene2~clusters, family="poisson"))$coefficients[2,4],3)
Xlong2$pval[Xlong2$gene==3] <- round(summary(glm(X$gene3~clusters, family="poisson"))$coefficients[2,4],3)
Xlong2$pval[Xlong2$gene==4] <- round(summary(glm(X$gene4~clusters, family="poisson"))$coefficients[2,4],3)
Xlong2$pval[Xlong2$gene==5] <- round(summary(glm(X$gene5~clusters, family="poisson"))$coefficients[2,4],3)

Xlong2$pval[Xlong2$pval!=0] <- paste("=", Xlong2$pval[Xlong2$pval!=0])
Xlong2$pval[Xlong2$pval==0] <- "< 0.001"

library(tidyverse)
library(patchwork)
p3 <- ggplot(data=Xlong, aes(x=count))+geom_histogram(binwidth=0.5)+
  facet_grid(rows=vars(gene2))
p4 <- ggplot(data=Xlong2, aes(x=count, group=as.factor(cluster), fill=as.factor(cluster)))+geom_histogram(position="dodge", binwidth=0.5)+
  facet_grid(rows=vars(gene2))+labs(col="Cluster")+
  geom_text(aes(x=11, y=85, label=paste("p", pval)))+ labs(col="cluster")+labs(fill="Cluster")
p3+p4 & theme_bw()&  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) & ylab("Frequency") &
  ylim(0,100)
ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat_Resubmit_October/intro_hists.eps")
#ggsave("~/Dropbox/Pseudotime : PCA NEW/Paper/Biostat reviewer responses Aug 2022/intro_hists.png")
