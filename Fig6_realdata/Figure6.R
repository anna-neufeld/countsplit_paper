setwd("~/countsplit_paper/Fig6_realdata")
load("monocle_workspace_June13.RData")

library(tidyverse)
### Really affected by ridiculuously small numbers.
p1 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_full_naive_hand, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_test_naive_hand, col="Test Data"))+
  geom_abline(a=0,b=1)+
  scale_x_log10()+scale_y_log10()+
  theme_bw()+coord_fixed()+xlab("Count Splitting p-values")+ylab("Naive p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Naive gets to use")

#### This ignores all p-values that are less than 2e-16.
p2 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_full_naive_hand, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand, y=pvals_test_naive_hand, col="Test Data"))+
  geom_abline(a=0,b=1)+
  xlim(2e-16,1)+ylim(2e-16,1)+
  scale_x_log10(limits=c(2e-16,1))+scale_y_log10(limits=c(2e-16,1))+
  theme_bw()+coord_fixed()+xlab("Count Splitting p-values")+ylab("Naive p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Naive gets to use")
p2

#### This jitters all p-values that are less than numerical precision!!!
pvals_countsplit_full_hand2 <- pvals_countsplit_full_hand
pvals_full_naive_hand2 <- pvals_full_naive_hand
pvals_test_naive_hand2 <- pvals_test_naive_hand

pvals_countsplit_full_hand2[pvals_countsplit_full_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_countsplit_full_hand < 2e-16)))
pvals_full_naive_hand2[pvals_full_naive_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_full_naive_hand < 2e-16)))
pvals_test_naive_hand2[pvals_test_naive_hand < 2e-16] <- jitter(rep(2e-16, sum(pvals_test_naive_hand< 2e-16)))
p3 <- ggplot(data=NULL)+
  geom_point(aes(x=pvals_countsplit_full_hand2, y=pvals_full_naive_hand2, col="Full Data"))+
  geom_point(aes(x=pvals_countsplit_full_hand2, y=pvals_test_naive_hand2, col="Test Data"))+
  geom_abline(a=0,b=1)+
  scale_x_log10()+scale_y_log10()+
  theme_bw()+coord_fixed()+xlab("Count Splitting p-values")+ylab("Naive p-values")+
  ggtitle("All Days")+
  coord_fixed()+
  labs(col="Naive gets to use")
p3

#### NULL STORY!!!

### Make sure we don't include any genes that are just literally all 0s.
good_indices <- colSums(X0[, hvg.cm.var]) > 0
cols <- c("#00B81F", "#00A5FF", "#F8766D")
p4 <- ggplot(data=NULL)+
  geom_qq(aes(sample=pvals_countsplit_0[good_indices], col="Count Split"), distribution="qunif")+
  geom_qq(aes(sample=pvals_full_naive0[good_indices], col="Naive"), distribution="qunif")+
  geom_qq(aes(sample=pvals_test_naive0[good_indices], col="Test Naive"), distribution="qunif")+
  geom_abline()+theme_bw()+coord_fixed()+xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  labs(col="Method")+
  scale_color_manual(values=cols[c(3,1,2)])+
  ggtitle("Day 0 Only")
p4


p3+p4


## Type 1 errors!
mean(pvals_countsplit_0 < 0.05, na.rm=TRUE)
mean(pvals_full_naive0 < 0.05, na.rm=TRUE)
mean(pvals_test_naive0 < 0.05, na.rm=TRUE)

