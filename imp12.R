library(Matrix)
library(tidyverse)
library(MASS)
library(MVA)
library(maps)
library(MCMCpack)
library(doBy)
library(utilities)
library(irlba)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(latex2exp)
library(vegan3d)
library(car)
library(igraph)
library(Rdimtools)
library(combinat)

#loading the dataset
load("/cloud/project/df-forAranyak.RData")

X_hat<-cbind(df$X.1,df$X.2,df$X.3,df$X.4,df$X.5,df$X.6)
y<-df$dist


s<-10
d<-6
n_vec<-seq(50,100,5)

pred_diff_avg_vec<-vector()

for(n in n_vec)
{
  
  XX_hat<-X_hat[1:n,]
  KK<-n/2
  K<-as.integer(KK)
  z_hat<-isomap(vegdist(XX_hat, method="euclidean"),
                k=K,ndim=1,
                path="shortest")$points
  
  zs_hat<-z_hat[1:s]
  ys<-y[1:s]
  
  beta_naive<-cov(ys,zs_hat)/var(zs_hat)
  alpha_naive<-mean(ys)-beta_naive*mean(zs_hat)
  
  y0_pred<-alpha_naive+beta_naive*z_hat[(s+1)]
  y0<-y[(s+1)]
  pred_diff_avg<-abs(y0_pred-y0)
  
  dec<-c(n,pred_diff_avg)
  print(dec)
  
  pred_diff_avg_vec<-c(pred_diff_avg_vec,pred_diff_avg)

}

plot(n_vec,pred_diff_avg_vec,type="l",col="red")







