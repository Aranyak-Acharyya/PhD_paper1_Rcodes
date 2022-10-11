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

n<-nrow(df)

X_hat<-cbind(df$X.1,df$X.2,df$X.3,df$X.4,df$X.5,df$X.6)
y<-df$dist

plot(data.frame(X_hat))


d<-6


s_vec<-seq(10,91,1)





g<-function(S)
{
  pred_diff_vec<-vector()
  
  for(s in s_vec)
  {
    
    #determining number of nearest neighbour to run isomap
    K<-35
    
    #applying isomap to obtain 1d embeddings
    z_hat<-isomap(vegdist(X_hat, method="euclidean"),
                  k=K,ndim=1,
                  path="shortest")$points
    
    #extracting first s of them
    zs_hat<-z_hat[1:s]
    ys<-y[1:s]
    
    #estimating regression parameters treating isomap embeddings 
    #as regressors
    beta_naive<-cov(ys,zs_hat)/var(zs_hat)
    alpha_naive<-mean(ys)-beta_naive*mean(zs_hat)
    
    #predicting response
    y0_pred<-alpha_naive+beta_naive*z_hat[S]
    y0<-y[S]
    pred_diff<-abs(y0_pred-y0)/mean(y[(s+1):n])
    
    dec<-c(s,pred_diff)
    #print(dec)
    
    
    pred_diff_vec<-c(pred_diff_vec,pred_diff)
    
  }
  
  return(pred_diff_vec)
  
}

par(mar=c(2,2,2,2))


S_vec<-seq(s_vec[length(s_vec)]+1,100,1)
par(mfrow=c(3,3))
for(S in S_vec)
{
  plot(s_vec,g(S),type="l")
  #title(main = S)
  Sys.sleep(5)
}

print(g(61))



for(s in s_vec)
{
    
    #determining number of nearest neighbour to run isomap
    K<-35
    
    #applying isomap to obtain 1d embeddings
    z_hat<-isomap(vegdist(X_hat, method="euclidean"),
                  k=K,ndim=1,
                  path="shortest")$points
    
    #extracting first s of them
    zs_hat<-z_hat[1:s]
    ys<-y[1:s]
    
    #estimating regression parameters treating isomap embeddings 
      #as regressors
    beta_naive<-cov(ys,zs_hat)/var(zs_hat)
    alpha_naive<-mean(ys)-beta_naive*mean(zs_hat)
    
    #predicting response
    y0_pred<-alpha_naive+beta_naive*z_hat[S]
    y0<-y[S]
    pred_diff<-abs(y0_pred-y0)
    
    dec<-c(s,pred_diff)
    print(dec)
  
    
    pred_diff_vec<-c(pred_diff_vec,pred_diff)
  
}



print(pred_diff_vec)
plot(s_vec,pred_diff_vec,type="l",col="blue")






+
  #ggtitle("Consistency of regression parameter estimates on
  #      known non-linear manifold") +
  scale_colour_manual(values = c("red","orange","blue"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\beta}_{true}$",
                        "sample MSE of $\\hat{\\beta}_{naive}$",
                        "sample MSE of $\\hat{\\beta}_{adj,\\hat{\\sigma}}$"
                      ))))










