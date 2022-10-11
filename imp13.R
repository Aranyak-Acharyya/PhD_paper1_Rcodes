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


K<-35

#applying isomap to obtain 1d embeddings
z_hat<-isomap(vegdist(X_hat, method="euclidean"),
              k=K,ndim=1,
              path="shortest")$points

print(z_hat)
df1<-data.frame(z_hat,y)
ggplot(df1, aes(x=z_hat, y=y)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE) +
  ylab(TeX("responses($y_i$) $\\rightarrow$")) +
  xlab(TeX("isomap embeddings($\\hat{z}_i$) $\\rightarrow$")) 


linmod<-lm(y~z_hat)
print(linmod)
summary(linmod)

tot_vec<-seq(1,100,1)
trn_vec<-sample(seq(1,length(tot_vec),1),size = 85, replace = FALSE,
              prob = NULL)
print(trn_vec)
tst_vec<-setdiff(tot_vec,trn_vec)
print(tst_vec)


beta_naive<-cov(y[trn_vec],z_hat[trn_vec])/var(z_hat[trn_vec])
alpha_naive<-mean(y[trn_vec])-beta_naive*mean(z_hat[trn_vec])

y_tst_pred<-alpha_naive+beta_naive*z_hat[tst_vec]
plot(tst_vec,y_tst_pred,type="p")

y_tst<-y[tst_vec]
z_hat_tst<-z_hat[tst_vec]

plot(z_hat_tst,y_tst)

df<-data.frame(z_hat_tst,y_tst,y_tst_pred)
dfm<-melt(df, id.vars = 'z_hat_tst')
print(dfm)
ggplot(dfm, aes(x=z_hat_tst, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("sample MSEs $\\rightarrow$")) +
  xlab(TeX("number of nodes(n) $\\rightarrow$")) +
  #ggtitle("Consistency of regression parameter estimates on
  #      known non-linear manifold") +
  scale_colour_manual(values = c("red","orange","blue"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\beta}_{true}$",
                        "sample MSE of $\\hat{\\beta}_{naive}$",
                        "sample MSE of $\\hat{\\beta}_{adj,\\hat{\\sigma}}$"
                      ))))






#s_vec<-seq(10,91,1)

tot_vec<-seq(1,100,1)
trn_vec<-sample(seq(1,length(tot_vec),1),size = 91, replace = FALSE,
                prob = NULL)
print(trn_vec)
tst_vec<-setdiff(tot_vec,trn_vec)
print(tst_vec)


g<-function(S)
{
  pred_diff_vec<-vector()
  
  for(s in 1:length(trn_vec))
  {
    
    
    s_vec<-trn_vec[1:s]
    #extracting first s of them
    zs_hat<-z_hat[s_vec]
    ys<-y[s_vec]
    
    #estimating regression parameters treating isomap embeddings 
    #as regressors
    beta_naive<-cov(ys,zs_hat)/var(zs_hat)
    alpha_naive<-mean(ys)-beta_naive*mean(zs_hat)
    
    #predicting response
    y0_pred<-alpha_naive+beta_naive*z_hat[S]
    y0<-y[S]
    pred_diff<-abs(y0_pred-y0)/mean(y)
    
    #dec<-c(s,pred_diff)
    #print(dec)
    
    
    pred_diff_vec<-c(pred_diff_vec,pred_diff)
    
  }
  
  return(pred_diff_vec)
  
}

par(mar=c(2,2,2,2))


#S_vec<-seq(s_vec[length(trn_vec)]+1,100,1)
par(mfrow=c(3,3))
for(S in tst_vec)
{
  plot(trn_vec,g(S),type="p")
  #title(main = S)
  Sys.sleep(5)
}

print(tst_vec)

g(5)

print(g(5))

print(tst_vec)
plot(trn_vec,g(65),type="p")









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
