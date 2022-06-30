library(Matrix)
library(irlba)

n<-1000
T<-5
d<-4


gam<-function(P)
{
  pvec<-c(P[lower.tri(P,diag=TRUE)])
  avec<-rbinom(length(pvec),1,pvec)
  A<-matrix(nrow=nrow(P),ncol=ncol(P))
  A[lower.tri(A,diag=TRUE)]<-avec
  A[upper.tri(A)]<-t(A)[upper.tri(A)]
  return(A)
}

k_vec<-vector()

for(trial in 1:T)
{
  t<-runif(n,min=0,max=1)
  X<-cbind(t/2,t/2,t/2,t/2)
  P<-X%*%t(X)
  
  #generating adjacency matrix
  A<-gam(P)
  
  #finding raw (unrotated) ASE estimates of altent positions
  A_irlba<-irlba(A,d)
  X_hat_raw<-A_irlba$u%*%diag(A_irlba$d)^0.5
  
  #finding consistent adjacency spectral estimates of latent positions
  vals<-svd(t(X_hat_raw)%*%X)
  W<-vals$u%*%t(vals$v)
  X_hat<-X_hat_raw%*%W
  
  #finding the norm of the diff
  k<-norm(X_hat-X,type = "i")
  k_vec<-c(k_vec,k)
}

print(k_vec)
save(k_vec,file = "large.RData")