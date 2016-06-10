
### get pseudo data
get.pseudo.mix.multi <- function(x, para) {
  
  n.dim=2
  n.cluster=3
  nw=1000
  m_sd=3
  n.sample=nrow(x)
  
  quan.x <- matrix(0, nrow=n.sample, ncol=n.dim)
  
  for(i in 1:n.dim){
    w.cdf <- rep(0, nw)
    min.w <- min(para$mu[i,] - m_sd*para$sd[i,])
    max.w <- max(para$mu[i,] + m_sd*para$sd[i,])
    w <- seq(min.w, max.w, length=nw)  
    
    for(k in 1:n.cluster){
      w.cdf <- w.cdf + para$p[k]*pnorm(w, mean=para$mu[i,k], sd=para$sd[i,k])
    }
    
    for(j in 1:(nw-1)){
      index <- which(x[,i] >= w.cdf[j] & x[,i] < w.cdf[j+1])
      if(length(index)>0)
        quan.x[index,i] <- (x[index,i]-w.cdf[j])*(w[j+1]-w[j])/(w.cdf[j+1]-w.cdf[j]) +w[j]
    }
    
    index <- which(x[,i] < w.cdf[1])
    if(length(index)>0)
      quan.x[index,i] <- w[1]
    
    index <- which(x[,i] > w.cdf[nw])
    if(length(index)>0)
      quan.x[index,i] <- w[nw]  
  }
  return(quan.x)
}

