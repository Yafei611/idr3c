
### EM---E step
e.step.normal.multi <- function(z, para, n.sample, n.dim=2, n.cluster=3){
  e.z.temp <- matrix(0, nrow=n.sample, ncol=n.cluster)
  e.z.sum <- rep(0, n.sample)
  for(i in 1:n.cluster){
    e.z.temp[,i] <- para$p[i]*dmvnorm(z, para$mu[, i], para$covm[,,i]) 
  }
  e.z.sum <- rowSums(e.z.temp)
  e.z <- e.z.temp/e.z.sum
  return(e.z)
}



#### EM ---M step
m.step.normal.multi <- function(z, e.z, modeltype, n.dim=2, n.cluster=3){
  p <- colMeans(e.z)
  mu <- matrix(0, nrow=n.dim, ncol=n.cluster)
  sd <- matrix(0, nrow=n.dim, ncol=n.cluster)
  rho <- c(0,0,0)
  mu[,1] <- rep(0, n.dim)
  sd[,1] <- rep(1, n.dim)
  rho[1] <- 0
  
  if (modeltype=="EEE"){
    mu[,c(2,3)] <- (rep(1, n.dim) %*% (t(rowSums(z))%*%e.z/colSums(e.z)/n.dim))[,c(2,3)]   
    for(i in 2:n.cluster){
      sd2_p1 <- sum(e.z[,i]*rowSums((z-mu[,i])^2))
      sd2_p2 <- sum(e.z[,i])*n.dim
      rho_p1 <- sum(e.z[,i]*(z[,1]-mu[1,i])*(z[,2]-mu[2,i]))*n.dim
      sd[,i] <- sqrt(rep(sd2_p1/sd2_p2,n.dim))
      rho[i] <- rho_p1/sd2_p1
    }
    covm <- crt_covm3(sd,rho);
  }
  
  return(list(p=p, mu=mu, sd=sd, rho=rho, covm=covm))
}


#### get loglikelihood
loglik.normal.multi <- function(z, para, n.cluster=3){
  den.m <- 0
  for(i in 1:n.cluster){
    den.m <- den.m + para$p[i]*dmvnorm(z, para$mu[,i], para$covm[,,i])
  }
  l.m <- sum(log(den.m))
  return(l.m) 
}


