
#' Differential expression analysis by integrating data from multiple sources
#'
#' Mapping genes expression data from different sources to a normal mixture. Then performing model based cluster and computing irreproducibility rate (idr) of being differential expressed for each genes.
#'
#' @import mvtnorm 
#' @import stats
#'
#' @param x A matrix or data frame of log-fold change of normalized gene expression data, i.e., of non-negative values. The rows correspond to observations (e.g., expression level for a gene), the columns correspond to data sources (labs/platforms etc.). Note that technical replicates (i.e., several sequencing ruins/lanes from the same sample) should be took the average.
#' @param para A list that contains the parameter of the normal mixture where the pseudo data were generated.
#' \itemize{
#'  \item{$p} {The size of each components (proportion of non/up/down regulated genes)}
#'  \item{$mu} {Means for each normal component. para$mu should be a 2-by-3 matrix, each row are means for the same components.}
#'  \item{$sd} {Standard deviations for each normal component. para$sd should be a 2-by-3 matrix, each row are standard deviations for the same components.}
#'  \item{$rho} {Pearson correlation for each normal component.}
#' }

#' @param f_ctrl A list that contains the parameter controls pseudo-EM procedure. The default \code{f_ctrl} usually works well. 
#' \itemize{
#' \item{$max.t.inner} {The maximum iterations numbers of the inner EM algorithm}
#' \item{$max.t.outer} {The maximum number of times that re-generating pseudo data}
#' \item{$b.em} {Stop criterion for inner EM algorithm.}
#' \item{$b.outer} {Stop criterion for outer iteration.}
#' }
#' @param disp A list control the display of the trace of parameter estimation.
#' \itemize{
#'  \item{$z.rankplot} {Logical. Should the plot of pseudo data be displayed when updated?}
#'  \item{$z.rankplot.cut} {The proportion of samples with low idr that will be marked on the plot.}
#'  \item{$outer.trace} {Logical. Should the trace of the parameter estimation be displayed?}
#'  \item{$inner.trace} {Logical. Should the trace of the parameter estimation for the inner EM be displayed?}
#' }
#' @return \code{IDR.3component} produces a list that contain idr value, parameter estimation and other necessary components.
#' \itemize{
#'  \item{$para} {The final estimation of the normal mixture}
#'  \item{$loglik.trace.in} {Trace of log-likelihood of EM produce}
#'  \item{$loglik.trace.ot} {Trace of log-likelihood of outer iterations}
#'  \item{$idr} {the posterior probability that a gene is in the non-DEGs group}
#'  \item{$IDR} {computed irreproducibility rate}
#' }
#' @details Please refer \bold{Lyu, Y., & Li, Q. (2016). A semi-parametric statistical model for integrating gene expression profiles across different platforms. BMC Bioinformatics, 17(1), 51.}
#' @export
#' @examples
#' x = expr
#' plot(x[,1],x[,2],xlim=c(-8,8),ylim=c(-8,8),cex = .4);
#' idr.out = IDR.3component(x = x)



######main
IDR.3component <- function(x, 
                           para = list(p = c(0.4,0.3,0.3), mu = cbind(c(0,0),c(3,3),c(-3,-3)), sd = cbind(c(1,1),c(1,1),c(1,1)), rho = c(0,0.84,0.84)),
                           f_ctrl = list(max.t.inner = 10000,max.t.outer = 10000,b.em = 1e-6,b.outer = 1e-3),
                           disp = list(z.rankplot = FALSE,z.rankplot.cut = 0.2,outer.trace = TRUE,inner.trace = FALSE)){
  
  # initialization

  modeltype="EEE"
  para$covm = crt_covm3(para$sd,para$rho)

  n.dim = ncol(x)
  n.sample = nrow(x)
  n.cluster = 3
  
  e.z=matrix(0,ncol=3,nrow=n.sample);
  
  outer.t = 0;
  inner.t = 0;
  loglik.outer.trace <- NULL;
  loglik.inner.trace <- NULL;
  loglik.ot = -9*10^9;
  loglik.ot0 = 0;
  
  # getting ecdf of x
  x.cdf = matrix(0, nrow=n.sample, ncol=n.dim)
  afactor = n.sample/(n.sample+1)
  
  for(i in 1:n.dim){
    x.cdf.func = ecdf(x[,i])
    x.cdf[,i] = x.cdf.func(x[,i])*afactor
  }
  
  while(1){
    outer.t = outer.t+1;
    loglik.ot0 = loglik.ot;
    
    z = get.pseudo.mix.multi(x.cdf, para)
    
    if (disp$z.rankplot) {
      plot(z[,1], z[,2], xlab="d1", ylab="d2", cex = .4, main = ("latent varibale z"),
           col = as.numeric(rank(e.z[,1])<n.sample*disp$z.rankplot.cut)+1);
    }
    
    loglik.ot = loglik.normal.multi(z, para)
    loglik.outer.trace[outer.t] = loglik.ot;
    
    if (disp$outer.trace) {
      cat("# t=",outer.t," loglik=",loglik.ot,"\n",sep=" ")
    }

    if ((abs(loglik.ot-loglik.ot0)<abs(f_ctrl$b.outer*loglik.ot0))||(outer.t+1>f_ctrl$max.t.outer)) {break}
    
    ## inner
    loglik.in = -9*10^9;
    loglik.in0 = 0;
    
    while(1){
      inner.t = inner.t+1;
      loglik.in0 = loglik.in;
      
      # EM for latent structure
      e.z = e.step.normal.multi(z, para, n.sample, n.dim=2, n.cluster=3);
      para = m.step.normal.multi(z, e.z, modeltype);  
      
      loglik.in = loglik.normal.multi(z,para);
      loglik.inner.trace[inner.t] = loglik.in;
      
      if (disp$inner.trace) {
        cat("t=",outer.t,inner.t," loglik=",loglik.in,
            "  p=",round(para$p,2),
            " mu=",round(para$mu[,c(2,3)],2),
            " sd=",round(para$sd[,c(2,3)],2),
            " rho=",round(para$rho[c(2,3)],2),
            "\n",sep=" ")
      }
      
      if ((abs(loglik.in-loglik.in0)<abs(f_ctrl$b.em*loglik.in0))||(inner.t+1>f_ctrl$max.t.inner)) {break}
    }
  }
  
  # probability of being in the random component
  idr <- e.z[,1]
  o <- order(idr)
  idr.o <- idr[o]
  idr.rank <- rank(idr.o, ties.method="max")
  
  top.mean <- function(index, x){mean(x[1:index])}
  
  # IDR of selected peaks
  IDR.o <- sapply(idr.rank, top.mean, idr.o)
  IDR <- rep(NA, length(IDR.o))
  IDR[o] <- IDR.o
  
  return(list(para=para,loglik.trace.in=loglik.inner.trace, loglik.trace.ot=loglik.outer.trace, idr=idr, IDR=IDR))
}




