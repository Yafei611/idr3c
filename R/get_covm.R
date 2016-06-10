
#' Create the covariance matrix for each component of the normal mixture.
#'
#' This function helps to create the 3 covariance matrix (2-by-2) for each component of the normal mixture. 
#'
#' @param sds The standard deviations for each of the 3 normal components. sd should be a 2-by-3 matrix, each row are standard deviations for the same components.
#' @param rhos The Pearson correlation between expression data sources for each component. Rho should be a vector with a length of 3. According to the idr3c model, the first element should be 0 under most conditions.
#' @export
#' @examples
#' sds = cbind(c(1,1),c(1,1),c(1,1))
#' rhos = c(0,0.84,0.84)
#' covms = crt_covm3(sds=sds,rhos=rhos)
#' print(covms)

### creating covariance matrix for 3 normal components
crt_covm3 <- function(sds,rhos){
  covm <- array(NA, c(2,2,3))
  covm[,,1] <- matrix(c(sds[1,1]^2,sds[1,1]*sds[2,1]*rhos[1],sds[1,1]*sds[2,1]*rhos[1],sds[2,1]^2), nrow=2);
  covm[,,2] <- matrix(c(sds[1,2]^2,sds[1,2]*sds[2,2]*rhos[2],sds[1,2]*sds[2,2]*rhos[2],sds[2,2]^2), nrow=2);
  covm[,,3] <- matrix(c(sds[1,3]^2,sds[1,3]*sds[2,3]*rhos[3],sds[1,3]*sds[2,3]*rhos[3],sds[2,3]^2), nrow=2);
  return(covm)
}





