#' @title Test for anomalies in wavelet space via conditional simulation
#' 
#' @description Test for anomalies using \code{EFDR} and conditional simulation. The noisy image 
#' can be partially observed, or/and aggregated at different resolutions
#' @param Zvec vector of observations such that \eqn{Ztilde = H.Z}
#' @param H matrix mapping the fine-resolution image \code{Z} in vector form to \code{Ztllde}. Must have 
#'          as many rows as \code{Ztilde} and \code{n1} x \code{n2} columns
#' @param n1 number of rows in fine-resolution image
#' @param n2 number of columns in fine-resolution image
#' @param rho_est_method method with which to estimate the level of exchangeability rho; can be either "CPL" (copula model) or "MOM" (method of moments)
#' @param iter.cs number of conditional simulations to carry out
#' @param wf type of wavelet to employ. Defaults to `la8', the Daubechies orthonormal compactly supported wavelet of length \code{L = 8} (Daubechies, 1992), least asymmetric family. Other options include `haar' (Haar wavelet), `fk8' (Fejer-Korovkin wavelet with \code{L=8}) and `mb8' (minimum-bandwidth wavelet with \code{L=8}). Please type `\code{waveslim::wave.filter}' in the console for a full list of wavelet names
#' @param J number of resolutions to employ in wavelet decomposition
#' @param alpha significance level at which tests are carried out
#' @param n.hyp number of hypotheses tests to carry out with EFDR. If a vector is supplied, the optimal one from 
#' the set of proposed number of tests is chosen
#' @param b the number of neighbours to consider in EFDR
#' @param iteration number of Monte Carlo iterations to employ when determining which of the proposed number of tests 
#' in \code{n.hyp} is the optimal number of tests
#' @param parallel number of cores to use with parallel backend; needs to be an integer less than or equal to the number of available cores
#' @return List with three fields:
#' \describe{
#'  \item{\code{filtered}}{the discrete wavelet transform containing the anomalous wavelet coefficients in the signal}
#'  \item{\code{Z}}{the image containing the anomalous wavelets in the signal}
#'  \item{\code{reject_coeff}}{indices of wavelets under which the null hypothesis of no anomaly was rejected}
#'  \item{\code{pvalue_ordered}}{ordered p-values under the null hypothesis. The column names indicate the 
#'                wavelet to which the p-value belongs}
#'  \item{\code{nhat}}{the number of tests carried out.}
#' }
#' @export
#' @references 
#' Daubechies, I. (1992) Ten Lectures on Wavelets, CBMS-NSF Regional Conference Series in Applied Mathematics, SIAM: Philadelphia.
#' 
#' Shen, X., Huang, H.-C., and Cressie, N. 'Nonparametric hypothesis testing for a spatial signal.' Journal of the American Statistical Association 97.460 (2002): 1122-1140.
#'
#' @examples
#' ## Set up experiment
#' n <- 32       # 32 x 32 images
#' r <- 10       # signal of size 10 x 10
#' h <- 5        # intensity of 5
#' grid <- 8     # aggregated to 8 x 8 image
#' parallel <- 4 # use 4 cores
#'
#' ## Simulate the pixel-level data
#' raw_grid <- expand.grid(x = seq(1, n), y = seq(1, n))
#' df <- data.frame(raw_grid)                        # spatial grid
#' dd <- as.matrix(dist(raw_grid, diag = TRUE))      # distance matrix
#' Sigma <- exp(-dd/5)                               # cov. fn.
#' diag(Sigma) <- 1                                  # fix diagonal
#' L <- t(chol(Sigma))                               # lower Cholesky factor
#' mu <- matrix(0, n, n)                             # zero mean
#' mu[(n/2-r/2):(n/2+r/2), (n/2-r/2):(n/2+r/2)] <- h # add signal
#' Z <- mu + matrix(L %*% rnorm(n^2), n, n)          # simulate data
#' 
#' ## Construct H (aggregation) matrix
#' H <- matrix(0, grid^2, n^2)
#' for(i in 1:grid^2) {
#'   ind <- rep(rep(c(0L,1L,0L),
#'              c((n/grid)*((i-1)%%grid),n/grid,(n-n/grid-n/grid*((i-1)%%grid)))),
#'           n/grid)
#'   H[i,which(c(rep(0L,(ceiling(i/grid)-1)*n^2/grid),ind) == TRUE)] <- 1/(n/grid)^2
#' }
#'
#' ## Aggregate the signal
#' z_tilde <- c(H %*% c(Z))
#'
#' ## Run EFDR using conditional simulation
#' \dontrun{out2 <- test.efdr.condsim(Zvec = z_tilde, H = H, n1 = n, n2 = n, 
#'                                    parallel = parallel)}
test.efdr.condsim <- function(Zvec, H, n1, n2, rho_est_method = c("CPL", "MOM"), iter.cs = 100,
                              wf = "la8", J = 2, alpha = 0.05, n.hyp = 100, b = 11,
                              iteration = 200, parallel = 1L) {
  
  if(is.matrix(Zvec))
    if(ncol(Zvec > 1)) stop("Zvec needs to be a vector")
  if(!length(Zvec) == nrow(H)) stop("H needs to have as many rows as there are elements in Zvec")
  if(!(n1*n2 == ncol(H))) stop("H needs to have n1 x n2 columns")
  if(!(n1 %in% 2^(1:10))) stop("n1 needs to be a dyadic power and less or equal to 1024")
    if(!(n2 %in% 2^(1:10))) stop("n2 needs to be a dyadic power and less or equal to 1024")
  rho_est_method <- match.arg(rho_est_method)
  
  theta_ML <- rep(NA,3*J+1)
  pvalue.all <- rep(NA,iter.cs)
  mu.hat.all <- array(NA,c(iter.cs,n1,n2))
  raw_grid <- expand.grid(x=seq(1,n1),y=seq(1,n2)) 
  df <- data.frame(raw_grid)
  dd <- as.matrix(dist(raw_grid,diag=T))
  dwt.mu <- dwt.2d(x=matrix(0,n1,n2), wf=wf, J=J)
  W <- matrix(0,n1*n2,n1*n2)
  
  for(i in 1:(n1*n2)) {
    indvar <- unlist(dwt.mu)*0
    indvar[i] <- 1
    indvar_list <- .relist.dwt(indvar,dwt.mu)
    W[i,] <- idwt.2d(indvar_list)
  }
  
  num <- c(0,cumsum(rep((n1/2^(seq(1,J)))*(n2/2^(seq(1,J))),e=3)),n1*n2)
  
  cat("Estimating spatial parameters ...\n")
  phi.hat <- optimize(negloglik_exp, interval=c(0,20), x=Zvec, dd=dd, H=H)$min
  sigmasq.hat <- c(t(Zvec) %*% solve(H %*% exp(-dd/phi.hat) %*% t(H)) %*% Zvec/nrow(H))
  theta.temp <- sigmasq.hat*colSums(t(W)*(exp(-dd/phi.hat) %*% t(W)))
  
  for(j in 1:(3*J+1)) 
    theta_ML[j] <- mean(theta.temp[seq(num[j]+1,num[j+1])])
  
  Sigma11 <- crossprod(chol(diag(rep(theta_ML,diff(num)))) %*% W)
  Sigma12 <- Sigma11 %*% t(H) 
  Sigma21 <- t(Sigma12)
  Sigma22 <- H %*% Sigma12
  Sigma22inv_Sigma21 <- solve(Sigma22,Sigma21)
  cholSigma11 <- chol(Sigma11)
  Temp1 <- crossprod(cholSigma11 %*% t(H))
  cholTemp1 <- chol(Temp1)
  Temp2 <- crossprod(t(solve(cholTemp1)) %*% H %*% Sigma11)
  ev <- eigen(Sigma11-Temp2, symmetric=T)
  LL <- ev$vectors %*% (t(ev$vectors)*sqrt(pmax(ev$values,0)))
  mean12 <- t(Sigma22inv_Sigma21) %*% Zvec
  
  cat("Conditionally simulating and doing EFDR on each simulation ...\n")
  for(i in 1:iter.cs) {
    cat(paste0("\n ----------------- Simulation ", i, " ----------------- \n"))
    df$z <- mean12 + LL%*%rnorm(n1*n2)
    Z_mat <- df.to.mat(df)
    out <- test.efdr(Z_mat, wf=wf, J=J, n.hyp=n.hyp, parallel=parallel)
    pvalue.all[i] <- min(out$pvalue_ordered*out$n.hyp/seq(1,out$n.hyp))
    mu.hat.all[i,,] <- out$Z
  }
  
  pvalue.all <- pmax(pvalue.all,1e-10)
  ti <- -2*log(pvalue.all)
  
  cat("Estimating rho ...\n")
  if(rho_est_method=="CPL") rho.hat <- rho_est_MCL(ti)
  if(rho_est_method=="MOM") rho.hat <- rho_est_MoM(ti)
  ahat <- iter.cs/(1+(iter.cs-1)*rho.hat)
  bhat <- 1/(2*(1+(iter.cs-1)*rho.hat))
  pvalue <- 1-pgamma(q=sum(ti), shape=ahat, rate=bhat)
  mu.hat <- apply(mu.hat.all,c(2,3),mean)
  list(pvalue = pvalue, mu.hat = mu.hat, phi.hat = phi.hat, 
       sigmasq.hat = sigmasq.hat, rho.hat = rho.hat, pvalue.ave = mean(pvalue.all))
}

negloglik_exp <- function(phi,x,dd,H) {
  sigma <- H %*% exp(-dd/phi) %*% t(H)
  determinant(sigma)$mod+nrow(H)*log(crossprod(solve(sigma,x),x))
}

negloglik_pair <- function(r,U) {
  M <- length(U)
  pair <- t(combn(1:M,2))
  U.pair <- cbind(U[pair[,1]],U[pair[,2]])
  norm.cop <- normalCopula(r, dim=2)
  -sum(dCopula(U.pair, norm.cop, log=T))
}

rho_est_MCL <- function(X) {
  M <- length(X)
  U <- pexp(X, rate=0.5)
  r_est <- optimize(negloglik_pair, interval=c(0,1), U=U)$minimum
  norm.cop <- normalCopula(r_est, dim=M)
  Usim <- rCopula(10000, norm.cop)
  Xexp <- qexp(Usim, rate=0.5)
  rho <- mean(cor(Xexp)[lower.tri(cor(Xexp))])
  rho
}

rho_est_MoM <- function(X) {
  M <- length(X)
  sqD <- outer(X, X, FUN='-')^2
  sumsqD <- sum(sqD)/2
  rho <- 1-(sumsqD/(M-1))/sum((X-2)^2)
  rho <- max(rho, -(1/(M-1))+1e-8)
  rho
}
