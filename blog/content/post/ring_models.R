library(mvtnorm)
library(colourvalues)
library(matrixStats)

logistic <- function(x) 1 / (1+exp(-x))
logit <- function(p) log(p / (1-p))


# Set up parameters
N <- 1e4
# Treatment dimensionality
m <- 10

# Latent class probability
pU <- 0.5

# Variance of treatments
sig <- 1
sigma2X_U <- function(U){ ifelse(U==0, sig^2, (2 * sig)^2) }
mean_mat <- matrix(logistic(c(2, 1, -1, -2)), nc=2)

# Mean of outcomes
mu_Y_factory <- function(mean_mat){
  # number of levels of each variable
  nU <- dim(mean_mat)[2]
  nX <- dim(mean_mat)[1]
  function(U, X){
    # Set up X as row index, U as column index
    Xd <- as.numeric(sqrt(rowMeans(X^2)) > 1.5)
    idx <- U * nU + Xd + 1
    mean_mat[idx]
  }
}
mu_Y <- mu_Y_factory(mean_mat)

# True counterfactual mean
EY_do_x <- function(x, mu_Y, pU){
  pU * mu_Y(1, x) + (1-pU) * mu_Y(0, x)
}

# Mean of proxy
mu_W <- function(U) ifelse(U==0, -1, 1)

gendat_rings <- function(N, m, pU, sigma2X_U, mu_Y, mu_W){ 
  # Simulate Data
  nU <- rbinom(1, N, pU)
  U <- c(rep(0, N-nU), rep(1, nU))
  X0 <- rmvnorm(N - nU, rep(0, m), diag(sigma2X_U(0), m))
  X1 <- rmvnorm(nU, rep(0, m), diag(sigma2X_U(1), m))
  
  X <- rbind(X0, X1)
  W <- rnorm(N, mu_W(U), 1)
  Y <- rbinom(N, 1, mu_Y(U, X))
  
  return(list(U=U, X=X, Y=Y, W=W, s2X=sigma2X_U(U)))
}



ring_ll <- function(pars, X, Y, W=NULL, debug=FALSE){
  m <- dim(X)[2]
  
  pU <- logistic(pars[1])
  sig2X_U <- function(U) ifelse(U==0, exp(pars[2]), exp(pars[3]))
  mu_Y <- mu_Y_factory(matrix(logistic(pars[4:7]), nc=2))
  
  if(!is.null(W))
    mu_W <- function(U) ifelse(U==0, pars[8], pars[9])
 
  # XY contribution 
  U0_XY_ll <-  log(1-pU) + dmvnorm(X, rep(0, m), diag(sig2X_U(0), m), log=TRUE) +
    dbinom(Y, 1, mu_Y(0, X), log=TRUE)
  U1_XY_ll <- log(pU) + dmvnorm(X, rep(0, m), diag(sig2X_U(1), m), log=TRUE) +
    dbinom(Y, 1, mu_Y(1, X), log=TRUE)
  
  # Proxy contribution
  U0_W_ll <- if(is.null(W)) 0 else dnorm(W, mu_W(0), 1, log=TRUE)
  U1_W_ll <- if(is.null(W)) 0 else dnorm(W, mu_W(1), 1, log=TRUE)
  
  if(debug) browser()
  
  sum(apply(cbind(U0_XY_ll + U0_W_ll, U1_XY_ll + U1_W_ll), 1, matrixStats::logSumExp))
}

display_res <- function(res, truth){
  pars <- res$par
  
  display_pars <- c(logistic(pars[1]), exp(pars[2:3]), logistic(pars[4:7]))
  display_truth <- c(logistic(truth[1]), exp(truth[2:3]), logistic(truth[4:7]))
  par_names <- c("pU", "s2X_0", "s2X_1", paste("muY", c('00', '01', '10', '11'), sep=''))
  
  if(length(res$par) > 7){
    display_pars <- c(display_pars, pars[8:9])
    display_truth <- c(display_truth, truth[8:9])
    par_names <- c(names, "muW0", "muW1")
  }
  
  names(display_pars) <- par_names
  
  print(rbind(display_pars, display_truth))
}

gd <- gendat_rings(N, m, pU, sigma2X_U, mu_Y, mu_W)
true_pars <- c(logit(pU), log(sig^2), log((2 * sig)^2), logit(as.vector(mean_mat)))
init_pars <- true_pars
init_pars <- runif(7, -1, 1)
#res <- optim(init_pars, function(pars) ring_ll(pars, gd$X, gd$Y), method='L-BFGS-B',
#      control=list(fnscale=-1, factr=1e3, maxit=1e4, trace=1))
display_res(res, init_pars)

## Plotting Below

# 2-D polar coordinates to cartesion coordinates  
polar2cartesian <- function(r, ang){
  cbind(r*cos(ang), r*sin(ang))
}

# Plot 2-D polar projection of high-dimensional spherical Gaussians from model
# with decision boundary for Zhat
plot_polar_proj <- function(sig, m=100, N=1e3, decision=1.5, ...){
  x0 <- rmvnorm(N, rep(0,m), diag(sig^2, m))
  x1 <- rmvnorm(N, rep(0,m), diag((2*sig)^2, m))
  U <- c(rep(0, N), rep(1, N))
  
  plot_polar_by_U(sig, m, rbind(x0, x1), U, ...)
}
 
plot_polar_by_U <- function(sig, m, X, U, Y=NULL, ...){
  # Polar projection of high-dimensional vector.
  proj_coords <- function(x1){
    d1 <- c(1, rep(0, m-1))
    d2 <- c(0, 1, rep(0, m-2))
    proj <- cbind(d1, d2)
    x1p <- x1 %*% proj
    
    norm1 <- sqrt(rowSums(x1^2)/m)
    raw_dir <- acos(x1p[,1]/sqrt(x1p[,1]^2+x1p[,2]^2))
    dir1 <- ifelse(x1p[,2] > 0, raw_dir, 2*base::pi-raw_dir)
    cbind(norm1*cos(dir1), norm1*sin(dir1))
  }
  
  ps <- proj_coords(X)
  boundary <- polar2cartesian(1.5*sig, seq(0, 2*base::pi, length.out=50))
  
  col0 <- 'red'
  col1 <- 'blue'
  shape0 <- 4 
  shape1 <- 16
  if(is.null(Y)){
    colors <- ifelse(U==0, col0, col1)
    pch <- 46
    cex <- 2
  } else {
    colors <- colour_values(Y, palette='viridis')
    pch <- ifelse(U==0, shape0, shape1)
    cex <- 0.75
  }
  plot(ps[,1], ps[,2], col=colors, pch=pch, cex=cex,
       xlab=NA, ylab=NA, main=sprintf("m = %d", m), ...)
  lines(boundary[,1], boundary[,2])
}
png('rings_positivity.png', width=5*150, height=floor(5/3*150), units='px')
par(mfrow=c(1, 3), mar=c(2,2,2,2), bty='n')
lims <- c(-3, 3)
plot_polar_proj(1, m=2, ylim=lims, xlim=lims)
plot_polar_proj(1, m=20, ylim=lims, xlim=lims)
plot_polar_proj(1, m=200, ylim=lims, xlim=lims)
legend('bottomright', c("U = 0", "U = 1", expression(hat(U)~boundary)),
       col=c("red", "blue", "black"), pch=c(46, 46, NA), lty=c(NA, NA, 1),
       bty='o', box.col='white', bg="#00000011", pt.cex=8)
dev.off()

par(mfrow=c(1,1))
gd <- gendat_rings(1e3, 200, pU, sigma2X_U, mu_Y, mu_W)
plot_polar_by_U(sig, 200, gd$X, gd$U, gd$Y, ylim=lims, xlim=lims)
  