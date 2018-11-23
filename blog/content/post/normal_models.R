library(mvtnorm)
library(magrittr)

## MVN model

# Model without proxy
N <- 1e4

# Start with 10 treatments
m <- 10
a_0 <- 2
b_0 <- 1
s2X <- 2
# Scale coefficients so that the variance of Y stays fixed even as m grows.
alpha_ <- true_alpha <- rep(a_0, m) / sqrt(m)
beta_ <- rep(b_0, m) / sqrt(m)
gamma_ <- 2

sigma2_U_ <- true_sigma2_U <- 3
sigma2_X_ <- rep(s2X, m)
sigma2_Y_ <- 10

gendat2 <- function(N, numproxies, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y){
  m <- length(alpha) + numproxies
  alpha <- c(alpha, rep(alpha[1], numproxies))
  beta <- c(beta, rep(0, numproxies))
  
  U <- rnorm(N, 0, sqrt(sigma2_U))
  Uvec <- rep(U, each=m)
  
  X <- matrix(rnorm(m * N, alpha*Uvec, sqrt(sigma2_X)), nr=m)
  Y <- colSums(as.vector(beta) * X) + gamma * U + rnorm(N, 0, sqrt(sigma2_Y))
  list(U = U, X = X, Y = Y)
}

emp_covmat <- function(gdat){
  cov(cbind(gdat$U, t(gdat$X), gdat$Y))
}

partial_cov <- function(full_covmat, numproxies, both=FALSE){
  m <- dim(full_covmat)[1] - 2 - numproxies
  XX <- full_covmat[1 + (1:(m+numproxies)), 1 + (1:(m+numproxies))]
  UY <- full_covmat[c(1, m+2+numproxies), c(1, m+2+numproxies)]
  UY_X <- full_covmat[1 + (1:m), c(1, m+2+numproxies)]
  
  UYdotX <- UY - t(UY_X) %*% solve(XX) %*% UY_X
  if(both) list(UY=UY, UYdotX = UYdotX) else UYdotX
}

make_covmat <- function(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y){
  XX_theory <- alpha %*% t(alpha) * sigma2_U + diag(sigma2_X)
  YY_theory <- (t(beta) %*% alpha + gamma)^2 * sigma2_U + t(beta) %*% diag(sigma2_X) %*% beta+ sigma2_Y
  XY_theory <- XX_theory %*% beta + sigma2_U * gamma * alpha
  Sigma_theory <- rbind(cbind(XX_theory, XY_theory), c(as.vector(XY_theory), YY_theory))
  Sigma_theory 
}

beta_C <- function(C, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y, check=TRUE){
  covmat <- make_covmat(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)
  m <- nrow(covmat)-1
  
  SYX <- covmat[1:m, m+1]
  beta <- diag(1/sigma2_X, m) %*% (SYX - alpha * C * sigma2_U)
  SYY <- covmat[m+1, m+1]
  sigma2_Y <- SYY - C^2 * sigma2_U - t(beta) %*% diag(sigma2_X, m) %*% beta
  
  if (check){
    sigma2_Y
  } else {
    list(beta=beta, sigma2_Y=sigma2_Y, gamma=as.vector(C-as.vector(beta) %*% alpha))
  }
}

sens_region <- function(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y){
  # Find upper and lower bounds on C = (alpha'beta + gamma). This is proportional to the
  # covariance between U and Y, which is unidentified in this model.
  new_sigma2 <- function(cc) beta_C(cc, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)
  
  lo <- uniroot(new_sigma2, c(-10 * abs(alpha %*% beta), as.vector(alpha %*% beta + gamma)))$root
  hi <- uniroot(new_sigma2, c(10 * abs(alpha %*% beta), as.vector(alpha %*% beta + gamma)))$root
  
  param_lo <- beta_C(lo, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y, check=FALSE)
  param_hi <- beta_C(hi, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y, check=FALSE)
  return(list(lo=param_lo, hi=param_hi, c_int=c(lo, hi)))
}

param_trace <- function(sens, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y){
  c_grid <- seq(sens$c_int[1], sens$c_int[2], length.out=50)
  
  params <- function(cc) {
    struct <- beta_C(cc, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y, check=FALSE)
    full <- full_covmat(alpha, struct$beta, struct$gamma, sigma2_U, sigma2_X, struct$sigma2_Y)
    part <- partial_cov(full, 0, both=TRUE)
    c(struct$sigma2_Y, part$UY[1,2], part$UYdotX[1,2], full[nrow(full), ncol(full)], struct$beta[1], struct$gamma)
  } 
  
  rbind(c_grid, sapply(c_grid, params))
}

new_gamma <- function(C, alpha, beta){
  as.vector(C - as.vector(beta) %*% alpha)
}

full_covmat <- function(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y){
  UU_theory <- sigma2_U
  XU_theory <- alpha * sigma2_U
  YU_theory <- (as.vector(beta) %*% alpha + gamma) * sigma2_U
  
  Sigma_theory <- rbind(
    c(UU_theory, XU_theory, YU_theory),
    cbind(c(XU_theory, YU_theory),
          make_covmat(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)))
  Sigma_theory
}


## Script doing full sensitivity analysis
sr <- sens_region(alpha_, beta_, gamma_, sigma2_U_, sigma2_X_, sigma2_Y_)
pt <- param_trace(sr, alpha_, beta_, gamma_, sigma2_U_, sigma2_X_, sigma2_Y_)

par(mfrow=c(1,3))
matplot(pt[1,], t(pt[2:5,]))
plot(pt[1,], t(pt[6,]))
plot(pt[1,], t(pt[7,]))

stop()


normal_ll_proxy <- function(par, dat, numproxies, alpha_prior=true_alpha[1]/3,
                            reg=10, deets=FALSE){
  m <- dim(dat)[2] - 1
  N <- dim(dat)[1]
  
  alpha <- rep(par[1], m)
  beta <- c(rep(par[2], m-numproxies), rep(0, numproxies))
  gamma <- par[3]
  sigma2_U <- exp(par[4])
  sigma2_X <- rep(exp(par[5]), m)
  sigma2_Y <- exp(par[6])
  
  cmat <- make_covmat(alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)
  # browser()
  
  if(deets){
    return(cmat)
  } 
  sum(dmvnorm(dat, rep(0, m+1), cmat, log=TRUE)) - reg*(alpha[1] - alpha_prior)^2
}


display_res <- function(res, numproxies, dat){
  par_names <- c("alpha", "beta", "gamma", "sigma2_U", "sigma2_X", "sigma2_Y")
  display_pars <- c(res$par[1:3], exp(res$par[4:6]))
  true_pars <- c(alpha_[1], beta_[1], gamma_, sigma2_U_, sigma2_X_[1], sigma2_Y_)
  display_pars <- rbind(display_pars, true_pars)
  colnames(display_pars) <- par_names
  print(display_pars %>% round(3))
  print(normal_ll_proxy(res$par, dat, numproxies, deets=TRUE) %>% round(3))
  print(res$value)
}



sim_fit <- function(N, numproxies, silent=TRUE, alpha=alpha_, beta=beta_, gamma=gamma_,
                    sigma2_U=sigma2_U_, sigma2_X=sigma2_X_, sigma2_Y=sigma2_Y_, ...){
  i <- 0
  while(TRUE){
    i <- i+1
    simdat <- gendat2(N, numproxies, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)
    dat <- as.matrix(data.frame(t(simdat$X), simdat$Y))
    #init_par <- c(alpha[1], beta[1], gamma, log(sigma2_U), log(sigma2_X[1]), log(sigma2_Y))
    init_par <- runif(6, min=-1, max=1)
    init_par <- c(alpha[1], beta[1], gamma, sigma2_U, sigma2_X[1], sigma2_Y)
    
    print(dim(normal_ll_proxy(init_par, dat, numproxies, deets=TRUE)))
    
    res <- tryCatch(optim(init_par, function(par) normal_ll_proxy(par, dat, numproxies, ...),
              method="L-BFGS-B", control=list(fnscale=-1, factr=1e3, maxit=1e4)),
              error=function(e){ if(i < 10) NULL else stop(e)})
    if(!is.null(res)) break
  }
  if(!silent){
    display_res(res, numproxies, dat)
  }
  return(res)
}

xx0 <- sim_fit(N, 0, silent=FALSE, reg=10)
xx1 <- sim_fit(N, 1, silent=FALSE, reg=10)
xx2 <- sim_fit(N, 2, silent=FALSE, reg=10)
stop()

ests_withproxy <- function(numproxies, endpoints, nreps=20){
  cc_range <- sort(c(10^seq(endpoints[1], endpoints[2], length.out=3), 1))
  ests <- sapply(cc_range,
                 function(alpha_prior_fac){
                   sapply(1:nreps, function(i){
                       sim_fit(N, numproxies, silent=TRUE,
                               alpha_prior=alpha_[1]*alpha_prior_fac,
                               reg=10)$par[2]
                     })
                 })
#       function(alpha_prior_fac){
#                 optim(init_par,
#                       function(par) {
#                         normal_ll_proxy(par, dat, numproxies,
#                                         alpha_prior=alpha_[1]*alpha_prior_fac)
#                       },
#                       method="L-BFGS-B",
#                       control=list(fnscale=-1, pgtol=1e-18, maxit=1e4))$par[2]
#               })
  list(cc_range=cc_range, ests=ests)
}

ests_noproxy <- ests_withproxy(0, c(-0.5, 1))
ests_proxy1 <- ests_withproxy(1, c(-0.5, 1))

#pdf('normal_estimates.pdf', width=7, height=5)
matplot(ests_noproxy$cc_range,
        cbind(colMeans(ests_noproxy$ests), colMeans(ests_proxy1$ests)) / beta_[1],
        type='l',
        bty='n',
        xlab="Unidentified Latent Scaling Factor c Targeted by Penalty",
        ylab="Effect Multiplier",
        main="Ignorance Region as Effect Multiplier",
        lwd=3)
#abline(h=1, lwd=3)
abline(h=c(0,1), col='gray', lty=c(1,2))
abline(v=1, col='gray')
#dev.off()

stop()


# Case with one proxy
N <- 5e4

# Start with 10 treatments
numproxies <- 1
m <- 10 + numproxies
a_0 <- 2
b_0 <- 1
s2X <- 2
# Scale coefficients so that the variance of Y stays fixed even as m grows.
# Add proxy as a unit with known beta = 0.
alpha <- rep(a_0, m) / sqrt(m)
beta <- c(rep(b_0, m-numproxies), rep(0, numproxies)) / sqrt(m-numproxies)
gamma <- 2

sigma2_U <- true_sigma2_U <- 3
sigma2_X <- rep(s2X, m)
sigma2_Y <- 10

simdat <- gendat2(N, alpha, beta, gamma, sigma2_U, sigma2_X, sigma2_Y)
dat <- as.matrix(data.frame(t(simdat$X), simdat$Y))
#init_par <- c(alpha[1], beta[1], gamma, log(sigma2_U), log(sigma2_X[1]), log(sigma2_Y))
init_par <- runif(6, min=-1, max=1)
#normal_ll(init_par, dat, deets=TRUE)


res <- optim(init_par, function(par) normal_ll_proxy(par, dat), method="L-BFGS-B",
             control=list(fnscale=-1, pgtol=1e-18, maxit=1e4))
display_res(res)

#stop()

# Make plots

ests_noproxy <- sapply(10^seq(-0.7, 1, length.out=15),
               function(reg_center_fac){
                 optim(init_par,
                       function(par) {
                         normal_ll_proxy(par, dat, log_scale=log(true_sigma2_U/reg_center_fac^2))
                       },
                       method="L-BFGS-B",
                       control=list(fnscale=-1, pgtol=1e-18, maxit=1e4))$par[2]
               })


ests_proxy1 <- sapply(10^seq(-0.5, 2, length.out=15),
       function(reg_center_fac){
         optim(init_par,
               function(par) {
                 normal_ll_proxy(par, dat, log_scale=log(true_sigma2_U/reg_center_fac^2)) },
               method="L-BFGS-B",
               control=list(fnscale=-1, pgtol=1e-18, maxit=1e4))$par[2]
       })

matplot(cbind(ests_noproxy, ests_proxy1))
abline(h = beta[1])
