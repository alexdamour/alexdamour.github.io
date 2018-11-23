WmidU <- matrix(c(0.3, 0.7, 0.4, 0.6), nr=2)

U <- c(0.25, 0.75)

W <- WmidU %*% U

solve(WmidU)

solve(WmidU) %*% W

BRmat <- t( WmidU / as.vector(W)) * as.vector(U)

# WmidU is the pushforward operator

# Example 1: Using the moment formula from Miao et al. in a
# binary treatment, binary U, normal emission example.

db <- function(var, prob) prob^var * (1-prob)^(1-var)
  
logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

run_proxy_sim <- function(){
  N <- 1e5

  pU <- 0.5
  beta_X <- 1
  beta_Z <- -2
  beta_W <- 1
  
  mean_mat <- matrix(c(-2, -1, 1, 2), nr=2) * 4
  
  pX_U <- function(U) logistic(logit(pU) + U*beta_X)
  pZ_U <- function(U) logistic(logit(pU) + U*beta_Z)
  pW_U <- function(U) logistic(logit(pU) + U*beta_W)
  
  U <- rbinom(N, 1, pU)
  X <- rbinom(N, 1, logistic(pX_U(U)))
  Z <- rbinom(N, 1, logistic(PZ_U(U)))
  W <- rbinom(N, 1, logistic(PW_U(U)))
  
  mu_Y_factory <- function(mean_mat){
    # number of levels of each variable
    nU <- 2
    nX <- 2
    function(U, X){
      # Set up x-dependence as columns, u-dependence as rows
      idx <- U * nU + X + 1
      mean_mat[idx]
    }
  }
  mu_Y <- mu_Y_factory(mean_mat)
  
  Y <- rnorm(N, mu_Y(U, X))
  
  miao_formula <- function(X, Z, W, Y, x){
    df <- data.frame(x=X, z=Z, w=W, y=Y)
    wz_x_counts <- table(df[,c('w', 'z', 'x')])
    pw_zx <- apply(wz_x_counts, 3, function(x) t(t(x) / rowSums(t(x))))
    WZmat <- matrix(pw_zx[,as.character(x)], nr=2)
    Wvec <- matrix(table(W) / length(W), nc=1)
    solve(WZmat) %*% Wvec
  }
  
  miao_formula_analytical <- function(pU, pX_U, pZ_U, pW_U, x){
    pzx <- function(z, x){
      db(x, pX_U(0)) * db(z, pZ_U(0)) * (1-pU) +
        db(x, pX_U(1)) * db(z, pZ_U(1)) * pU
    }
    
    pw_zx <- function(w, z, x){
      # P(w | z,x) = sum_u P(w | u) P(u | z, x) = sum_u P(w | u) p(z, x, | u) p(u) / p(z, x)
      db(w, pW_U(0)) * db(z, pZ_U(0)) * db(x, pX_U(0)) * (1-pU) / pzx(z, x) +
        db(w, pW_U(1)) * db(z, pZ_U(1)) * db(x, pX_U(1)) * pU / pzx(z, x)
    }
    
    pW_Zx_mat <- function(x) matrix(pw_zx(c(0, 1, 0, 1), c(0, 0, 1, 1), x), nc=2)
    
    pW_vec <- db(c(0, 1), pU * pW_U(1) + (1-pU) * pW_U(0))
    
    solve(pW_Zx_mat(x)) %*% pW_vec
  }
  
  EY_zx <- function(x) sapply(c(0, 1), function(z) mean(Y[Z == z & X == x]))
  
  ests <- sapply(c(0, 1), function(x)  EY_zx(x) %*% miao_formula(X, Z, W, Y, x))
  ests_ana <- sapply(c(0, 1), function(x) EY_zx(x) %*% miao_formula_analytical(pU, pX_U, pZ_U, pW_U, x))
  truth_sim <- sapply(c(0, 1), function(x) mean(mu_Y(U, x)))
  truth_ana <- sapply(c(0, 1), function(x) pU * mu_Y(1, x) + (1-pU) * mu_Y(0, x))
  
  cbind(ests,
        truth_sim,
        ests_ana,
        truth_ana)
}

#sims <- sapply(1:100, function(x) run_proxy_sim())
#print(rowMeans(sims))

print(run_proxy_sim())

stop()


library(ggplot2)
library(ggridges)
library(magrittr)

uxy_df <- data.frame(u=U, x=X, y=Y, z=Z, w=W)
obs_plot <- ggplot(uxy_df, aes(y, as.factor(x))) + geom_density_ridges(scale=0.9) +
  facet_wrap(~z+w)+
  #geom_density_ridges(
  #  aes(point_color = as.factor(w), point_fill = as.factor(w), point_shape = as.factor(w)),
  #  alpha = .2, point_alpha = 1, jittered_points = TRUE
  #) +
  #scale_point_color_hue(l = 40) +
  #scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23))
  scale_y_discrete(expand = c(0.01, 0))
print(obs_plot)


#hist(Y, breaks=100)

#idxs <- expand.grid(u=uidx, x=xidx)
#abline(v = mu_Y(idxs$u, idxs$x), col=idxs$u + 1, lty = idxs$x + 1, lwd=2)

if(FALSE){
  ## Some other stuff
  
  sep_scale_X <- 2
  sep_scale_Z <- -2
  sep_scale_W <- 4
  alpha <- 2
  beta <- 2
  gamma <- -2
  
  
  U <- rbinom
  
  X_c <- rnorm(N, sep_scale_X * (2*U - 1))
  X <- as.numeric(X_c > 0)
  Z_c <- rnorm(N, sep_scale_Z * (2*U - 1))
  Z <- as.numeric(Z_c > 0)
  W_c <- rnorm(N, sep_scale_W * (2*U - 1))
  W <- as.numeric(W_c > 0)
  
  alpha_v <- c(2, 1)
  gamma_v <- c(2, 1)
  
  # Y <- alpha * X_c + beta * U + gamma * W_c
  Y <- alpha_v[U+1] * X_c + gamma_v[U+1] * W_c
  
  hist(X_c, breaks=100)
  hist(W_c, breaks=100)
  hist(Y, breaks=100)
  
  cov(cbind(X_c, W_c, Z_c))
  par(mfrow=c(4, 1), mar=c(1,1,1,1))
  hist(X_c, breaks=100)
  hist(W_c, breaks=100)
  hist(Z_c, breaks=100)
  hist(Y, breaks=100)
}

