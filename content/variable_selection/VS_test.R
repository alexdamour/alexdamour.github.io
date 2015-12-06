library(mvtnorm)

N <- 1e2
p <- 25 
rho <- 0.75
#covmat <- rho^as.matrix(dist(1:p))
covmat <- diag(1,p)
X <- rmvnorm(N, rep(0,p), covmat)
sigma <- 3

E_Y1 <- X[,2] - rho*X[,1] + 0.2*X[,p]
Y1 <- X[,2] - rho*X[,1] + 0.2*X[,p] + rnorm(N, 0, sqrt(sigma))
covmat0 <- sapply(1:p, function(i){ covmat[i,2]-rho*covmat[i,1] + 0.2*covmat[i,p] })

covmat_full <- rbind(covmat0, covmat)
covmat_full <- cbind(c(1,covmat0),covmat_full)

numreps <- 100
statmat <- matrix(NA, nc=5, nr=0)
y_MSE_obs_mat <- matrix(NA, nc=numreps, nr=N)
for(i in 1:p){
    X_cur <- X[,1:i]
    proj_inv <- solve(covmat[1:i,1:i])
    out_proj <- covmat0[1:i]

    for(r in 1:numreps){
        y_hat <- out_proj%*%proj_inv%*%t(X_cur)
        y_var <- as.vector(2*sigma-out_proj%*%proj_inv%*%out_proj)
        y_bias2 <- (y_hat-E_Y1)^2
        y_MSE <- mean(y_var + y_est_var)
        y_SE_obs <- (y_hat-Y1)^2
        #print(y_MSE)
        y_SE_obs_mat[,r] <- y_SE_obs
    }

    y_est_var <- diag(X_cur%*%(solve(t(X_cur)%*%X_cur)*sigma)%*%t(X_cur))
    statmat <- rbind(statmat, c(y_var, mean(y_bias2), mean(y_est_var), y_MSE, y_MSE_obs))
}

par(mfrow=c(2,2))
titles <- c("Resid Variance", "Bias^2", "Estimation Variance", "MSE")
for(i in 1:4){
    plot(statmat[,i], type='l', main=titles[i])
    if(i ==4)
        lines(statmat[,5], col='red')
}

m1 <- lm(Y1 ~ X -1)
