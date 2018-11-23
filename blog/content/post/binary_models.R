# Utilities

logistic <- function(x) 1 / (1 + exp(-x))
logit <- function(p) p / (1 - p)
odds2prob <- function(o) o / (1 + o)

# Global values

m_ <- 3
alpha0_ <- -1
alpha1_ <- 2
beta0_ <- -2
beta1_ <- 0.5
gamma0_ <- 0.2
piU_ <- 0.3

# Conditional distributions

pAk_U <- function(U){
  logistic(alpha0_ + alpha1_ * U)
}

oddsU_A <- function(A){
  piU_ / (1-piU_) * (pAk_U(1) / pAk_U(0))^sum(A) * ((1-pAk_U(1)) / (1-pAk_U(0)))^(m_ - sum(A))
}

pU_A <- function(A){
  odds2prob(oddsU_A(A))
}

pY_AU <- function(A, U){
  logistic(beta0_ + sum(beta1_ * A) + gamma0_ * U)
}

pY_A <- function(A){
  pU_A(A) * pY_AU(A, 1) + (1-pU_A(A)) * pY_AU(A, 0)
}

pY_doA <- function(A, condfun=pY_AU){
  piU_ * condfun(A, 1) + (1-piU_) * condfun(A, 0)
}




# Plot conditional distributions of U given A
A_classes <- list(c(0,0,0), c(0,0,1), c(0,1,1), c(1,1,1))
pU_A_vec <- sapply(A_classes, pU_A)
pY_A_vec <- sapply(A_classes, pY_A)
pY_doA_vec <- sapply(A_classes, pY_doA)

make_row <- function(p){
  c(1-p, p)
}

true_tabs <- function(A_classes){
 sapply(A_classes, function(A){ 
    row1prop <- make_row(pY_AU(A, 0))
    row2prop <- make_row(pY_AU(A, 1))
    rbind(row1prop * (1-pU_A(A)), row2prop * (pU_A(A)))
  })
}

plot(pU_A_vec)
points(pY_A_vec, col='red')

# Contingency table utilities

p11_bound <- function(pi1, pi2){
  c(max(0, pi1 + pi2 - 1), min(pi1, pi2))
}

rho <- function(p11, pi1, pi2) (p11 - pi1 * pi2) / sqrt(pi1 * (1-pi1) * pi2 * (1-pi2))

ctable_fill <- function(pi1, pi2, p11){
  p10 <- pi1 - p11
  p01 <- pi2 - p11
  p00 <- 1 - (p11 + p10 + p01)
  matrix(c(p00, p10, p01, p11), nr=2)
}

make_condfun <- function(piU, piY, p11){
  ct <- ctable_fill(piU, piY, p11)
  normalizer <- rowSums(ct)
  function(A, U){
    ct[U+1,2] / normalizer[U+1]
  }
}

# For each p(U, Y | A), calculate bounds on p11

num_sens <- 3
pY_doA_sensmat <- matrix(nc=0, nr=num_sens)
for(A_class in A_classes){
  piU <- pU_A(A_class)
  piY <- pY_A(A_class)
  pUY_bound <- p11_bound(piU, piY)
  print(pUY_bound)
  
  p11_seq <- seq(pUY_bound[1], pUY_bound[2], length.out=num_sens)
  tabs <- sapply(p11_seq, function(p11) ctable_fill(piU, piY, p11))
  rhos <- sapply(p11_seq, function(p11) rho(piU, piY, p11))
  print(rhos)
  pY_doA_sweep <- sapply(p11_seq, function(p11) pY_doA(A_class, make_condfun(piU, piY, p11)))
  
  pY_doA_sensmat <- cbind(pY_doA_sensmat, pY_doA_sweep)
  
  #mincor_condfun <- function(A, U) ctable_fill(piU, piY, pUY_bound[1])[U+1,2]
  #maxcor_condfun <- function(A, U) ctable_fill(piU, piY, pUY_bound[2])[U+1,2]
  #pY_doA_mincor <- c(pY_doA_mincor, pY_doA(A_class, mincor_condfun))
  #pY_doA_maxcor <- c(pY_doA_maxcor, pY_doA(A_class, maxcor_condfun))
}

matplot(t(rbind(pY_doA_sensmat, pY_doA_vec)))
#points(pY_doA_vec, pch=15, cex=2)


stop()

pp <- pY %*% t(pU)

print(rowSums(pp))
print(colSums(pp))

# Joe et al p 210
# max{0, p1 + p2 - 1} \leq p11 \leq min{p1, p2}

p11_seq <- seq(0, 0.2, length.out=20)



tabs <- sapply(p11_seq, function(pp11) ctable_fill(pY[2], pU[2], pp11))
rhos <- rho(p11_seq, pY[2], pU[2])
plot(p11_seq, rhos)

indep_rat <- (pY[2] * pU[2]) / (pY[1] * pU[1])
p_min <- min(pY[2], pU[2])
p_max <- max(pY[2], pU[2])

lo_cor <- max(-sqrt(indep_rat), -sqrt(1/indep_rat))
hi_cor <- sqrt((p_min * (1-p_max)) / (p_max * (1-p_min)))

print(c(lo_cor, hi_cor))