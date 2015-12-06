Stat 300
A Design-Based Perspective on Variable Selection
========================================================
author: Alex D'Amour
date: November 25, 2014
incremental: true

Perspectives on Model Selection
========================================================
**Why?**
* Stealth regularization (AIC, CV, LASSO)
* Legitimate regularization (nonparametric regression)
* Infer "relevant" parameters (vaguely causal?)

**This talk**: Find a model that provides the best predictive performance for our given sample size. Note that predictive performance includes estimation uncertainty, bias, and residual variation.

**Holy grail**: Eliminate high-dimensional nuisance without high-dimensional priors.

Perspectives on Model Selection
========================================================
**Is the truth...**
* Finite dimensional and sparse? (Donoho et al. Compressed sensing, fundamentally parametric.)
* Infinite dimensional and dense? (Meng 2014, Everything is variation, fundamentally nonparametric.)

**This talk**: The latter. There may exist a sparse set of predictors, but no reason to believe that the predictors as collected define the proper basis for such sparsity.

<small>**More reading**: Liu and Yang, 2009. "Parametric or nonparametric? A parametricness index for model selection."</small>

Aside: What is truth?
========================================================
**Example**
Let $X_i$  be a $p$-dimensional multivariate normal with covariance matrix $\Sigma$ defined so that $\Sigma_{k,l} = \rho^{|k-l|}$, $0 < \rho < 1$.

Consider:
$$
Y_i \sim X_{i,2} - \rho X_{i,1} + 0.2 X_{i,p} + \mathcal N(0,3).
$$

"True"" model includes covariates $(1,2,p)$.
But for any subset $A \subset \{1,\cdots ,p\}$,
$$
Y_i | X_{i,A} \sim \beta_A^{\top} X_{i,A}\mathcal + N(0, \sigma_A).
$$
because $(Y_i,X_i)$ are jointly multivariate normal. 

"Truth" only has special status because it has minimal residual variance.

Aside: What is truth?
========================================================

Simulation: $N = 100$, $p=25$, $\rho = 0.75$.

For simplicity, consider only growing models $A_k = \{1, \cdots, k\}$.





```
Error in library(mvtnorm) : there is no package called 'mvtnorm'
```
