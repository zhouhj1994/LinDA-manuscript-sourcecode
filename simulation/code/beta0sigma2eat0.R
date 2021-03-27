para.fun <- function(otu.tab, m, D = NULL, model = c('log-normal', 'Gamma')) {
  otu.tab <- as.matrix(otu.tab)
  colnames(otu.tab) = NULL
  rownames(otu.tab) = NULL
  
  has.read <- rowSums(otu.tab > 0)
  ind.taxa <- sort(order(has.read, decreasing = TRUE)[1 : m])
  otu.tab.sel <- otu.tab[ind.taxa, ]
  
  if(model == 'log-normal'){
    W <- D %*% log(otu.tab.sel + 0.1)
    n <- ncol(otu.tab.sel)
    
    beta.diff.est <- rowMeans(W)
    sigma2.add.est <- rowSums((W - beta.diff.est) ^ 2) / (n - 1)
    fit <- cv.glmnet(D, beta.diff.est, alpha = 1, standardize = FALSE, intercept = FALSE)
    beta.est <- coef(fit, s = 'lambda.min')[-1]
    
    s <- sum(sigma2.add.est) / (m - 1)
    ind <- t(t(matrix(Matrix::which(D != 0), nrow = m - 1)) - 
               (0 : (m - 1)) * m * (m - 1) /  2)
    tmp <- matrix(sigma2.add.est[ind], nrow = m - 1)
    sigma2.est <- (colSums(tmp) - s) / (m - 2)
    out <- list(otu.tab.sel = otu.tab.sel, beta0 = beta.est, sigma2 = sigma2.est)
  } else if(model == 'Gamma') {
    res <- dirmult(t(otu.tab.sel))
    pi0 <- res$pi
    theta0 <- res$theta
    eta0 <- res$gamma
    out <- list(otu.tab.sel = otu.tab.sel, pi0 = pi0, theta0 = theta0, eta0 = eta0)
  }
  return(out)
}

D.fun <- function (m) {
  D <- foreach(k = 1 : (m - 1), .combine = 'rbind') %do% {
    
    i <- rep(1 : (m - k), 2)
    j <- c(rep(k, m - k), (k + 1) : m)
    x <- c(rep(1, m - k), rep(-1, m - k))
    
    sparseMatrix(i, j, x = x, dims = c(m - k, m))
    
  }
  return(D)
}

#####################################
load("combo.RData")
library(Matrix)
library(foreach)
library(glmnet)
library(dirmult)

m <- 500
D <- D.fun(m)

res <- para.fun(otu.tab, m, D, model = 'log-normal') 
log.normal.para <- list(beta0 = res[[2]], sigma2 = res[[3]])

res <- para.fun(otu.tab, m, model = 'Gamma')
theta0 <- 0.003 
pi0 <- res$pi0
eta0.sum <- 1 / theta0 - 1
eta0 <- pi0 * eta0.sum
Gamma.para <- list(pi0 = pi0, theta0 = theta0, eta0 = eta0)

saveRDS(log.normal.para, "log.normal.para.rds")
saveRDS(Gamma.para, "Gamma.para.rds")

otu.tab.sel <- res$otu.tab.sel
s <- colSums(otu.tab.sel)
ave <- mean(s)
sd2 <- var(s)

## r: The number of successes before terminating.
## p: The probability of a success.
p <- ave / sd2
r <- ave * p / (1 - p)
p; r # 0.0007; 5.3

## use mean and dispersion parameterization
r <- ave ^ 2 / (sd2 - ave)
ave; r # 7645; 5.3

## for runtime m = 5000
m <- 5000
D <- D.fun(m)
 
res <- para.fun(otu.tab, m, D, model = 'log-normal') 
log.normal.para <- list(beta0 = res[[2]], sigma2 = res[[3]])
saveRDS(log.normal.para, "log.normal.para5000.rds")


