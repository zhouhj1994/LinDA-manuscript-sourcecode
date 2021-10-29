simulate.data <- function(m, n, gamma, mu, beta0 = NULL, sigma2 = NULL, eta0 = NULL, 
                          theta0 = NULL, kappa0 = NULL, S1.zero.prop = NULL, 
                          model = c('S0', 'S1', 'S2', 'S3', 'S4', 'S5', 
                                    'S6', 'S7.1', 'S7.2', 'S8'), 
                          case = c('C0', 'C1', 'C2'), strong = FALSE) {
  if (model == 'S4') {
    ind <- sample(1 : 500, m)
    beta0 <- beta0[ind]
    sigma2 <- sigma2[ind]
  }
  if (model %in% c('S0', 'S1', 'S2', 'S4', 'S5', 'S6', 'S7.1', 'S7.2')) {
    X0 <- matrix(exp(rnorm(m * n, beta0, sqrt(sigma2))), nrow = m)
  } else if (model == 'S3'){
    X0 <- matrix(rgamma(m * n, shape = eta0, rate = 1), nrow = m)
  } else if (model == 'S8') {
    X0 <- matrix(rnbinom(m * n, size = theta0, mu = exp(kappa0 * 7645)), nrow = m)
  }
  pi0 <- t(t(X0) / colSums(X0))
  pi0.ave <- rowMeans(pi0)
  if (any(pi0.ave == 0)) {
    ind <- which(pi0.ave == 0)
    pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
    pi0.ave <- pi0.ave / sum(pi0.ave)
  }
  tmp <- (pi0.ave > 0.005)
  mu <- 2 * mu * (n <= 50) + mu * (n > 50)
  mu.1 <- log(mu * tmp + mu * (0.005 / pi0.ave) ^ (1 / 3) * (1 - tmp))
  
  if(strong) {
    H <- rep(0, m)
    abun <- order(pi0.ave, decreasing = TRUE)[1 : floor(m * 0.25)]
    H[sample(abun, m * gamma)] <- 1
  } else {
    H <- rbinom(m, 1, gamma)
  }
  alpha <- mu.1 * H
  
  if (model %in% c('S7.1', 'S7.2')) {
    if (model == 'S7.1') {
      n.id <- n / 2
      u <- rep(0 : 1, n.id)
      id <- rep(1 : n.id, each = 2)
    } else if (model == 'S7.2') {
      if (n == 50) {
        n.id <- 25
        u <- c(rep(0, 24), rep(1, 26))
        id <- rep(1 : n.id, each = 2)
      } else if (n == 200) {
        n.id <- 50
        u <- c(rep(0, 100), rep(1, 100))
        id <- rep(1 : n.id, each = 4)
      }
    }
    tau2 <- runif(m, 0, 1) * sigma2
    r <- matrix(rnorm(m * n.id, 0, sqrt(tau2)), nrow = m)[, id]
    Z <- cbind(u, id)
    beta <- alpha
    tmp <- beta %*% t(Z[, 1]) + beta0
    X <- matrix(exp(rnorm(m * n, tmp, rep(sqrt(sigma2), n))), nrow = m) * exp(r)
  } else {
    if (case == 'C0') {
      u <- rbinom(n, 1, 0.5)
      while(length(unique(u)) == 1) {
        u <- rbinom(n, 1, 0.5)
      }
      Z <- cbind(u)
      beta <- alpha
    } else if (case == 'C1') {
      u <- rnorm(n, 0, 1)
      Z <- cbind(u)
      beta <- alpha
    } else if (case == 'C2') {
      z1 <- rbinom(n, 1, 0.5)
      z1[which(z1 == 0)] <- -1
      z2 <- rnorm(n, 0, 1)
      beta1 <- rnorm(m, 1, 1)
      beta2 <- rnorm(m, 2, 1)
      
      u <- rbinom(n, 1, 1 / (1 + exp(-0.5 * z1 - 0.5 * z2)))
      while (length(unique(u)) == 1) {
        u <- rbinom(n, 1, 1 / (1 + exp(-z1 - z2)))
      }
      Z <- cbind(u, z1, z2)
      beta <- cbind(alpha, beta1, beta2)
    }
  }
  
  if (model %in% c('S0', 'S1', 'S4', 'S5', 'S6')) {
    tmp <- beta %*% t(Z) + beta0
    X <- matrix(exp(rnorm(m * n, tmp, rep(sqrt(sigma2), n))), nrow = m)
  } else if (model == 'S3') {
    tmp <- exp(beta %*% t(Z)) * eta0
    X <- matrix(rgamma(m * n, shape = tmp, rate = 1), nrow = m)
  } else if (model == 'S2') {
    tmp1 <- matrix(0.5, 10, 10)
    diag(tmp1) <- 1
    tmp2 <- matrix(-0.5, 10, 10)
    corr <- rbind(cbind(tmp1, tmp2), cbind(tmp2, tmp1))
    Corr <- list()
    for(i in 1 : 25) Corr[[i]] <- corr
    Sig2 <- bdiag(Corr) * (sqrt(sigma2) %*% t(sqrt(sigma2)))
    tmp <- eigen(Sig2)
    Sig <- tmp$vectors %*% diag(sqrt(tmp$values)) %*% t(tmp$vectors)
    tmp <- beta %*% t(Z) + beta0
    X <- exp(Sig %*% matrix(rnorm(m * n), nrow = m) + tmp)
  } else if (model == 'S8') {
    N <- rnbinom(n, size = 5.3, mu = 7645)
    X <- exp(kappa0 %*% t(N) + beta %*% t(Z))
    Y <- matrix(rnbinom(m * n, size = theta0, mu = X), nrow = m)
  }
  
  pi <- t(t(X) / colSums(X))
  if (model %in% c('S0', 'S2', 'S3', 'S5', 'S7.1', 'S7.2')) {
    N <- rnbinom(n, size = 5.3, mu = 7645)
    Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
  } else if (model == 'S1') {
    ind <- sample(1 : (m * n), m * n * S1.zero.prop)
    Y <- floor(X)
    Y[ind] <- 0
  } else if (model == 'S4') {
    N <- rnbinom(n, size = 5.3, mu = 1500)
    Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
  } else if (model == 'S6') {
    N <- rep(NA, n)
    ind <- which(u == 0)
    N[ind] <- rnbinom(length(ind), size = 5.3, mu = 5000)
    N[-ind] <- rnbinom(n - length(ind), size = 5.3, mu = 50000)
    Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
  }
  
  Z <- as.data.frame(Z)
  if (case %in% c('C0', 'C2')) Z$u <- factor(Z$u)
  if (case == 'C2') Z$z1 <- factor(Z$z1)
  if (model %in% c('S7.1', 'S7.2')) Z$id <- factor(Z$id)
  Y <- as.data.frame(Y)
  colnames(Y) <- rownames(Z) <- paste0('sample', 1 : n)
  rownames(Y) <- paste0('taxon', 1 : m)
  return(list(Y = Y, Z = Z, H = H, X = X, beta = beta))
}
