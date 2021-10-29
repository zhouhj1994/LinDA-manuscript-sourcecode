para.fun <- function(otu.tab, m) {
  otu.tab <- as.matrix(otu.tab)
  colnames(otu.tab) = NULL
  rownames(otu.tab) = NULL
  
  has.read <- rowSums(otu.tab > 0)
  ind.taxa <- sort(order(has.read, decreasing = TRUE)[1 : m])
  otu.tab.sel <- otu.tab[ind.taxa, ]
  
  theta0 <- kappa0 <- rep(NA, m)
  N <- colSums(otu.tab.sel)
  for(i in 1 : m) {
    fit <- glm.nb(otu.tab.sel[i,] ~ N - 1)
    theta0[i] <- fit$theta
    kappa0[i] <- as.vector(coef(fit))
  }
  out <- list(otu.tab.sel = otu.tab.sel, theta0 = theta0, kappa0 = kappa0)
  return(out)
}

load("combo.RData")
library(MASS)

m <- 500
res <- para.fun(otu.tab, m) 
negBinom.para <- list(theta0 = res[[2]], kappa0 = res[[3]])
saveRDS(negBinom.para, "negBinom.para.rds")
