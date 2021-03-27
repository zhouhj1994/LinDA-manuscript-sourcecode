run.fun <- function(m, n, gamma, mu,  
                    beta0 = NULL, sigma2 = NULL, eta0 = NULL, 
                    S1.zero.prop = NULL, 
                    model = c('S7.1', 'S7.2'), 
                    case = 'C0', alpha, n.met) {
  data <- simulate.data(m, n, gamma, mu, beta0, sigma2, eta0, S1.zero.prop, model, case)
  Y <- data$Y
  Z <- data$Z
  H <- data$H
  
  res <- linda(Y, Z, '~u')
  pval.ols.bc <- res$output[[1]]$pvalue
  
  lfc <- res$output[[1]]$log2FoldChange
  bias <- res$bias[1]
  lfcSE <- res$output[[1]]$lfcSE
  df <- res$output[[1]]$df
  
  stat <- (lfc + bias) / lfcSE
  pval.ols <- 2 * pt(-abs(stat), df)
  
  res <- linda(Y, Z, '~u+(1|id)')
  pval.lmm.bc <- res$output[[1]]$pvalue
  
  lfc <- res$output[[1]]$log2FoldChange
  bias <- res$bias[1]
  lfcSE <- res$output[[1]]$lfcSE
  df <- res$output[[1]]$df
  
  stat <- (lfc + bias) / lfcSE
  pval.lmm <- 2 * pt(-abs(stat), df)
  
  pval.mat <- cbind(pval.ols.bc, pval.ols, pval.lmm.bc, pval.lmm)
  rej.list <- sapply(1 : 4, function(i){
    qval <- p.adjust(pval.mat[, i], method = 'BH')
    which(qval <= alpha)
  })

  rej.list[[n.met + 1]] <- 10086
  fdp <- rep(NA, n.met)
  power <- rep(NA, n.met)
  
  if (sum(H) == 0) {
    for (i in 1 : n.met) {
      rej <- rej.list[[i]]
      if(length(rej) == 0) {
        fdp[i] <- 0
      } else if (!is.na(rej)){
        fdp[i] <- 1
      }
    }
  } else {
    for (i in 1 : n.met) {
      rej <- rej.list[[i]]
      if(length(rej) == 0) {
        fdp[i] <- 0
        power[i] <- 0
      } else if(!is.na(rej)){
        fdp[i] <- sum(H[rej] == 0) / length(rej)
        power[i] <- sum(H[rej] == 1) / sum(H)
      }
    }
  }
  return(c(fdp, power))
}

###############################################################
## Simulation setups
###############################################################
source("simulate_data.R")
para0 <- readRDS("log.normal.para.rds")

setwd("/scratch/user/zhouhuijuan/LinDA/simulation/result")

beta0 <- para0$beta0
sigma2 <- para0$sigma2 

n.met <- 4

sample.size.vec <- c(50, 200)

m <- 500
nsim <- 100

sig.density.vec <- c(0.05, 0.2)

sig.strength.vec <- seq(1.05, 2, length.out = 6)

s1 <- 2; s2 <- 2; s3 <- 6

sample.size <- rep(sample.size.vec, each = s2 * s3)
sig.density <- rep(rep(sig.density.vec, each = s3), s1)
sig.strength <- rep(sig.strength.vec, s1 * s2)

setting <- cbind(sample.size, sig.density, sig.strength)

s <- s1 * s2 * s3 

###############################################################
## Simulation runs
###############################################################
library(doSNOW)
cl <- makeCluster(20, type = "SOCK") 
registerDoSNOW(cl)

output.raw <- foreach (i = 1 : nsim,.combine = 'rbind') %dopar% {
  library(LinDA)
  result <- foreach(j = 1 : s, .combine = 'rbind') %do% {
    para <- setting[j, ]
    n <- para[1]
    gamma <- para[2]
    mu <- para[3]
    set.seed((i - 1) * s + j)
    run.fun(m, n, gamma, mu, 
            beta0 = beta0, sigma2 = sigma2,
            model = model, alpha = 0.05, n.met = n.met) 
  }
  write.table(result, paste0("output_raw_", model, "C0_", i, ".txt"))
  result
}

write.table(output.raw, paste0("output_raw_", model, "C0.txt"))

output <- foreach (i = 1 : s, .combine = 'rbind') %do% {
  res <- output.raw[seq(i, nsim * s, s), ]
  est <- colMeans(res, na.rm = TRUE)
  num.nan <- colSums(matrix(!is.na(res), nrow = nsim))
  est.sd <- sqrt(rowSums((t(res) - est) ^ 2, na.rm = TRUE) / ((num.nan - 1) * num.nan))
  
  fdr.est <- est[1 : n.met]
  fdr.est.sd <- est.sd[1 : n.met]
  power.est <- est[(n.met + 1) : (2 * n.met)]
  power.est.sd <- est.sd[(n.met + 1) : (2 * n.met)]
  
  res <- c(fdr.est, power.est, fdr.est.sd, power.est.sd)
  res[which(is.nan(res))] <- NA
  res
}

write.table(output, paste0("output_", model, "C0.txt"))

stopCluster(cl)

