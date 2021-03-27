run.fun <- function(m, n, gamma, mu,  
                    beta0 = NULL, sigma2 = NULL, eta0 = NULL, 
                    S1.zero.prop = NULL, 
                    model = c('S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'), 
                    case = c('C0', 'C1', 'C2'), formula, alpha, n.met) {
  data <- simulate.data(m, n, gamma, mu, beta0, sigma2, eta0, S1.zero.prop, model, case)
  Y <- data$Y
  Z <- data$Z
  H <- data$H
  
  rej.linda <- linda.fun(Y, Z, formula, alpha)
  
  err <- try(rej.ancombc <- ancombc.fun(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.ancombc <- NA
  
  err <- try(rej.ancombc.2 <- ancombc.fun2(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.ancombc.2 <- NA
  
  err <- try(rej.aldex2 <- aldex2.fun(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.aldex2 <- NA
  
  err <- try(rej.deseq2 <- deseq2.fun(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.deseq2 <- NA
  
  err <- try(rej.edgeR <- edgeR.fun(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.edgeR <- NA
  
  err <- try(rej.metagenomeSeq <- metagenomeSeq.fun(Y, Z, formula, alpha), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.metagenomeSeq <- NA
  
  err <- try(rej.metagenomeSeq.2 <- metagenomeSeq.fun2(Y, Z, formula, alpha), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.metagenomeSeq.2 <- NA
  
  err <- try(rej.wilcox <- wilco.spear.fun(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.wilcox <- NA
  
  err <- try(rej.wilcox.2 <- wilco.spear.fun2(Y, Z, formula, alpha), silent = TRUE)
  if(class(err) == 'try-error') rej.wilcox.2 <- NA
  
  rej.list <- list(rej.linda, rej.ancombc, rej.ancombc.2,
                   rej.aldex2, rej.deseq2, rej.edgeR,
                   rej.metagenomeSeq, rej.metagenomeSeq.2, rej.wilcox, rej.wilcox.2)
  
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
source("competing_methods.R")
para0 <- readRDS("log.normal.para.rds")
para1 <- readRDS("Gamma.para.rds")

setwd("/scratch/user/zhouhuijuan/LinDA/simulation/result")

beta0 <- para0$beta0
sigma2 <- para0$sigma2 
eta0 <- para1$eta0

n.met <- 10

if(case == 'C2') {
  formula <- 'u+z1+z2'
} else {
  formula <- 'u'
}

if (model == 'S5') {
  sample.size.vec <- c(20, 30)
} else {
  sample.size.vec <- c(50, 200)
}

if (model == 'S4') {
  m <- 50
  nsim <- 1000
} else {
  m <- 500
  nsim <- 100
}

sig.density.vec <- c(0.05, 0.2)

sig.strength.vec <- seq(1.05, 2, length.out = 6)

s1 <- 2; s2 <- 2; s3 <- 6

sample.size <- rep(sample.size.vec, each = s2 * s3)
sig.density <- rep(rep(sig.density.vec, each = s3), s1)
sig.strength <- rep(sig.strength.vec, s1 * s2)

setting <- cbind(sample.size, sig.density, sig.strength)

s <- s1 * s2 * s3 

S1.zero.prop <- 0.3

###############################################################
## Simulation runs
###############################################################
library(doSNOW)
cl <- makeCluster(20, type = "SOCK") 
registerDoSNOW(cl)

output.raw <- foreach (i = 1 : nsim,.combine = 'rbind') %dopar% {
  library(LinDA)
  library(ANCOMBC)
  library(phyloseq)
  library(ALDEx2)
  library(GMPR)
  library(DESeq2)
  library(edgeR)
  library(metagenomeSeq)
  
  result <- foreach(j = 1 : s, .combine = 'rbind') %do% {
    para <- setting[j, ]
    n <- para[1]
    gamma <- para[2]
    mu <- para[3]
    set.seed((i - 1) * s + j)
    run.fun(m, n, gamma, mu, 
            beta0 = beta0, sigma2 = sigma2, eta0 = eta0, 
            S1.zero.prop = S1.zero.prop,
            model = model, case = case,
            formula = formula, alpha = 0.05, n.met = n.met) 
  }
  write.table(result, paste0("output_raw_", model, case, "_", i, ".txt"))
  result
}

write.table(output.raw, paste0("output_raw_", model, case, ".txt"))

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

write.table(output, paste0("output_", model, case, ".txt"))

stopCluster(cl)

