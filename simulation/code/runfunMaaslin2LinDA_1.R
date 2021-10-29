run.fun <- function(m, n, gamma, mu,  
                    beta0 = NULL, sigma2 = NULL, eta0 = NULL, S1.zero.prop = NULL, 
                    model = c('S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'), 
                    case = c('C0', 'C1', 'C2'), strong,
                    fix.eff, formula, alpha, n.met) {
  data <- simulate.data(m, n, gamma, mu, beta0, sigma2, eta0, theta0 = NULL, kappa0 = NULL,
                        S1.zero.prop, model, case, strong)
  Y <- data$Y
  Z <- data$Z
  H <- data$H
  
  res <- linda(Y, Z, paste('~', formula))
  pval.ols.bc <- res$output[[1]]$pvalue
  
  lfc <- res$output[[1]]$log2FoldChange
  bias <- res$bias[1]
  lfcSE <- res$output[[1]]$lfcSE
  df <- res$output[[1]]$df
  
  stat <- (lfc + bias) / lfcSE
  pval.ols <- 2 * pt(-abs(stat), df)
  
  qval.ols.bc <- p.adjust(pval.ols.bc, method = 'BH')
  qval.ols <- p.adjust(pval.ols, method = 'BH')
  rej.linda.1 <- which(qval.ols.bc <= alpha)
  rej.linda.2 <- which(qval.ols <= alpha)
  
  maaslin2.fun <- function(otu.tab, meta, fix.eff, alpha, 
                           normalization, transform) {
    res <- Maaslin2(
      input_data = otu.tab,
      input_metadata = meta,
      output = "masslin2", 
      min_prevalence = 0,
      fixed_effects = fix.eff,
      random_effects = NULL,
      normalization = normalization,
      transform = transform,
      plot_heatmap = FALSE,
      plot_scatter = FALSE)
    ind <- which(res$results$metadata == fix.eff[1])
    ord <- match(rownames(otu.tab), res$results$feature[ind])
    pval <- res$results$pval[ind][ord]
    qval <- p.adjust(pval, method = 'BH')
    rej <- which(qval <= alpha)
    return(rej)
  }
  
  err <- try(rej.maaslin2.1 <- maaslin2.fun(Y, Z, fix.eff, alpha,
                                            normalization = 'TSS', transform = 'LOG'), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.maaslin2.1 <- NA
  
  err <- try(rej.maaslin2.2 <- maaslin2.fun(Y, Z, fix.eff, alpha,
                                            normalization = 'TMM', transform = 'LOG'), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.maaslin2.2 <- NA
  
  err <- try(rej.maaslin2.3 <- maaslin2.fun(Y, Z, fix.eff, alpha,
                                            normalization = 'CSS', transform = 'LOG'), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.maaslin2.3 <- NA
  
  err <- try(rej.maaslin2.4 <- maaslin2.fun(Y, Z, fix.eff, alpha,
                                            normalization = 'CLR', transform = 'NONE'), 
             silent = TRUE)
  if(class(err) == 'try-error') rej.maaslin2.4 <- NA
  
  rej.list <- list(rej.linda.1, rej.linda.2, rej.maaslin2.1, rej.maaslin2.2,
                   rej.maaslin2.3, rej.maaslin2.4)
  
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
para1 <- readRDS("Gamma.para.rds")

beta0 <- para0$beta0
sigma2 <- para0$sigma2 
eta0 <- para1$eta0

n.met <- 6

if(case == 'C2') {
  formula <- 'u+z1+z2'
} else {
  formula <- 'u'
}

if(case == 'C2') {
  fix.eff <- c('u', 'z1', 'z2')
} else {
  fix.eff <- 'u'
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
  library(Maaslin2)
  library(foreach)
  library(Matrix)
  
  result <- foreach(j = 1 : s, .combine = 'rbind') %do% {
    para <- setting[j, ]
    n <- para[1]
    gamma <- para[2]
    mu <- para[3]
    set.seed((i - 1) * s + j)
    run.fun(m, n, gamma, mu, 
            beta0 = beta0, sigma2 = sigma2, eta0 = eta0, 
            S1.zero.prop = S1.zero.prop,
            model = model, case = case, strong = strong, 
            fix.eff = fix.eff, formula = formula, alpha = 0.05, n.met = n.met) 
  }
  if(strong) {
    sink(paste0(model, case,'Strong_progress_1.txt'), append = TRUE)
  } else {
    sink(paste0(model, case,'_progress_1.txt'), append = TRUE) 
  }
  print(i)
  sink()
  result
}

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

if(strong) {
  write.table(output, paste0("output_", model, case, "Strong_Maaslin2LinDA_1", ".txt"))
} else {
  write.table(output, paste0("output_", model, case, "_Maaslin2LinDA_1", ".txt"))
}

stopCluster(cl)

## S0C0, S0C1, S0C2, S0C0Strong, S1C0-S6C0
