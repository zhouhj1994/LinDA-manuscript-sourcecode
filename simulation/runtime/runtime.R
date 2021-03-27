setwd("/scratch/user/zhouhuijuan/LinDA/simulation/code")
source('simulate_data.R')
source('competing_methods.R')
para <- readRDS("log.normal.para.rds")
para5000 <- readRDS("log.normal.para5000.rds")
setwd("/scratch/user/zhouhuijuan/LinDA/simulation/runtime")

fun <- function(m, n, case, formula) {
  data <- simulate.data(m = m, n = n, gamma = 0.05, mu = 1.05, beta0 = beta0, 
                        sigma2 = sigma2, eta0 = NULL, S1.zero.prop = NULL, 
                        model = 'S0', case = case)
  otu.tab <- data$Y
  meta <- data$Z
  H <- data$H
  
  t1 <- system.time(rej1 <- linda.fun(otu.tab, meta, formula, alpha))
  print(t1)
  print(sum(H[rej1] == 0) / length(rej1))
  print(sum(H[rej1] == 1) / sum(H)) 
  
  t2 <- system.time(rej2 <- ancombc.fun(otu.tab, meta, formula, alpha))
  print(t2)
  print(sum(H[rej2] == 0) / length(rej2))
  print(sum(H[rej2] == 1) / sum(H)) 
  return(c(t1[3], t2[3]))
}

alpha <- 0.05
sink('output.txt', split = TRUE)
library(benchmarkme)
get_ram()
get_cpu()

beta0 <- para$beta0
sigma2 <- para$sigma2 

set.seed(666)
res11 <- fun(m = 500, n = 200, case = 'C0', formula = 'u') 
set.seed(666)
res12 <- fun(m = 500, n = 200, case = 'C1', formula = 'u') 
set.seed(666)
res13 <- fun(m = 500, n = 200, case = 'C2', formula = 'u+z1+z2') 
set.seed(666)
res21 <- fun(m = 500, n = 10000, case = 'C0', formula = 'u') 
set.seed(666)
res22 <- fun(m = 500, n = 10000, case = 'C1', formula = 'u') 
set.seed(666)
res23 <- fun(m = 500, n = 10000, case = 'C2', formula = 'u+z1+z2') 

beta0 <- para5000$beta0
sigma2 <- para5000$sigma2 

set.seed(666)
res31 <- fun(m = 5000, n = 200, case = 'C0', formula = 'u') 
set.seed(666)
res32 <- fun(m = 5000, n = 200, case = 'C1', formula = 'u') 
set.seed(666)
res33 <- fun(m = 5000, n = 200, case = 'C2', formula = 'u+z1+z2') 
set.seed(666)
res41 <- fun(m = 5000, n = 10000, case = 'C0', formula = 'u') 
set.seed(666)
res42 <- fun(m = 5000, n = 10000, case = 'C1', formula = 'u') 
set.seed(666)
res43 <- fun(m = 5000, n = 10000, case = 'C2', formula = 'u+z1+z2') 
sink()

res <- rbind(c(res11, res12, res13),
             c(res21, res22, res23),
             c(res31, res32, res33),
             c(res41, res42, res43))

write.table(res, 'runtime.txt')

#######################################
# library(stargazer)
# res <- as.matrix(read.table('runtime.txt'))
# stargazer(res, single.row = TRUE)



