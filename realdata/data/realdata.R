load("CDI.RData")
ind <- data.obj$meta.dat$disease_stat %in% c('Case', 'DiarrhealControl')
CDI.otu <- as.data.frame(data.obj$otu.tab[, ind])
CDI.meta <- data.frame(Disease = factor(data.obj$meta.dat$disease_stat[ind]))
rownames(CDI.meta) <- rownames(data.obj$meta.dat)[ind]
dim(CDI.otu)

# case <- factor(data.obj$meta.dat$disease_stat[ind])
# anti <- factor(data.obj$meta.dat$antibiotics[ind])
# fisher.test(table(case, anti))
# age <- data.obj$meta.dat$age[ind]
# summary(glm(case~age+anti, family = 'binomial'))
# 1822 183

load("EmilyIBD.RData")
ind <- meta$ULCERATIVE_COLIT_OR_CROHNS_DIS %in% c("Crohn's disease", "Healthy")
IBD.otu <- as.data.frame(otu[, ind])
IBD.meta <- cbind.data.frame(Disease = factor(meta$ULCERATIVE_COLIT_OR_CROHNS_DIS[ind]),
                             Antibiotic = factor(meta$ANTIBIOTICS[ind]))
rownames(IBD.meta) <- rownames(meta)[ind]

crohn <- factor(meta$ULCERATIVE_COLIT_OR_CROHNS_DIS[ind])
anti <- factor(meta$ANTIBIOTICS[ind])
fisher.test(table(crohn, anti)) #p=0.0.03; OR=0
dim(IBD.otu)
# 2687   82

load("RA_elife.RData")
ind <- data.obj$meta.dat$Disease %in% c('HLT', 'NORA')
RA.otu <- as.data.frame(data.obj$otu.tab[, ind])
RA.meta <- data.frame(Disease = factor(data.obj$meta.dat$Disease[ind]))
rownames(RA.meta) <- rownames(data.obj$meta.dat)[ind]
dim(RA.otu)
# 991  72

load("smoker_qiita_full.RData")
ind <- smokers$meta$AIRWAYSITE == 'Throat'
SMOKE.otu <- as.data.frame(smokers$otu[, ind])
SMOKE.meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
                               Sex = factor(smokers$meta$SEX[ind]),
                               Site = factor(smokers$meta$SIDEOFBODY[ind]),
                               SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
rownames(SMOKE.meta) <- rownames(smokers$meta)[ind]

smoke <- smokers$meta$SMOKER[ind]
sex <- smokers$meta$SEX[ind]
fisher.test(table(smoke, sex)) #p=0.0.02; OR=2.26
dim(SMOKE.otu)
# 2156 145

data.list <- list(CDI.otu = CDI.otu, CDI.meta = CDI.meta,
                  IBD.otu = IBD.otu, IBD.meta = IBD.meta,
                  RA.otu = RA.otu, RA.meta = RA.meta,
                  SMOKE.otu = SMOKE.otu, SMOKE.meta = SMOKE.meta)

saveRDS(data.list, "CDI_IBD_RA_SMOKE.rds")
