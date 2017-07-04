############################################################################

### Code for assessing the factor structure and reliability of items (e.g,
### for measuring positive and negative affect) in typical ESM studies (or
### more generally, in longitudinal research, where subjects repeatedly
### provide ratings to the items).
###
### Written by Wolfgang Viechtbauer (wvb at wvbauer.com)
### Code is licensed under the following license:
### Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
### See: https://creativecommons.org/licenses/by-sa/4.0/

############################################################################

library(nlme)
library(sem)

### load example data (3 items for positive affect, 4 items for negative affect)

dat <- read.table("data.dat", header=TRUE, sep="\t")

### change into very long format

dat <- reshape(dat, direction="long", varying=c("pa1", "pa2", "pa3", "na1", "na2", "na3", "na4"), v.names="y", idvar="obs", timevar="item")
dat$obs <- NULL
dat <- dat[order(dat$id, dat$beep, dat$item),]
rownames(dat) <- 1:nrow(dat)
dat$itemnum <- dat$item
dat$item <- factor(dat$item)

############################################################################

### fit multivariate multilevel model

res <- lme(y ~ item - 1,
           random = ~ item - 1 | id,
           weights = varIdent(form = ~ 1 | item),
           correlation = corSymm(form = ~ itemnum | id/beep),
           data = dat, na.action = na.omit,
           control = list(msVerbose=TRUE, maxIter=10000,
                          msMaxIter=10000, msMaxEval=10000))
summary(res)
#save(res, file="res_fac_rel_esm.rdata")

### function to extract the var-cov matrix at the beep level

extractR <- function(object) {
   mC  <- attr(summary(object$modelStruct)$corStruct, "maxCov")
   Rcor <- diag(mC)
   dimnames(Rcor) <- list(1:mC, 1:mC)
   Rcor[lower.tri(Rcor)] <- nlme:::coef.corSymm(summary(object$modelStruct)$corStruct, FALSE)
   Rcor[upper.tri(Rcor)] <- t(Rcor)[upper.tri(Rcor)]
   x <- summary(object$modelStruct)$varStruct
   class(x) <- attr(x, "oClass")
   s.fact <- coef(x, uncons = FALSE, allCoef = TRUE)
   sds <- object$sigma * s.fact
   R <- diag(sds) %*% Rcor %*% diag(sds)
   return(R)
}

### extract var-cov and correlation matrices at the person and beep levels

G <- getVarCov(res)[1:7,1:7]
rownames(G) <- colnames(G) <- c("pa1", "pa2", "pa3", "na1", "na2", "na3", "na4")
R <- extractR(res)
rownames(R) <- colnames(R) <- c("pa1", "pa2", "pa3", "na1", "na2", "na3", "na4")

Gr <- cov2cor(G)
Rr <- cov2cor(R)

round(Gr, 3)
round(Rr, 3)

### fit 2-factor (CFA) models with Gr and Rr matrices

model <- matrix(c(
'PA -> pa1',   'PA.to.pa1', NA,
'PA -> pa2',   'PA.to.pa2', NA,
'PA -> pa3',   'PA.to.pa3', NA,
'NA -> na1',   'NA.to.na1', NA,
'NA -> na2',   'NA.to.na2', NA,
'NA -> na3',   'NA.to.na3', NA,
'NA -> na4',   'NA.to.na4', NA,
'NA <-> NA',    NA,         1,
'PA <-> PA',    NA,         1,
'NA <-> PA',   'cor.PANA',  NA,
'pa1 <-> pa1', 's2.pa1',    NA,
'pa2 <-> pa2', 's2.pa2',    NA,
'pa3 <-> pa3', 's2.pa3',    NA,
'na1 <-> na1', 's2.na1',    NA,
'na2 <-> na2', 's2.na2',    NA,
'na3 <-> na3', 's2.na3',    NA,
'na4 <-> na4', 's2.na4',    NA),
ncol=3, byrow=TRUE)

res.sem.pers <- sem(model, Gr, N=res$dims$ngrps[1])
res.sem.beep <- sem(model, Rr, N=res$dims$N, optimizer=optimizerNlminb)

summary(res.sem.pers, fit.indices=c("RMSEA","CFI","NNFI","SRMR"))
summary(res.sem.beep, fit.indices=c("RMSEA","CFI","NNFI","SRMR"))

### compute reliabilities of the PA and NA scales at the person and beep levels

getvals <- function(x, pattern)
   unname(x[grepl(pattern, names(x), fixed=TRUE)])

rel.pa.pers <- sum(getvals(coef(res.sem.pers), "PA.to"))^2 / (sum(getvals(coef(res.sem.pers), "PA.to"))^2 + sum(getvals(coef(res.sem.pers), "s2.pa")))
rel.na.pers <- sum(getvals(coef(res.sem.pers), "NA.to"))^2 / (sum(getvals(coef(res.sem.pers), "NA.to"))^2 + sum(getvals(coef(res.sem.pers), "s2.na")))
rel.pa.beep <- sum(getvals(coef(res.sem.beep), "PA.to"))^2 / (sum(getvals(coef(res.sem.beep), "PA.to"))^2 + sum(getvals(coef(res.sem.beep), "s2.pa")))
rel.na.beep <- sum(getvals(coef(res.sem.beep), "NA.to"))^2 / (sum(getvals(coef(res.sem.beep), "NA.to"))^2 + sum(getvals(coef(res.sem.beep), "s2.na")))

cat("Rel PA Person:", round(rel.pa.pers, 2), "\n")
cat("Rel NA Person:", round(rel.na.pers, 2), "\n")
cat("Rel PA Beep:  ", round(rel.pa.beep, 2), "\n")
cat("Rel NA Beep:  ", round(rel.na.beep, 2), "\n")

############################################################################
