############################################################################

### Code for assessing the correlation between two outcomes in typical ESM
### studies (or more generally, in longitudinal research, where the two
### outcomes are repeatedly measured).
###
### Written by Wolfgang Viechtbauer (wvb at wvbauer.com)
### Code is licensed under the following license:
### Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
### See: https://creativecommons.org/licenses/by-sa/4.0/

############################################################################

rm(list=ls())

library(nlme)

### load example data
dat <- read.table("data.dat", header=TRUE, sep="\t")

### calculate PA and NA
dat$pa <- rowMeans(dat[, grepl("pa", names(dat))])
dat$na <- rowMeans(dat[, grepl("na", names(dat))])

### keep only variables that are needed
dat <- dat[, c("id", "sex", "beep", "pa", "na")]

### change into very long format
dat <- reshape(dat, direction="long", varying=c("pa", "na"), v.names="y", idvar="obs", timevar="outcome")
dat$obs <- NULL
dat <- dat[order(dat$id, dat$beep, dat$outcome),]
rownames(dat) <- 1:nrow(dat)
dat$outnum <- dat$outcome
dat$outcome <- factor(dat$outcome)

############################################################################

### fit multivariate multilevel model

res <- lme(y ~ outcome - 1,
           random = ~ outcome - 1 | id,
           weights = varIdent(form = ~ 1 | outcome),
           correlation = corSymm(form = ~ outnum | id/beep),
           data = dat, na.action = na.omit,
           control = list(msVerbose=TRUE, maxIter=10000,
                          msMaxIter=10000, msMaxEval=10000))
summary(res)

### get CIs for correlations at the person and beep levels with intervals()

intervals(res)

### get CIs for correlations at the person and beep levels manually (using Fisher transformation)

getvals <- function(x, pattern)
   unname(x[grepl(pattern, names(x), fixed=TRUE)])

vals <- coef(res$modelStruct, unconstrained=FALSE)
r <- getvals(vals, "cov") / sqrt(getvals(vals, "var(outcome1)") * getvals(vals, "var(outcome2)"))
c(lower=tanh(atanh(r) - qnorm(.975) * sqrt(1 / (res$dims$ngrps[[1]] - 3))), estimate=r, upper=tanh(atanh(r) + qnorm(.975) * sqrt(1 / (res$dims$ngrps[[1]] - 3))))
r <- getvals(coef(res$modelStruct, unconstrained=FALSE), "corStruct")
c(lower=tanh(atanh(r) - qnorm(.975) * sqrt(1 / (res$dims$N/2 - 3))), estimate=r, upper=tanh(atanh(r) + qnorm(.975) * sqrt(1 / (res$dims$N/2 - 3))))

############################################################################

### compare two groups

### fit multivariate multilevel model in each group

levels(dat$sex) ### below i=1 for females, i=2 for males

sav <- list()

for (i in 1:2) {

   sav[[i]] <- lme(y ~ outcome - 1,
                   random = ~ outcome - 1 | id,
                   weights = varIdent(form = ~ 1 | outcome),
                   correlation = corSymm(form = ~ outnum | id/beep),
                   data = dat, na.action = na.omit,
                   subset = c(sex == levels(sex)[i]),
                   control = list(msVerbose=TRUE, maxIter=10000,
                                  msMaxIter=10000, msMaxEval=10000))

}

summary(sav[[1]])
summary(sav[[2]])

vals1 <- coef(sav[[1]]$modelStruct, unconstrained=FALSE)
vals2 <- coef(sav[[2]]$modelStruct, unconstrained=FALSE)

### test for significant difference between correlations at the person level (using Fisher transformation)

r1 <- getvals(vals1, "cov") / sqrt(getvals(vals1, "var(outcome1)") * getvals(vals1, "var(outcome2)"))
r2 <- getvals(vals2, "cov") / sqrt(getvals(vals2, "var(outcome1)") * getvals(vals2, "var(outcome2)"))
z <- (atanh(r1) - atanh(r2)) / sqrt((1 / (sav[[1]]$dims$ngrps[[1]] - 3)) + (1 / (sav[[2]]$dims$ngrps[[1]] - 3)))
z
2*pnorm(abs(z), lower.tail=FALSE)

### test for significant difference between correlations at the beep level (using Fisher transformation)

r1 <- getvals(coef(sav[[1]]$modelStruct, unconstrained=FALSE), "corStruct")
r2 <- getvals(coef(sav[[2]]$modelStruct, unconstrained=FALSE), "corStruct")
z <- (atanh(r1) - atanh(r2)) / sqrt((1 / (sav[[1]]$dims$N/2 - 3)) + (1 / (sav[[2]]$dims$N/2 - 3)))
z
2*pnorm(abs(z), lower.tail=FALSE)

############################################################################
