############################################################################

### Code for testing whether a predictor is differentially related to two
### outcomes in typical ESM studies (or more generally, in longitudinal
### research, where the two outcomes are repeatedly measured).
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

res <- lme(y ~ outcome*sex - 1,
           random = ~ outcome - 1 | id,
           weights = varIdent(form = ~ 1 | outcome),
           correlation = corSymm(form = ~ outnum | id/beep),
           data = dat, na.action = na.omit,
           control = list(msVerbose=TRUE, maxIter=10000,
                          msMaxIter=10000, msMaxEval=10000))
summary(res)

############################################################################
