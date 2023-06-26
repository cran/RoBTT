## ----setup, include=FALSE-----------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/Introduction/fit_3.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev      = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}
print_conditional <- function(fit, coef = "delta"){
  with(fit, 
       sprintf("%1$s, 95%% CI [%2$s, %3$s]",
               format(round(mean(RoBTT$posterior[[coef]]), 2), nsmall = 3),
               format(round(quantile(RoBTT$posterior[[coef]], 0.025), 2), nsmall = 2),
               format(round(quantile(RoBTT$posterior[[coef]], 0.975), 2), nsmall = 2)
       ))
} 

## ----include = FALSE, eval = FALSE--------------------------------------------
#  # pre-fit all models (easier to update the code on package update)
#  library(RoBTT)
#  
#  data("fertilization", package = "RoBTT")
#  
#  fit_01 <- RoBTT(
#    x1       = fertilization$Self,
#    x2       = fertilization$Crossed,
#    parallel = TRUE,
#    prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
#    prior_rho   = NULL, # this indicates no prior on the variance allocation factor -> equal variance test
#    prior_nu    = NULL, # this indicates no prior on the degrees of freedom -> normal distribution test
#    seed        = 0
#  )
#  fit_10 <- RoBTT(
#    x1       = fertilization$Self,
#    x2       = fertilization$Crossed,
#    parallel = TRUE,
#    prior_delta    = prior("cauchy", list(0, 1/sqrt(2))),
#    prior_rho      = prior("beta",   list(3, 3)), #prior on variance allocation
#    prior_rho_null = NULL, # remove models assuming equal variance
#    prior_nu       = NULL,
#    seed           = 0
#  )
#  fit_1 <- RoBTT(
#    x1       = fertilization$Self,
#    x2       = fertilization$Crossed,
#    parallel = TRUE,
#    prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
#    prior_rho   = prior("beta",   list(3, 3)),
#    prior_nu    = NULL,
#    seed        = 1,
#    control     = set_control(adapt_delta = 0.95)
#  )
#  fit_2 <- RoBTT(
#    x1       = fertilization$Self,
#    x2       = fertilization$Crossed,
#    parallel = TRUE,
#    prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
#    prior_rho   = prior("beta",   list(3, 3)),
#    prior_nu    = prior("exp",    list(1)), # prior on degrees of freedom
#    seed        = 2
#  )
#  fit_3 <- RoBTT(
#    x1       = fertilization$Self,
#    x2       = fertilization$Crossed,
#    parallel = TRUE,
#    prior_delta      = prior("cauchy", list(0, 1/sqrt(2)), list(0, Inf)),
#    prior_rho        = prior("beta",   list(3, 3)),
#    prior_nu         = prior("exp",    list(1)),
#    prior_delta_null = prior("normal", list(0, 0.15), list(0, Inf)), #prior distribution and truncation
#    seed             = 3
#  )
#  
#  saveRDS(fit_01, file = "../models/Introduction/fit_01.RDS")
#  saveRDS(fit_10, file = "../models/Introduction/fit_10.RDS")
#  saveRDS(fit_1,  file = "../models/Introduction/fit_1.RDS")
#  saveRDS(fit_2,  file = "../models/Introduction/fit_2.RDS")
#  saveRDS(fit_3,  file = "../models/Introduction/fit_3.RDS")

## ----include = FALSE----------------------------------------------------------
# pre-load the fitted models to save time on compilation
library(RoBTT)

fit_01 <- readRDS(file = "../models/Introduction/fit_01.RDS")
fit_10 <- readRDS(file = "../models/Introduction/fit_10.RDS")
fit_1  <- readRDS(file = "../models/Introduction/fit_1.RDS")
fit_2  <- readRDS(file = "../models/Introduction/fit_2.RDS")
fit_3  <- readRDS(file = "../models/Introduction/fit_3.RDS")

## -----------------------------------------------------------------------------
library(RoBTT)
data("fertilization", package = "RoBTT")
head(fertilization)

## -----------------------------------------------------------------------------
# get overall summary
summary(fit_01)

## -----------------------------------------------------------------------------
# get individual model summaries
summary(fit_01, type = "models")

## -----------------------------------------------------------------------------
summary(fit_01, conditional = TRUE)

## -----------------------------------------------------------------------------
summary(fit_10)

## -----------------------------------------------------------------------------
summary(fit_10, type = "models")

## -----------------------------------------------------------------------------
summary(fit_1)

## -----------------------------------------------------------------------------
summary(fit_1, type = "models")

## -----------------------------------------------------------------------------
summary(fit_2, type = "models")

## -----------------------------------------------------------------------------
summary(fit_3, type = "models")

