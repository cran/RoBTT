## ----setup, include=FALSE-----------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/Truncated/fit1_trunc.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev      = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## ----include = FALSE, eval = FALSE--------------------------------------------
#  # pre-fit all models (easier to update the code on package update)
#  library(RoBTT)
#  
#  set.seed(42)
#  x1 <- rnorm(100, 0, 1)
#  x2 <- rnorm(100, 0, 1)
#  
#  stats1 <- boxplot.stats(x1)
#  lower_whisker1 <- stats1$stats[1]
#  upper_whisker1 <- stats1$stats[5]
#  
#  
#  stats2 <- boxplot.stats(x2)
#  lower_whisker2 <- stats2$stats[1]
#  upper_whisker2 <- stats2$stats[5]
#  
#  # Exclude outliers based on identified whiskers
#  x1_filtered <- x1[x1 >= lower_whisker1 & x1 <= upper_whisker1]
#  x2_filtered <- x2[x2 >= lower_whisker2 & x2 <= upper_whisker2]
#  
#  # Define whiskers for truncated likelihood application
#  whisker1 <- c(lower_whisker1, upper_whisker1)
#  whisker2 <- c(lower_whisker2, upper_whisker2)
#  
#  fit1_trunc <- RoBTT(
#    x1 = x1_filtered, x2 = x2_filtered,
#    truncation = list(x1 = whisker1, x2 = whisker2),
#    seed = 1, parallel = FALSE)
#  
#  
#  cut_off <- c(-2,2)
#  
#  x1 <- x1[x1 >= -2 & x1 <= 2]
#  x2 <- x2[x2 >= -2 & x2 <= 2]
#  
#  fit2_trunc  <- RoBTT(
#    x1 = x1, x2 = x2,
#    truncation = list(x = cut_off),
#    seed = 1, parallel = FALSE)
#  
#  saveRDS(fit1_trunc, file = "../models/Truncated/fit1_trunc.RDS")
#  saveRDS(fit2_trunc, file = "../models/Truncated/fit2_trunc.RDS")

## ----include = FALSE----------------------------------------------------------
# pre-load the fitted models to save time on compilation
library(RoBTT)

fit1_trunc <- readRDS(file = "../models/Truncated/fit1_trunc.RDS")
fit2_trunc <- readRDS(file = "../models/Truncated/fit2_trunc.RDS")

## ----echo=T,error=F, message=F, warning=F, results="hide"---------------------
# Install RoBTT from CRAN
# install.packages("RoBTT")

# Load the RoBTT package
library(RoBTT)

## ----echo=T,error=F, message=F, warning=F, results="hide"---------------------
set.seed(42)
x1 <- rnorm(100, 0, 1)
x2 <- rnorm(100, 0, 1)

## ----echo=T,error=F, message=F, warning=F, results="hide"---------------------
# Identify outliers using boxplot statistics for each group
stats1 <- boxplot.stats(x1)
lower_whisker1 <- stats1$stats[1]
upper_whisker1 <- stats1$stats[5]


stats2 <- boxplot.stats(x2)
lower_whisker2 <- stats2$stats[1]
upper_whisker2 <- stats2$stats[5]

# Exclude outliers based on identified whiskers
x1_filtered <- x1[x1 >= lower_whisker1 & x1 <= upper_whisker1]
x2_filtered <- x2[x2 >= lower_whisker2 & x2 <= upper_whisker2]

# Define whiskers for truncated likelihood application
whisker1 <- c(lower_whisker1, upper_whisker1)
whisker2 <- c(lower_whisker2, upper_whisker2)

## -----------------------------------------------------------------------------
summary(fit1_trunc, group_estimates = TRUE)

## -----------------------------------------------------------------------------
summary(fit1_trunc, group_estimates = TRUE, type = "models")

## ----echo=T,error=F, message=F, warning=F, results="hide"---------------------
cut_off <- c(-2,2)

x1 <- x1[x1 >= -2 & x1 <= 2]
x2 <- x2[x2 >= -2 & x2 <= 2]

