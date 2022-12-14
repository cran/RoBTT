% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check-input-and-settings.R
\name{RoBTT_control}
\alias{RoBTT_control}
\alias{set_control,}
\alias{set_convergence_checks}
\alias{set_control}
\title{Convergence checks of the fitting process}
\usage{
set_convergence_checks(max_Rhat = 1.05, min_ESS = 500)

set_control(adapt_delta = 0.8, max_treedepth = 15, bridge_max_iter = 1000)
}
\arguments{
\item{max_Rhat}{maximum value of the R-hat diagnostic.
Defaults to \code{1.05}.}

\item{min_ESS}{minimum estimated sample size.
Defaults to \code{500}.}

\item{adapt_delta}{tuning parameter of HMC.
Defaults to \code{0.80}.}

\item{max_treedepth}{tuning parameter of HMC.
Defaults to \code{15}.}

\item{bridge_max_iter}{maximum number of iterations for the
\link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}}
}
\value{
\code{set_control} returns a list of control settings
and \code{set_convergence_checks} returns a list of convergence checks settings.
}
\description{
Set values for the convergence checks of the fitting process.
}
