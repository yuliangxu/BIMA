% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{M_regression_GS}
\alias{M_regression_GS}
\title{Image on scalar regression with GP prior for alpha instead of STGP}
\usage{
M_regression_GS(
  data,
  init,
  region_idx_cpp,
  kernel,
  n_mcmc,
  display_progress = TRUE
)
}
\arguments{
\item{data}{A List that contains
\itemize{
\item M The image predictor, n by p
\item X The scalar exposure variable, n by 1
\item C The q confounders, n by q
}}

\item{init}{A list object that contains the following element
\itemize{
  \item zetam A vector of length q, initial value for zetam
  \item sigma_alpha A numeric scalar, intial value for sigma_alpha
  \item sigma_M A numeric scalar, intial value for sigma_M
  \item theta_eta A n by L matrix, initial value for theta_eta
  \item sigma_eta A numeric scalar, intial value for sigma_eta
  \item theta_alpha A vector of length L, intial value for theta_alpha
  \item a A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
  \item b A numeric variable for the Inverse-Gamma(a,b), priors for \eqn{\sigma^2_Y,\sigma^2_\beta}
  \item sigma_zetam A numeric scalar, intial value for sigma_zetam
}}

\item{region_idx_cpp}{A list object of length num_region, each element is a vector of
the indices of each voxel in that region. Note that this index starts from 0.}

\item{kernel}{A List object
\itemize{
  \item D A list of length num_region, each element contains the eigen-values in one region
  \item Q A list of length num_region, each element contains the basis function (p_r by L_r) in one region
}}

\item{n_mcmc}{An integer to indicate the total number of MCMC iterations}

\item{display_progress}{True for displaying progress bar}
}
\value{
A list of
\itemize{
  \item theta_eta
  \item zetam_mcmc
  \item sigma_M2_inv_mcmc
  \item sigma_alpha2_inv_mcmc
  \item sigma_zetam2_inv_mcmc
  \item sigma_eta2_inv_mcmc
  \item logLL_mcmc
  \item Timer
}
}
\description{
This is not a main function in BIMA package. It only provides an alternative to the STGP prior.
A basis decomposition is used. The main coefficient beta follows STGP prior.
Kernel matrices need to be prespeficified
}
