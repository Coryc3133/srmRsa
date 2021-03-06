% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{srm_rsa_params_tbl}
\alias{srm_rsa_params_tbl}
\title{RSA SRM Parameter Table}
\usage{
srm_rsa_params_tbl(fitted_srm_rsa, caption = "SRM RSA Parameters")
}
\arguments{
\item{fitted_srm_rsa}{The fitted SRM RSA object from which you want the
parameter table.}

\item{caption}{The caption for the table. Defaults to generic label.}
}
\value{
a wrapper for knitr::kable that extracts and tables
the parameters from the SRMRSA model.
}
\description{
This function extracts the SRM components, the slopes from the regression
used to define the response surface, the surface parameters (a1 through a5),
and the intercepts and variances of x and y variables from RSA.
}
\keyword{RSA}
