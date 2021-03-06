% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mxPenaltyFunctions.R
\name{mxRegularizeLASSO}
\alias{mxRegularizeLASSO}
\title{MxRegularizeLASSO}
\usage{
mxRegularizeLASSO(
  what,
  name,
  lambda = 0,
  lambda.step = 0.01,
  lambda.max = NA,
  lambda.min = NA
)
}
\arguments{
\item{what}{A list of parameters to regularize}

\item{name}{Name of the regularizer object}

\item{lambda}{strength of the penalty to be applied at starting values (default 0)}

\item{lambda.step}{step function for lambda step (default .01)}

\item{lambda.max}{end of lambda range (default .4)}

\item{lambda.min}{minimum lambda value (default lambda)}
}
\description{
Least Absolute Selection and Shrinkage Operator regularization
}
