\name{make.prior}
\alias{make.prior}
\alias{make.prior.exponential}
\alias{make.prior.uniform}

\title{Simple Prior Functions}

\description{Functions for generating prior functions for use with
  \code{\link{mcmc}}, etc.}

\usage{
make.prior.exponential(r)
make.prior.uniform(lower, upper)
}

\arguments{
  \item{r}{Scalar or vector of rate parameters}
  \item{lower}{Lower bound of the parameter}
  \item{upper}{Upper bound of the parameter}
}

\details{The exponential prior probability distribution has probability
  density \deqn{\sum_i r_i e^{-r_i x_i}}{sum(r[i]*exp(-r[i]*x[i]))}
  where the \eqn{i} denotes the \eqn{i}th parameter.  If \code{r} is a
  scalar, then the same rate is used for all parameters.

  These functions each return a function that may be used as the
  \code{prior} argument to \code{mcmc()}.  }

\keyword{models}
\author{Richard G. FitzJohn}
