\name{combine}
\alias{combine}

\title{Combine Several Likelihood Functions Multiplicatively}

\description{
  Combine several likelihood functions, so that the new functions gives
  the product of all likelihoods (the sum of the log likelihoods).  This
  assumes that all likelihoods are independent from one another!

  This function is little tested.  Use at your own risk!
}

\usage{
combine(liks)
}

\arguments{
  \item{liks}{A list of likelihood functions.  All must be of  the same
    type, with the same argnames, and not constrained.}
}

\author{Richard G. FitzJohn}
\keyword{programming}
