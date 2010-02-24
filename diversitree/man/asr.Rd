\name{asr}
\alias{asr}
\alias{asr.marginal}
\alias{asr.joint}
\alias{asr.stoch}

\title{Ancestral State Reconstruction}

\description{Perform ancestral state reconstruction.  These functions
  are all generic and will dispatch on the class of the given likelihood
  functions.  Currently methods exist for all generics for Mk2, and
  marginal ancestral state reconstructions are supported for BiSSE.}

\usage{
asr.marginal(lik, pars, nodes=NULL, ...)
asr.joint(lik, pars, n=1, simplify=TRUE, ...)
asr.stoch(lik, pars, n=1, ...) 
}

\arguments{
  \item{lik}{A likelihood function.}
  \item{pars}{A vector of parameters, suitable for \code{lik}.}
  \item{nodes}{For \code{asr.marginal} only; an optional vector of
    nodes to return ancestral states for (using ape's index).  By
    default, all nodes are returned.}
  \item{n}{The number of samples to draw from the joint distribution, or
    number of stochastic reconstructions to make.}
  \item{simplify}{For \code{asr.joint}, return a vector of states when
    \code{n=1}? Otherwise leaves as a one row matrix.}
  \item{...}{Additional arguments passed through to future methods}
}


\details{
  These three functions compute marginal, joint, and stochastic
  ancestral reconstructions.
}

\value{
  The return values of the functions are likely to change in the near
  future.  Watch out!
}

\seealso{
  \link{asr.mkn} and \link{asr.bisse} for methods specific to particular
  classes, with examples of use.
}

\author{Richard G. FitzJohn}
\keyword{model}
