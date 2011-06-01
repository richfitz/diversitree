\name{asr}
\alias{asr}
\alias{asr.marginal}
\alias{asr.joint}
\alias{asr.stoch}

\alias{make.asr.marginal}
\alias{make.asr.joint}
\alias{make.asr.stoch}

\title{Ancestral State Reconstruction}

\description{Perform ancestral state reconstruction.  These functions
  are all generic and will dispatch on the class of the given likelihood
  functions.  Currently methods exist for all generics for Mk2, and
  marginal ancestral state reconstructions are supported for BiSSE.}

\usage{
asr.marginal(lik, pars, nodes=NULL, ...)
asr.joint(lik, pars, n=1, ...)
asr.stoch(lik, pars, n=1, ...) 

make.asr.marginal(lik, ...)
make.asr.joint(lik, ...)
make.asr.stoch(lik, ...)
}

\arguments{
  \item{lik}{A likelihood function.}
  
  \item{pars}{A vector of parameters, suitable for \code{lik}.}

  \item{nodes}{For \code{asr.marginal} only; an optional vector of
    nodes to return ancestral states for (using ape's index).  By
    default, all nodes are returned.}
  
  \item{n}{The number of samples to draw from the joint distribution, or
    number of stochastic reconstructions to make.}
  
  \item{...}{Additional arguments passed through to future methods}
}


\details{
  These three functions compute marginal, joint, and stochastic
  ancestral reconstructions.  The \code{make} versions return functions
  that can efficiently be used many times over.
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
