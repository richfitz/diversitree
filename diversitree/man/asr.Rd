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
asr.joint(lik, pars, n=1, root.state=NA, simplify=TRUE, ...)
asr.stoch(lik, pars, n=1, root.state=NA, ...)
}

\arguments{
  \item{lik}{A likelihood function.}
  \item{pars}{A vector of parameters, suitable for \code{lik}.}
  \item{nodes}{For \code{asr.marginal} only; an optional vector of
    nodes to return ancestral states for (using ape's index).  By
    default, all nodes are returned.}
  \item{n}{The number of samples to draw from the joint distribution, or
    number of stochastic reconstructions to make.}
  \item{root.state}{For \code{asr.marginal} and \code{asr.joint}; the
    state to start the root in.  If \code{NA} is given (the default),
    the root state is drawn from the conditional likelihoods at the
    root.}
  \item{simplify}{Should the states be simplified into a matrix or
    left unsimplified as a list?}
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
