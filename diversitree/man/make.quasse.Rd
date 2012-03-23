\name{make.quasse}
\alias{make.quasse}
\alias{starting.point.quasse}

\title{Quantitative State Speciation and Extinction Model}

\description{Prepare to run QuaSSE (Quantitative State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.quasse(tree, states, states.sd, lambda, mu, control,
            sampling.f=NULL)
starting.point.quasse(tree, states, states.sd=NULL)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be a
    numeric real values.  Missing values (\code{NA}) are not yet handled.
    This vector must have names that correspond to the tip labels in the
    phylogenetic tree  (\code{tree$tip.label}).}

  \item{states.sd}{A scalar or vector corresponding to the standard
    error around the mean in states (the initial probability
    distribution is assumed to be normal).}

  \item{lambda}{A function to use as the speciation function.  The first
    argument of this must be \code{x} (see Details).}

  \item{mu}{A function to use as the extinction function.  The first
    argument of this must be \code{x} (see Details.)}

  \item{control}{A list of parameters for tuning the performance of the
    integrator.  A guess at reasonble values will be made here.  See
    Details for possible entries.}

  \item{sampling.f}{Scalar with the estimated proportion of extant
    species that are included in the phylogeny.  A value of \code{0.75}
    means that three quarters of extant species are included in the
    phylogeny.  By default all species are assumed to be known.}
}

\details{
  The \code{control} list may contain the following elements:

  \itemize{
    \item \code{method}: one of \code{fftC} or \code{fftR} to switch
    between \code{C} (fast) and \R (slow) backends for the integration.
    Both use non-adaptive fft-based convolutions.  Eventually, an
    adaptive methods-of-lines approach will be available.

    \item \code{dt.max}: Maximum time step to use for the integration.
    By default, this will be set to 1/1000 of the tree depth.  Smaller
    values will slow down calculations, but improve accuracy.

    \item \code{nx}: The number of bins into which the character space
    is divided (default=1024).  Larger values will be slower and more
    accurate.  For the \code{fftC} integration method, this should be an
    integer power of 2 (512, 2048, etc).

    \item \code{r}: Scaling factor that multiplies \code{nx} for a "high
    resolution" section at the tips of the tree (default=4, giving a
    high resolution character space divided into 4096 bins).  This helps
    improve accuracy while possibly tight initial probability
    distributions flatten out as time progresses towards the root.
    Larger values will be slower and more accurate.  For the \code{fftC}
    integration method, this should be a power of 2 (2, 4, 8, so that
    \code{nx*r} is a power of 2).
    
    \item \code{tc}: where in the tree to switch to the low-resolution
    integration (zero corresponds to the present, larger numbers moving
    towards the root).  By default, this happens at 10\% of the tree
    depth.  Smaller values will be faster, but less accurate.

    \item \code{xmid}: Mid point to center the character space.  By
    default this is at the mid point of the extremes of the character
    states.

    \item \code{tips.combined}: Get a modest speed-up by simultaneously
    integrating all tips?  By default, this is \code{FALSE}, but
    speedups of up to 25\% are possible with this set to \code{TRUE}.

    \item \code{w}: Number of standard deviations of the normal
    distribution induced by Brownian motion to use when doing the
    convolutions (default=5).  Probably best to leave this one alone.
  }
}

\section{Warning}{
  In an attempt at being computationally efficient, a substantial amount
  of information is cached in memory so that it does not have to be
  created each time.  However, this can interact poorly with the
  \code{multicore} package.  In particular, likelihood functions should
  not be made within a call to \code{mclapply}, or they will not share
  memory with the main \R thread, and will not work (this will cause an
  error, but should no longer crash \R).

  The method has less general testing than BiSSE, and is a little more
  fragile.  In particular, because of the way that I chose to implement
  the integrator, there is a very real chance of likelihood calculation
  failure when your data are a poor fit to the model; this can be
  annoyingly difficult to diagnose (you will just get a \code{-Inf} log
  likelihood, but the problem is often just caused by two sister species
  on short branches with quite different states).  There are also a
  large number of options for fine tuning the integration, but these
  aren't really discussed in any great detail anywhere.
}

\examples{
## Example showing simple integration with two different backends,
## plus the splits.
lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)

set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
                   single.lineage=FALSE, verbose=TRUE)

nodes <- c("nd13", "nd9", "nd5")
split.t <- Inf

pars <- c(.1, .2, 0, 2.5, .03, 0, .01)
pars4 <- unlist(rep(list(pars), 4))

sd <- 1/200
control.C.1 <- list(dt.max=1/200)
% control.C.2 <- c(control.C.1, tips.combined=TRUE)

control.R.1 <- list(dt.max=1/200, method="fftR")
lik.C.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars)) # -62.06409

% lik.C.2 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x, control.C.2)
% (ll.C.2 <- lik.C.2(pars)) # -62.06409

\dontrun{
## slow...
lik.R.1 <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.R.1 <- lik.R.1(pars)) # -62.06409
}

lik.s.C.1 <- make.quasse.split(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                               nodes, split.t, control.C.1)
(ll.s.C.1 <- lik.s.C.1(pars4)) # -62.06409
}

\author{Richard G. FitzJohn}
\keyword{models}
