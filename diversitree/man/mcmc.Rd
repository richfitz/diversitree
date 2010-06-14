\name{mcmc}
\alias{mcmc}
\title{Simple Markov Chain Monte Carlo with Slice Sampling}

\description{Run a simple-minded MCMC using slice samples (Neal 2003)
  for independent updating of each variable.}

\usage{
mcmc(f, x0, nsteps, w, lower=-Inf, upper=Inf, fail.value=-Inf,
     prior=NULL, print.every=1, ...)
}

\details{The parameter \code{w} affects how many function evaluations
are required between sample updates, but it does not affect how fast the
MCMC ``mixes'' (Neal 2003).  In particular, \code{w} is not analagous to
the step sizes used in conventional Metropolis-Hastings updaters that
use some fixed kernel for updates.  Ideally, \code{w} would be set to
approximately the width of the high probability region.  I find that
chosing the distance between the 5\% and 95\% quantiles of the marginal
distributions of each parameter works well, computed from this
preliminary set of samples (see Examples).

If the MCMC is stopped by an interrupt (Escape on GUI versions of R,
Control-C on command-line version), it will return a truncated chain
with as many points as completed so far.
}

\arguments{
  \item{f}{Function to run MCMC on.}
  \item{x0}{Initial x location (vector).}
  \item{nsteps}{Number of MCMC steps to take.}
  \item{w}{(Approximate) width of the high probability region; this is
    the step size.  Must be same length as \code{x0}.  See Details for
  interpretation.}
  \item{lower}{Lower bounds on parameter space (scalar or vector of same
    length as \code{x0}).}
  \item{upper}{Upper bounds on parameter space (scalar or vector of same
    length as \code{x0}).}
  \item{fail.value}{Value to use where function evaluation fails.  The
    default (negative infinity) corresponds to zero probability.  Most
    values that fail are invalid for a given model (negative rates, etc)
    or have negligble probability, so this is reasonable.  Set to
    \code{NULL} to turn off checking.}
  \item{print.every}{The position and its probability will be printed
    every \code{print.every} generations.  Set this to 0 to prevent
    printing.}
  \item{prior}{An optional prior probability distribution function.
    This will be evaluated at each point that the sampler visits.  See
  \code{\link{make.prior}} for more information.}
  \item{...}{Arguments passed to the function \code{f}}
}

\examples{
## This example uses a 203 species phylogeny evolved under the BiSSE
## model.  This takes a substantial amount of time, so is not computed.

\dontrun{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

## First, create a likelihood function:
lik <- make.bisse(phy, phy$tip.state)
lik(pars)

## This produces about a sample a second, so takes a while.  The "upper"
## limit is a hard upper limit, above which the sampler will never let
## the parameter going (in effect, putting a uniform prior on the range
## lower..upper, and returning the joint distribution conditional on the
## parameters being in this range).
tmp <- mcmc(lik, pars, nsteps=100, w=rep(.1, 6), lower=0)

## The argument 'w' works best when it is about the width of the "high
## probability" region for that parameter.  This takes the with of the
## 90% quantile range.  The resulting widths are only slightly faster
## than the first guess.  Samples are generated about 1/s; allow 15
## minutes to generate 1000 samples.
w <- diff(sapply(tmp[2:7], quantile, c(.05, .95)))
out <- mcmc(lik, pars, nsteps=1000, w=w, lower=0)

## You can do several things with this.  Look for candidate ML points:
out[which.max(out$p),]

## Or look at the marginal distribution of parameters
profiles.plot(out["lambda0"], col.line="red")

## Or look at the joint marginal distribution of pairs of parameters
plot(lambda0 ~ mu0, out)
}
}

\seealso{
  \code{\link{make.bd}}, \code{\link{make.bisse}}, and
  \code{\link{make.mkn}}, all of which provide likelihood functions that
  are suitable for use with this function.  The help page for
  \code{\link{make.bd}} has further examples of using MCMC, and
  \code{\link{make.bisse}} has examples of using priors with MCMC.
}

\references{
  Neal R.M. 2003. Slice sampling. Annals of Statistics 31:705-767.
}

\keyword{models}
\author{Richard G. FitzJohn}