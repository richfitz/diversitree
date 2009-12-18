\name{mcmc}
\alias{mcmc}
\title{Simple Markov Chain Monte Carlo with Slice Sampling}

\description{Run a simple-minded MCMC using slice samples for
  independent updating of each variable.}

\usage{
mcmc(f, x0, nsteps, w, lower, upper, type=TYPE.LOG, fail.value=-Inf,
     print.every=1, ...)
}

\arguments{
  \item{f}{Function to run MCMC on.}
  \item{x0}{Initial x location (vector).}
  \item{nsteps}{Number of MCMC steps to take.}
  \item{w}{(Approximate) width of the high probability region; this is
    the step size.  Must be same length as \code{x0}.}
  \item{lower}{Lower bounds on parameter space (vector).}
  \item{upper}{Upper bounds on parameter space (vector).}
  \item{type}{Type of value that the function \code{f} returns.  The
    default is to assume it returns a log-value, \code{bisse:::TYPE.RAW}
    assumes a raw probability and \code{bisse:::TYPE.NLOG} assumes a
    negative log-probability.}
  \item{fail.value}{Value to use where function evaluation fails.  The
    default (negative infinity) corresponds to zero probability.  Most
    values that fail are invalid for a given model (negative rates, etc)
    or have negligble probability, so this is reasonable.  Set to
    \code{NULL} to turn off checking.}
  \item{print.every}{The position and its probability will be printed
    every \code{print.every} generations.  Set this to 0 to prevent
    printing.}
  \item{...}{Arguments passed to the function \code{f}}
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

f <- make.bisse(phy, phy$tip.state)
names(pars) <- argnames(f)
f(pars)

## This produces about a sample a second, so takes a while.  The "upper"
## limit is a hard upper limit, above which the sampler will never let
## the parameter going (in effect, putting a uniform prior on the range
## lower..upper, and returning the joint distribution conditional on the
## parameters being in this range).
out <- mcmc(f, pars, nsteps=100, w=c(.1, .1, .03, .03, .01, .01),
            lower=rep(0, 6), upper=rep(2, 6), fail.value=-Inf)


## The argument 'w' works best when it is about the width of the "high
## probability" region for that parameter.  This takes the with of the
## 90% quantile range.  The resulting widths are about 10% slightly
## faster than the first guess.
w <- diff(sapply(out[2:7], quantile, c(.05, .95)))
out <- mcmc(f, pars, nsteps=100, w=w,
            lower=rep(0, 6), upper=rep(2, 6), fail.value=-Inf)

## You can do several things with this.  Look for ML-type values:
out[which.max(out$p),] # stochastic MLE

## Or look at the marginal distribution of parameters
hist(out$lambda0)

## Or look at the joint marginal distribution of pairs of parameters
plot(lambda0 ~ mu0, out)
}

\keyword{models}
