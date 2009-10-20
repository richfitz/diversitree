\name{find.mle}
\alias{find.mle}
\alias{find.mle.bisse}
\alias{logLik.mle.bisse}
\alias{starting.point}
\alias{anova.mle.bisse}

\title{Maximimum Likelihood Inference}

\description{Find the maximum likelihood point of a model by nonlinear
  optimisation.}

\usage{
find.mle(func, x.init, ...)
\method{find.mle}{bisse}(func, x.init,
               method=c("L-BFGS-B", "Nelder-Mead", "subplex"),
               control=list(), lower=NULL, upper=NULL,
               fail.value=NULL, hessian=FALSE, ...)
starting.point(tree, q.div=5)
\method{logLik}{mle.bisse}(object, ...)
\method{anova}{mle.bisse}(object, ...)
}

\arguments{
  \item{func}{A likelihood function.  This is assumed to return the log
    likelihood (see Details).  The function must take a vector of
    parameters as the first argument.}
  \item{x.init}{Initial starting point for the optimisation.}
  \item{method}{Algorithm used for the optimisation.  \code{L-BFGS-B}
    (the default) and \code{Nelder-Mead} use \code{\link{optim}}.
    Because \code{L-BFGS-B} uses finite-differencing to approximate the
    gradient, it is prone to getting stuck where the likelihood surface
    is fairly flat.  Both of these algorithms are strict hill-climbers.
    The method \code{subplex} is much slower (possibly 10x the number of
    evaluations) but tends explores more of parameter space.}
  \item{control}{List of options to \code{\link{optim}} or
    \code{subplex}}
  \item{lower}{Lower limits for the parameters.  If omitted, the
    \code{bisse} method adds sensible lower limits in the absence of
    additional information (see Details).}
  \item{upper}{Upper limits for the parameters. If omitted, no upper
    bound is assumed (this is usually fine).}
  \item{fail.value}{Value to use if the likelihood calculation fails.
    \code{\link{optim}} requires finite values for all evaluations, so if you
    have a problem with failed evaluations set this
    to be an unlikely value (but not too much smaller than observed
    numbers).  For \code{subplex}, it is better to leave this
    alone, as undefined likelihoods are treated as impossibly likely by
    \code{subplex} already.}
  \item{hessian}{Logical, indicating if the hessian should be computed;
    this uses the algorithm in the \code{numDeriv} package and can add a
    significant number of function evaluations to approximate the
    hessian (however, the answer tends to be fairly accurate).}
  \item{...}{For \code{find.mle}, arguments passed to the function
    \code{f}.  For \code{logLik} this is ignored.  For \code{anova},
    this is one or more models to compare against the model
    \code{object} (either submodels or supermodels or the test is
    meaningless).}
  \item{tree}{A phylogeny}
  \item{q.div}{Ratio of diversification rate to character change rate}
  \item{object}{ML object returned by \code{find.mle}}
}

\details{
  Arguments for constrained BiSSE models are processed a little
  differently.  Constrained BiSSE models have fewer than six parameters;
  if a constrained model is detected and \code{x.init} has six
  arguments, then the arguments \code{lower} and \code{x.init} are
  subset so that only the appropriate elements are used.  This does not
  work for other vectorised arguments to \code{\link{optim}} (e.g.,
  \code{upper}).

  Lower limits are required for optimisation to work well.  Being rates,
  all parameters have a natural lower limit of zero, but the likelihood
  will be undefined with polymorphic state data if both \code{q} values
  are zero, and undefined if both \code{lambda} values are zero.  The
  lower limit is set to be \code{1e-4} for both \code{lambda} and
  \code{q} values.  A warning will be given if the final point lies on
  any given constraint; it would be worth rerunning the analysis with
  contraints relaxed if this warning is given.

  \code{starting.point} produces a poor starting point to start from,
  based on the character-independent birth-death model.

  See \code{\link{find.mle}} for an example of \code{anova}.
}

\keyword{models}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

plot(phy)

f <- make.bisse(phy, phy$tip.state)
f(pars)

## Each ML search will take about 5-10s, so is not run by default.  The
## warnings are harmless, and just indicate points where the search
## tried parameters that do not produce valid likelihoods.
\dontrun{
fit.mle <- find.mle(f, pars)
fit.mle

## Subplex can sometimes find better solutions, though not this time:
fit.mle.subplex <- find.mle(f, pars, method="subplex")
fit.mle.subplex$lnLik - fit.mle$lnLik
}
}