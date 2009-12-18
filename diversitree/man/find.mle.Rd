\name{find.mle}
\alias{find.mle}

%% TODO: I need a demonstration of both of these in examples
\alias{logLik.fit.mle}
\alias{anova.fit.mle}

%% TODO: This all needs moving into its own file.
\alias{find.mle.bisse}
\alias{starting.point.bisse}
\alias{starting.point} % DEP
\alias{bisse.starting.point} % DEP

\title{Maximimum Likelihood Inference}

\description{Find the maximum likelihood point of a model by nonlinear
  optimisation.}

\usage{
find.mle(func, x.init, method, ...)
\method{find.mle}{bisse}(func, x.init, method, fail.value=NA, ...)
starting.point.bisse(tree, q.div=5)
bisse.starting.point(tree, q.div=5)
starting.point(tree, q.div=5)
\method{logLik}{fit.mle}(object, ...)
\method{anova}{fit.mle}(object, ...)
}

\arguments{
  \item{func}{A likelihood function.  This is assumed to return the log
    likelihood (see Details).  The function must take a vector of
    parameters as the first argument.}
  \item{x.init}{Initial starting point for the optimisation.}
  \item{method}{WARNING: this is out of date, and applies only to
    bisse.
    Algorithm used for the optimisation.  \code{L-BFGS-B}
    (the default) and \code{Nelder-Mead} use \code{\link{optim}}.
    Because \code{L-BFGS-B} uses finite-differencing to approximate the
    gradient, it is prone to getting stuck where the likelihood surface
    is fairly flat.  Both of these algorithms are strict hill-climbers.
    The method \code{subplex} is much slower (possibly 10x the number of
    evaluations) but tends explores more of parameter space.}
%   \item{control}{List of options to \code{\link{optim}} or
%     \code{subplex}}
  \item{fail.value}{Value to use if the likelihood calculation fails.
    \code{\link{optim}} requires finite values for all evaluations, so if you
    have a problem with failed evaluations set this
    to be an unlikely value (but not too much smaller than observed
    numbers).  For \code{subplex}, it is better to leave this
    alone, as undefined likelihoods are treated as impossibly likely by
    \code{subplex} already.}
%   \item{hessian}{Logical, indicating if the hessian should be computed;
%     this uses the algorithm in the \code{numDeriv} package and can add a
%     significant number of function evaluations to approximate the
%     hessian (however, the answer tends to be fairly accurate).}
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
  arguments, it is subseted so that only the appropriate elements are
  used.  This does not work for other vectorised arguments to
  \code{\link{optim}} (e.g., \code{upper}).

  \code{starting.point.bisse} produces a heuristic starting point to
  start from, based on the character-independent birth-death model.  You
  can probably do better than this; see the vignette, for example.
  \code{bisse.starting.point} is the same code, but deprecated in favour
  of \code{startng.point.bisse} - it will be removed in a future
  version.

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
fit.mle[1:2]

## Subplex can sometimes find better solutions, though not this time:
fit.mle.subplex <- find.mle(f, pars, method="subplex")
logLik(fit.mle.subplex) - logLik(fit.mle)
}
}