\name{find.mle}
\alias{find.mle}
\alias{coef.fit.mle}
\alias{logLik.fit.mle}
\alias{anova.fit.mle}

% \alias{find.mle.bisse}
% \alias{starting.point.bisse}
% \alias{starting.point} % DEP
% \alias{bisse.starting.point} % DEP

\title{Maximimum Likelihood Inference}

\description{Find the maximum likelihood point of a model by nonlinear
  optimisation.  \code{find.mle} is generic, and allows different
  default behaviour for different likelihood functions.}

\usage{
find.mle(func, x.init, method, ...)
\method{coef}{fit.mle}(object, ...)
\method{logLik}{fit.mle}(object, ...)
\method{anova}{fit.mle}(object, ...)
%\method{find.mle}{bisse}(func, x.init, method, fail.value=NA, ...)
%bisse.starting.point(tree, q.div=5)
%starting.point(tree, q.div=5)
}

\arguments{
  \item{func}{A likelihood function.  This is assumed to return the log
    likelihood (see Details).  The function must take a vector of
    parameters as the first argument.}
  \item{x.init}{Initial starting point for the optimisation.}
  \item{method}{Method to use for optimisation.  May be one of "optim",
    "subplex", "nlminb", "nlm" (partial unambigious string is allowed).}
  \item{...}{For \code{find.mle}, additional arguments passed through to
    the methods, optimisation routines, or to the likelihood function
    \code{func} - see Details.  For \code{anova}, this is one or more
    models to compare against the model \code{object} (either submodels
    or supermodels or the test is meaningless).}
  \item{object}{A fitted model, returned by \code{find.mle}.}
}

\details{
  \code{find.mle} starts a search for the maximum likelihood (ML)
  parameters from a starting point \code{x.init}.  \code{x.init} should
  be the correct length for \code{func}, so that \code{func(x.init)}
  returns a valid likelihood.  However, if \code{func} is a constrained
  function (via \code{\link{constrain}}) and \code{x.init} is the
  correct length for the unconstrained function then an attempt will be
  made to guess a valid starting point.  This will often do poorly and a
  warning will be given.

  Different methods will be dispatched for different types of likelihood
  functions.  Currently all models in \code{diversitree} are supported
  (\code{bisse}, \code{mk2}, \code{mkn}, \code{bd}, and \code{yule}).
  With the exception of the Yule pure-birth process, these methods just
  specify different default arguments for the underlying optimisation
  routines (the Yule model has an analytical solution, and no
  optimisation step is required).  Generally, it will not be necessary
  to specify the \code{method} argument to \code{find.mle} as a sensible
  method is chosen during dispatch.

  The \code{...} argument may contain additional arguments for the
  function \code{func}.  This includes things like \code{condition.surv}
  for conditioning on survival in BiSSE, birth-death, and Yule models.
  Specify this as
  \preformatted{
    find.mle(lik, x.init, condition.surv=TRUE)
  }
  (see the Examples).
  
  Different \code{method} arguments take different arguments passed
  through \code{...} to control their behaviour:

  \code{method="optim"}: Uses \R's \code{\link{optim}} function for the
  optimisation.  This allows access to a variety of general purpose
  optimisation algorithms.  The method \emph{within} \code{optim} can be
  chosen via the argument \code{optim.method}, which is set to
  "L-BFGS-B" by default (box constrained quasi-Newton optimisation).
  This should be suitable for most uses.  See the \code{method} argument
  of \code{\link{optim}} for other possibilities.  If \code{"L-BFGS-B"}
  is used, then upper and lower bounds may be specified by the arguments
  \code{lower} and \code{upper}.  The argument \code{control} can be
  used to specify other control parameters for the algorithms - see
  \code{\link{optim}} for details.  Most of the \code{optim} algorithms
  require finite values be returned at every evaluated point.  This is
  often not possible (extreme values of parameters or particular
  combinations may have zero likelihood and therefore -Inf
  log-likelihood).  To get around this, the argument \code{fail.value}
  can be used to specify a fallback value.  By default this is set to
  \code{func(x.init) - 1000}, which should work reasonably well for most
  cases.

  \code{method="subplex"}: Uses the "subplex" algorithm (a variant of
  the downhill simplex/Nelder-Mead algorithm that uses Nelder-Mead on a
  sequence of subspaces).  This algorithm generally requires more
  evaluations than \code{optim}-based optimisation, but does not require
  approximation of derivatives and seems to find the global optimum more
  reliably (though often less precisely).  Additional arguments are
  \code{control} to control aspects of the search (see
  \code{\link{subplex}} for details).  The argument \code{fail.value}
  can be used as in \code{method="optim"}, but by default \code{-Inf}
  will be used on failure to evaluate, which is generally appropriate.

  \code{method="nlminb"}: Uses the function \code{\link{nlminb}} for
  optimisation, so that optimising a Mk2/Mkn likelihood function behaves
  as similarly as possible to \code{ape}'s \code{\link{ace}} function.
  As for \code{method="optim"}, lower and upper bounds on parameters may
  be specified via \code{lower} and \code{upper}.  \code{fail.value} can
  be used to control behaviour on evaluation failure, but like
  \code{method="subplex"}, \code{-Inf} is used which should work in most
  cases.  Additional control parameters may be passed via \code{control}
  - see \code{link{nlminb} for details}.  This function is not generally
  recommended for use.

  \code{method="nlm"}: Uses the function \code{\link{nlm}} for
  optimisation, so that optimising a birth-death likelihood function
  behaves as similarly as possible to \code{ape}'s
  \code{\link{birthdeath}} function.  Takes the same additional
  arguments as \code{method="nlminb"} (except that \code{fail.value}
  behaves as for \code{method="optim"}).  Like \code{method="nlminb"},
  this is not recommended for general use.

  \code{code} and \code{logLik} methods exist for \code{fit.mle} objects
  so that parameters and log-likelihoods may be extracted.  This also
  allows use with \code{\link{AIC}}.

  Simple model comparison by way of likelihood ratio tests can be
  performed with \code{\link{anova}}.  See Examples for usage.
}

\value{
  A list of class \code{fit.mle}, with components
  \item{par}{The estimated parameters}
}


\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

## Here is the 203 species tree with the true character history coded.
## Red is state '1', which has twice the speciation rate of black (state
## '0').
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, cex=.5, show.node.state=FALSE)

## Make a BiSSE likelihood function
lik <- make.bisse(phy, phy$tip.state)
lik(pars)

## This takes ~30s to run, so is not enabled by default
\dontrun{
## Fit the full six-parameter model
fit <- find.mle(lik, pars)
fit[1:2]

coef(fit)   # Named vector of six parameters
logLik(fit) # -659.92
AIC(fit)    # 1331.85

## Subplex can sometimes find better solutions, though not this time:
fit.s <- find.mle(lik, pars, method="subplex")
logLik(fit.s) - logLik(fit)

## find.mle works with constrained models (see \link{constrain}).  Here
## the two speciation rates are constrained to be the same as each
## other.  Optimisation here is done with subplex because it *does*
## return a better point than the default.
lik.l <- constrain(lik, lambda0 ~ lambda1)
fit.l <- find.mle(lik.l, pars[-2], method="subplex")
logLik(fit.l) # 663.40

## Compare the models with \link{anova} - this shows that the more
## complicated model with two separate speciation rates fits
## significantly better than the simpler model with equal rates
## (p=0.008).
anova(fit, equal.lambda=fit.l)
}
}

\keyword{models}
\author{Richard G. FitzJohn}