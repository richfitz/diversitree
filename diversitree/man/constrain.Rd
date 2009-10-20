\name{constrain}
\alias{constrain}
\title{Constrain Parameters of a Model}

\description{Constrain a model to make submodels with fewer parameters.
  If \code{f} is a function that takes a vector \code{x} as its first
  argument, this function returns a new function that takes a
  shorter vector \code{x} with some elements constrained in some way;
  parameters can be fixed to particular values, constrained to be the
  same as other parameters, or arbitrary expressions of free
  parameters.
}

\usage{
constrain(f, ..., names=argnames(f))
}
\arguments{
  \item{f}{A function to constrain.}
  \item{...}{Formulae indicating how the function should be constrained
    (see Details and Examples).}
  \item{names}{Character vector of names, the same length as the number
    of parameters in \code{x}.  Use this only if \code{\link{argnames}}
    does not return a vector for your function.}
}

\details{
  The relationships are specified in the form \code{target ~ rel}, where
  \code{target} is the name of a vector to be constrained, and
  \code{rel} is some relationship.  For example \code{lambda0 ~ lambda1}
  would have the effect of making the parameters \code{lambda0} and
  \code{lambda1} take the same value.

  The \code{rel} term can be a constant (e.g., \code{target ~ 0}),
  another parameter (as above) or some expression of the parameters
  (e.g., \code{lambda0 ~ 2 * lambda1} or
  \code{lambda0 ~ lambda1 - mu1}).

  Terms that appear on the right hand side of an expression may not be
  constrained in another expression, and no term may be constrained
  twice.
}

\section{Warning}{
  Only a few checks are done to ensure that the resulting function makes
  any sense; it is possible that I have missed some cases.  It is
  probably dangerous to constrain a constrained function (just start
  again from the raw likelihood function) and there is currently no way
  of modifying constrained functions to remove the constraints.  These
  weaknesses will be addressed in a future version.
}

\value{
  This function returns a constrained function that can be passed
  through to \code{\link{find.mle}} and \code{mcmc}.  It will behave
  like any other function.  However, it has a modified \code{class}
  attribute so that some methods will dispatch differently
  (\code{argnames}, for example).  All arguments in addition to \code{x}
  will be passed through to the original function \code{f}.

  For help in designing constrained models, the returned function has
  an additional argument \code{pars.only}, when this is \code{TRUE} the
  function will return a named vector of arguments rather than evaluate
  the function (see Examples).
}

\examples{
## Same example likelihood function as for \link{find.mle}:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)
f <- make.bisse(phy, phy$tip.state)

argnames(f) # Canonical argument names (set by default)

## Equal speciation rates
g <- constrain(f, lambda0 ~ lambda1)
f(c(.1, .1, pars[3:6])) == g(c(.1, pars[3:6]))

## More complicated:
h <- constrain(f, lambda0 ~ lambda1 / 2, mu0 ~ mu1, q01 ~ 0)
argnames(h) # lambda1, mu1, q10
f(c(.1, .2, .03, .03, 0, .01)) == h(c(.2, .03, .01))
h(c(.2, .03, .01), pars.only=TRUE)

\dontrun{
## find.mle knows how to do a search on these reduced models.  Do a
## likelihood search on both the full model and the equal-lambda model:
fit.full <- find.mle(f, pars)
logLik(fit.full) # -659.92
fit.lambda <- find.mle(g, pars)
logLik(fit.lambda) # -663.40

## Hypothesis testing:
anova(fit.full, equal.lambda=fit.lambda)
}
}

\keyword{programming}
