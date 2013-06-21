\name{make.bd.t}
\alias{make.bd.t}

\title{Time-varing Birth-Death Models}

\description{Create a likelihood function for the birth-death model,
  where birth and/or death rates are arbitrary functions of time.}

\usage{
make.bd.t(tree, functions, sampling.f=NULL, unresolved=NULL,
          control=list(), truncate=FALSE, spline.data=NULL)
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format.}

  \item{functions}{A named list of functions of time.  See details.}

  \item{sampling.f}{Probability of an extant species being included in
    the phylogeny (sampling fraction).  By default, all extant species are
    assumed to be included.}

  \item{unresolved}{Not yet included: present in the argument list for
    future compatibility with \code{\link{make.bd}}.}
  
  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}

  \item{truncate}{Logical, indicating if functions should be truncated
    to zero when negative (rather than causing an error).  May be
    scalar (applying to all functions) or a vector (of length 2).}

  \item{spline.data}{List of data for spline-based time functions.  See
    details}.
}

\examples{
## First, show equivalence to the plain Birth-death model.  This is not
## a very interesting use of the functions, but it serves as a useful
## check.

## Here is a simulated 25 species tree for testing.
set.seed(1)
pars <- c(.1, .03)
phy <- trees(pars, "bd", max.taxa=25)[[1]]

## Next, make three different likelihood functions: a "normal" one that
## uses the direct birth-death calculation, an "ode" based one (that
## uses numerical integration to compute the likelihood, and is
## therefore not exact), and one that is time-varying, but that the
## time-dependent functions are constant.t().
lik.direct <- make.bd(phy)
lik.ode <- make.bd(phy, control=list(method="ode"))
lik.t <- make.bd.t(phy, c("constant.t", "constant.t"))

lik.direct(pars) # -22.50267

## ODE-based likelihood calculations are correct to about 1e-6.
lik.direct(pars) - lik.ode(pars)

## The ODE calculation agrees exactly with the time-varying (but
## constant) calculation.
lik.ode(pars) - lik.t(pars)

## Next, make a real case, where speciation is a linear function of
## time.
lik.t2 <- make.bd.t(phy, c("linear.t", "constant.t"))

## Confirm that this agrees with the previous calculations when the
## slope is zero
pars2 <- c(pars[1], 0, pars[2])
lik.t2(pars2) - lik.t(pars)

## The time penalty comes from moving to the ODE-based solution, not
## from the time dependence.
system.time(lik.direct(pars)) # ~ 0.000
system.time(lik.ode(pars))    # ~ 0.003
system.time(lik.t(pars))      # ~ 0.003
system.time(lik.t2(pars2))    # ~ 0.003

fit <- find.mle(lik.direct, pars)
fit.t2 <- find.mle(lik.t2, pars2)

## No significant improvement in model fit:
anova(fit, time.varying=fit.t2)
}

\author{Richard G. FitzJohn}
\keyword{models}
