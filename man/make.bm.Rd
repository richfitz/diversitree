\name{make.bm}
\alias{make.bm}
\alias{make.ou}
\alias{make.eb}
\alias{make.lambda}

\title{Brownian Motion and Related Models of Character Evolution}

\description{Create a likelihood function for models of simple Brownian
  Motion (BM), Ornstein-Uhlenbeck (OU), or Early Burst (EB) character
  evolution, or BM on a \dQuote{lambda} rescaled tree.  This function
  creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.bm(tree, states, states.sd=0, control=list())
make.ou(tree, states, states.sd=0, with.optimum=FALSE, control=list())
make.eb(tree, states, states.sd=0, control=list())
make.lambda(tree, states, states.sd=0, control=list())
}

\arguments{
  \item{tree}{A bifurcating phylogenetic tree, in \code{ape}
    \dQuote{phylo} format.}
  
  \item{states}{A vector of continuous valued character states.  This
    vector must be named with the tip labels of \code{tree}.}
  
  \item{states.sd}{An optional vector of measurement errors, as standard
    deviation around the mean.  If a single value is given it is used
    for all tips, otherwise the vector must be named as for
    \code{states}.}

  \item{with.optimum}{Should the optimum (often "theta") be considered
    a free parameter?  The default, \code{FALSE}, makes this behave like
    geiger's \code{fitContinuous}.  Setting \code{TRUE} leaves the
    optimim to be a free parameter to be estimated.  This setting can
    (currently) only be set to \code{TRUE} with
    \code{method="pruning"}.}

  \item{control}{A list of control parameters.  See details below.}
}

\details{The \code{control} argument is a named list of options.

  The main option is \code{method}.  Specifying
  \code{control=list(method="vcv")} uses a variance-covariance matrix
  based approach to compute the likelihood.  This is similar to the
  approach used by geiger, and is the default.

  Two alternative approaches are available.
  \code{control=list(method="pruning")} uses the transition density
  function for brownian motion along each branch, similar to how most
  methods in diversitree are computed.  This second approach is much
  faster for very large trees.  \code{control=list(method="contrasts")}
  uses Freckleton (2012)'s contrasts based approach, which is also much
  faster on large trees.

  When \code{method="pruning"} is specified, \code{backend="R"} or
  \code{backend="C"} may also be provided, which switch between a slow
  (and stable) R calculator and a fast (but less extensively tested) C
  calculator.  \code{backend="R"} is currently the default.

  The VCV-based functions are heavily based on \code{fitContinuous} in
  the \code{geiger} package.

  For non-ultrametric trees with OU models, computed likelihoods may
  differ because of the different root treatments.  This is particularly
  the case for models where the optimum is estimated.

  For the EB model, the parameter intepretation follows geiger; the 'a'
  parameter is equivalent to -log(g) in Bloomberg et al. 2003; when
  negative it indicates a decelerating rate of trait evolution over
  time.  When zero, it reduces to Brownian motion.

  % TODO: Describe root treatments!
}

\seealso{
  See
  \url{http://www.zoology.ubc.ca/prog/diversitree/examples/ou-nonultrametric/}
  for a discussion about calculations on non-ultrametric trees.
}

% TODO: Add example section and show equivalence to geiger?
\examples{
## Random data (following APE)
data(bird.orders)
set.seed(1)
x <- structure(rnorm(length(bird.orders$tip.label)),
               names=bird.orders$tip.label)

## With the VCV approach
fit1 <- find.mle(make.bm(bird.orders, x), .1)

## With the pruning calculations
lik.pruning <- make.bm(bird.orders, x, control=list(method="pruning"))
fit2 <- find.mle(lik.pruning, .1)

## All the same (need to drop the function from this though)
all.equal(fit1[-7], fit2[-7])

## If this is the same as the estimates from Geiger, to within the
## tolerances expected for the calculation and optimisation:
\dontrun{
fit3 <- fitContinuous(bird.orders, x)
all.equal(fit3$Trait1$lnl, fit1$lnLik)
all.equal(fit3$Trait1$beta, fit1$par, check.attributes=FALSE)
}

}

\author{Richard G. FitzJohn}
\keyword{models}
