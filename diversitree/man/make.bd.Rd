\name{make.bd}
\alias{make.bd}
\alias{make.yule}

\title{Constant Rate Birth-Death Models}

\description{Prepare to run a constant rate birth-death model on a
  phylogenetic tree.  This fits the Nee et al. 1994 equation,
  duplicating the \code{birthdeath} function in ape.  Differences with
  that function include (1) the function is not constrained to positive
  diversification rates (mu can exceed lambda), (2) [eventual] support
  for both random taxon sampling and unresolved terminal clades (but see
  \code{bd.ext}), and (3) run both MCMC and MLE fits to birth death
  trees.}

\usage{
make.bd(tree, times=branching.times(tree),
        sampling.f=NULL, unresolved=NULL)
make.yule(tree, times=branching.times(tree),
          sampling.f=NULL, unresolved=NULL)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{times}{Vector of branching times, as returned by
    \code{branching.times}.  You don't need to use this unless you know
    that you need to use this.  Don't use it at the same time as
    \code{tree}.}
  \item{unresolved}{(Not yet used.)  Unresolved clade information.}
  \item{sampling.f}{(Not yet used.)  Probability of an extant species
    being included in the phylogeny (sampling fraction).  By default,
    all extant species are assumed to be included.}
}

\details{
  \code{make.bd} returns a function of class \code{bd}.
  This function has argument list (and default values)
  \preformatted{
    f(pars, prior=NULL, fail.value=NULL)
  }
  This will change in future versions.
    The arguments are interpreted as
    \itemize{
      \item \code{pars} A vector of two parameters, in the order \code{lambda},
      \code{mu}.
      \item \code{prior} A vector of rates for an exponential prior.
      The prior probability distribution has probability density
      \deqn{r_i e^{-r_i x_i}}{r_i*exp(-r_i*x_i)}
      where the \eqn{i} denotes the \eqn{i}th parameter.  There is no
      ability to specify a form of prior other than an exponential
      function here (yet), and this part of the inteface is subject to
      change!
      \item \code{fail.value} Optional value to use if the likelihood
      calculation fails.
    }
  }

\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bisse}} for state-dependent birth-death models.
}
\examples{
## Simulate a tree under a constant rates birth-death model and look at
## the marginal likelihood of the speciation/extinction parameters:
set.seed(1)
phy <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]
lik <- make.bd(phy)
fit <- find.mle(lik, c(.1, .03))

samples <- mcmc(lik, fit$par, nsteps=5000,
                lower=c(-Inf, -Inf), upper=c(Inf, Inf), w=c(.1, .1),
                fail.value=-Inf, print.every=200)
samples$r <- with(samples, lambda - mu)
samples$a <- with(samples, mu / lambda)

col <- c("red", "blue", "green3")
profiles.plot(samples[c("lambda", "mu", "r")], col.line=col, las=1,
             ylab="Probability density", xlab="Parameter estimate")
legend("topright", c("lambda", "mu", "r"), fill=col)
abline(v=0, lty=2)

## Special case methods are worked out for the Yule model, for which
## analytic solutions are available:
lik.yule <- make.yule(phy)
lik.mu0 <- constrain(lik, mu ~ 0)

## The same to a reasonable tolerance:
fit.yule.1 <- find.mle(lik.yule)
fit.mu0.1 <- find.mle(lik.mu0, .1)
all.equal(fit.yule.1[1:2], fit.mu0.1[1:2], tolerance=1e-6)

## Without conditioning on survival
fit.yule.2 <- find.mle(lik.yule, condition.surv=FALSE)
fit.mu0.2 <- find.mle(lik.mu0, .1, condition.surv=FALSE)
all.equal(fit.yule.2[1:2], fit.mu0.2[1:2], tolerance=1e-6)
}
  
 
\author{Richard G. FitzJohn}
\keyword{models}
