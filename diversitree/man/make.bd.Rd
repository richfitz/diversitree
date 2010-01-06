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
make.bd(tree, times=NULL, sampling.f=NULL, unresolved=NULL)
make.yule(tree, times=NULL, sampling.f=NULL, unresolved=NULL)
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
    f(pars, prior=NULL, condition.surv=TRUE)
  }
  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of two parameters, in the order
    \code{lambda}, \code{mu}.
    \item \code{prior}: a valid prior.  See \code{\link{make.prior}} for
    more information.
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.
  }
}

\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bisse}} for state-dependent birth-death models.
}

\examples{
## Simulate a tree under a constant rates birth-death model and look at
## the maximum likelihood speciation/extinction parameters:
set.seed(1)
phy <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]
lik <- make.bd(phy)

## By default, optimisation gives a lambda close to 0.1 and extremely
## small mu:
fit <- find.mle(lik, c(.1, .03))
coef(fit)

## The above optimisation uses the algorithm \link{nlm} for
## compatibility with ape's \link{birthdeath}.  This can be slightly
## improved by using \link{optim} for the optimisation, which allows
## bounds to be specified:
fit.o <- find.mle(lik, c(.1, .03), method="optim", lower=0)
coef(fit.o)

logLik(fit.o) - logLik(fit) # slight improvement

## Special case methods are worked out for the Yule model, for which
## analytic solutions are available.  Compare a direct fit of the Yule
## model with one where mu is constrained to be zero:
lik.yule <- make.yule(phy)
lik.mu0 <- constrain(lik, mu ~ 0)

## The same to a reasonable tolerance:
fit.yule <- find.mle(lik.yule)
fit.mu0 <- find.mle(lik.mu0, .1)
all.equal(fit.yule[1:2], fit.mu0[1:2], tolerance=1e-6)

## There is no significant improvement in the fit by including the mu
## parameter (unsurprising as the ML value was zero)
anova(fit.o, yule=fit.yule)

## Optimisation can be done without conditioning on survival:
fit.nosurv <- find.mle(lik, c(.1, .03), method="optim", lower=0,
                       condition.surv=FALSE)
coef(fit.nosurv) # higher lambda than before

## Look at the marginal likelihoods, computed through MCMC (see
## \link{mcmc} for details, and increase \code{nsteps} for smoother
## plots [takes longer]).
samples <- mcmc(lik, fit$par, nsteps=2000,
                lower=c(-Inf, -Inf), upper=c(Inf, Inf), w=c(.1, .1),
                fail.value=-Inf, print.every=200)
samples$r <- with(samples, lambda - mu)

## Plot the profiles (see \link{profiles.plot}).
## The vertical lines are the simulated parameters, which match fairly
## well with the estimated ones.
col <- c("red", "blue", "green3")
profiles.plot(samples[c("lambda", "mu", "r")], col.line=col, las=1,
             ylab="Probability density", xlab="Parameter estimate")
legend("topright", c("lambda", "mu", "r"), fill=col)
abline(v=0, lty=2)
abline(v=c(.1, .03, .07), col=col)
}

\references{
  Nee S., May R.M., and Harvey P.H. 1994. The reconstructed
  evolutionary process. Philos.  Trans. R. Soc. Lond. B
  Biol. Sci. 344:305-311.
}
 
\author{Richard G. FitzJohn}
\keyword{models}
