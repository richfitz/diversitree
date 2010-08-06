\name{make.bisse}
\alias{make.bisse}
\alias{starting.point.bisse}

\title{Binary State Speciation and Extinction Model}

\description{Prepare to run BiSSE (Binary State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.bisse(tree, states, unresolved=NULL, sampling.f=NULL, nt.extra=10,
           safe=FALSE, strict=TRUE)
starting.point.bisse(tree, q.div=5, yule=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{states}{A vector of character states, each of which must be 0 or
    1, or \code{NA} if the state is unknown.  This vector must have
    names that correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips
    corresponding to unresolved clades, the state should be \code{NA}.}
  \item{unresolved}{Unresolved clade information: see section below for
    structure.}
  \item{sampling.f}{Vector of length 2 with the estimated proportion of
    extant species in state 0 and 1 that are included in the phylogeny.
    A value of \code{c(0.5, 0.75)} means that half of species in state 0
    and three quarters of species in state 1 are included in the
    phylogeny.  By default all species are assumed to be known.}
  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}
  \item{safe}{Use a safer, slower, integrator?  Leave this alone unless
    crashes are happening.}
  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0 and 1 only.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present.  The likelihood models tend to be poorly behaved
    where states are missing.}
  \item{q.div}{Ratio of diversification rate to character change rate.
    Eventually this will be changed to allow for Mk2 to be used for
    estimating q parameters.}
  \item{yule}{Logical: should starting parameters be Yule estimates
    rather than birth-death estimates?}
}

\details{
  \code{make.bisse} returns a function of class \code{bisse}.  This
  function has argument list (and default values)

  \preformatted{
    f(pars, prior=NULL, condition.surv=TRUE, root=ROOT.OBS,
    root.p=NULL, intermediates=FALSE)
  }
  
  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of six parameters, in the order
    \code{lambda0}, \code{lambda1}, \code{mu0}, \code{mu1},
    \code{q01}, \code{q10}.
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.
    \item \code{root}: Behaviour at the root (see Maddison et al. 2007,
    FitzJohn et al. 2009).  The possible options are
    \itemize{
      \item \code{ROOT.FLAT}: A flat prior, weighting
      \eqn{D_0}{D0} and \eqn{D_1}{D1} equally.
      \item \code{ROOT.EQUI}: Use the equilibrium distribution
      of the model, as described in Maddison et al. (2007).
      \item \code{ROOT.OBS}: Weight \eqn{D_0}{D0} and
      \eqn{D_1}{D1} by their relative probability of observing the
      data, following FitzJohn et al. 2009:
      \deqn{D = D_0\frac{D_0}{D_0 + D_1} + D_1\frac{D_1}{D_0 + D_1}}{
 	D = D0 * D0/(D0 + D1) + D1 * D1/(D0 + D1)}
      \item \code{ROOT.GIVEN}: Root will be in state 0
      with probability \code{root.p[1]}, and in state 1 with
      probability \code{root.p[1]}.
      \item \code{ROOT.BOTH}: Don't do anything at the root,
      and return both values.  (Note that this will not give you a
      likelihood!).
    }
    \item \code{root.p}: Root weightings for use when
    \code{root=ROOT.GIVEN}.  \code{sum(root.p)} should equal 1.
    \item \code{intermediates}: Add intermediates to the returned value as
    attributes:
    \itemize{
      \item \code{cache}: Cached tree traversal information.
      \item \code{intermediates}: Mostly branch end information.
      \item \code{vals}: Root \eqn{D} values.
    }
    At this point, you will have to poke about in the source for more
    information on these.
  }
  
  \code{starting.point.bisse} produces a heuristic starting point to
  start from, based on the character-independent birth-death model.  You
  can probably do better than this; see the vignette, for example.
  \code{bisse.starting.point} is the same code, but deprecated in favour
  of \code{starting.point.bisse} - it will be removed in a future
  version.
}

\section{Unresolved clade information}{
  This must be a \code{\link{data.frame}} with at least the four columns
  \itemize{
    \item \code{tip.label}, giving the name of the tip to which the data
    applies
    \item \code{Nc}, giving the number of species in the clade
    \item \code{n0}, \code{n1}, giving the number of species known to be
    in state 0 and 1, respectively.
  }
  These columns may be in any order, and additional columns will be
  ignored.  (Note that column names are case sensitive).

  An alternative way of specifying unresolved clade information is to
  use the function \code{\link{make.clade.tree}} to construct a tree
  where tips that represent clades contain information about which
  species are contained within the clades.  With a \code{clade.tree},
  the \code{unresolved} object will be automatically constructed from
  the state information in \code{states}.  (In this case, \code{states}
  must contain state information for the species contained within the
  unresolved clades.)
}


\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bd}} for state-independent birth-death models.

  The help pages for \code{\link{find.mle}} has further examples of ML
  searches on full and constrained BiSSE models.
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)

## Here is the 52 species tree with the true character history coded.
## Red is state '1', which has twice the speciation rate of black (state
## '0').
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)

lik <- make.bisse(phy, phy$tip.state)
lik(pars) # -159.71

## Heuristic guess at a starting point, based on the constant-rate
## birth-death model (uses \link{make.bd}).
p <- starting.point.bisse(phy)

%% This is code for better getting at the q's
%% system.time(find.mle(constrain(make.mk2(phy, phy$tip.state), q10 ~
%% q01), .1)$par)

## Start an ML search from this point.  This takes some time (~7s)
fit <- find.mle(lik, p, method="subplex")
logLik(fit) # -158.69

## The estimated parameters aren't too far away from the real ones, even
## with such a small tree
rbind(real=pars,
      estimated=round(coef(fit), 2))

## Test a constrained model where the speciation rates are set equal
## (takes ~4s).
lik.l <- constrain(lik, lambda1 ~ lambda0)
fit.l <- find.mle(lik.l, p[-1], method="subplex")
logLik(fit.l) # -158.74

## Despite the difference in the estimated parameters, there is no
## statistical support for this difference:
anova(fit, equal.lambda=fit.l)

## Run an MCMC.  Because we are fitting six parameters to a tree with
## only 50 species, priors will be needed.  I will use an exponential
## prior with rate 1/(2r), where r is the character independent
## diversificiation rate:
prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))

## This takes quite a while to run, so is not run by default
\dontrun{
tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,
            lower=rep(0, 6), upper=rep(Inf, 6), w=rep(.1, 6))
w <- diff(sapply(tmp[2:7], range))
samples <- mcmc(lik, fit$par, nsteps=1000, prior=prior,
                lower=rep(0, 6), upper=rep(Inf, 6), w=rep(.1, 6))

## See \link{profiles.plot} for more information on plotting these
## profiles.
col <- c("blue", "red")
profiles.plot(samples[c("lambda0", "lambda1")], col.line=col, las=1,
              xlab="Speciation rate",
              ylab="Posterior probability density")
legend("topright", c("lambda0", "lambda1"), col=col, lty=1)
}

## BiSSE reduces to the birth-death model and Mk2 when diversification
## is state independent (i.e., lambda0 ~ lambda1 and mu0 ~ mu1).
lik.mk2 <- make.mk2(phy, phy$tip.state)
lik.bd <- make.bd(phy)

## 1. BiSSE / Birth-Death
## Set the q01 and q10 parameters to arbitrary numbers (need not be
## symmetric), and constrain the lambdas and mus to be the same for each
## state.  The likelihood function now has just two parameters and
## will be proprtional to Nee's birth-death based likelihood:
lik.bisse.bd <- constrain(lik,
                          lambda1 ~ lambda0, mu1 ~ mu0,
                          q01 ~ .01, q10 ~ .02)
pars <- c(.1, .03)
## These differ by -167.3861 for both parameter sets:
lik.bisse.bd(pars)   - lik.bd(pars)
lik.bisse.bd(2*pars) - lik.bd(2*pars)

## 2. BiSSE / Mk2
## Same idea as above: set all diversification parameters to arbitrary
## values (but symmetric this time):
lik.bisse.mk2 <- constrain(lik,
                           lambda0 ~ .1, lambda1 ~ .1,
                           mu0 ~ .03, mu1 ~ .03)
## Differ by -150.4740 for both parameter sets.
lik.bisse.mk2(pars)   - lik.mk2(pars)
lik.bisse.mk2(2*pars) - lik.mk2(2*pars)

## 3. Sampled BiSSE / Birth-Death
## Pretend that the tree is only .6 sampled:
lik.bd2 <- make.bd(phy, sampling.f=.6)
lik.bisse2 <- make.bisse(phy, phy$tip.state, sampling.f=c(.6, .6))
lik.bisse2.bd <- constrain(lik.bisse2,
                           lambda1 ~ lambda0, mu1 ~ mu0,
                           q01 ~ .01, q10 ~ .01)

## Difference of -167.2876
lik.bisse2.bd(pars)   - lik.bd2(pars)
lik.bisse2.bd(2*pars) - lik.bd2(2*pars)

## 4. Unresolved clade BiSSE / Birth-Death
set.seed(1)
unresolved <- data.frame(tip.label=I(sample(phy$tip.label, 5)),
                         Nc=sample(2:10, 5), n0=0, n1=0)
unresolved.bd <- structure(unresolved$Nc, names=unresolved$tip.label)
lik.bisse3 <- make.bisse(phy, phy$tip.state, unresolved)
lik.bisse3.bd <- constrain(lik.bisse3,
                           lambda1 ~ lambda0, mu1 ~ mu0,
                           q01 ~ .01, q10 ~ .01)
lik.bd3 <- make.bd(phy, unresolved=unresolved.bd)

## Difference of -167.1523
lik.bisse3.bd(pars) - lik.bd3(pars)
lik.bisse3.bd(pars*2) - lik.bd3(pars*2)
}

\references{
  FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
  trait-dependent speciation and extinction rates from incompletely
  resolved phylogenies. Syst. Biol. 58:595-611.
  
  Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating
  a binary character's effect on speciation and
  extinction. Syst. Biol. 56:701-710.

  Nee S., May R.M., and Harvey P.H. 1994. The reconstructed
  evolutionary process. Philos.  Trans. R. Soc. Lond. B
  Biol. Sci. 344:305-311.
}

\author{Richard G. FitzJohn}
\keyword{models}

