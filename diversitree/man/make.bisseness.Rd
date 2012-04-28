\name{make.bisseness}
\alias{make.bisseness}

\title{Binary State Speciation and Extinction (Node Enhanced State Shift) Model}

\description{Prepare to run BiSSE-ness (Binary State Speciation and
  Extinction (Node Enhanced State Shift)) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.bisseness(tree, states, unresolved=NULL, sampling.f=NULL,
               nt.extra=10, strict=TRUE, control=list())
}

\arguments{
  \item{tree}{An ultrametric bifurcating phylogenetic tree, in
    \code{ape} \dQuote{phylo} format.}
  
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

  \item{nt.extra}{The number of "extra" species to include in the
    unresolved clade calculations.  This is in addition to the largest
    included unresolved clade.}

  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}
  
  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0 and 1 only.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present at least once in the tree.  The likelihood models
    tend to be poorly behaved where a state is not represented on the
    tree.}
}

\details{
  \code{make.bisse} returns a function of class \code{bisse}.  This
  function has argument list (and default values) [RICH:  Update to BiSSEness?]

  \preformatted{
    f(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
      intermediates=FALSE)
  }
  
  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of 10 parameters, in the order
    \code{lambda0}, \code{lambda1}, \code{mu0}, \code{mu1},
    \code{q01}, \code{q10}, \code{p0c}, \code{p0a}, \code{p1c}, \code{p1a}.
    \item \code{condition.surv} (logical): should the likelihood
    calculation condition on survival of two lineages and the speciation
    event subtending them?  This is done by default, following Nee et
    al. 1994.  For BiSSE-ness, equation (A5) in Magnuson-Ford and Otto
    describes how conditioning on survival alters the likelihood of
    observing the data.
    % Modification here:
    \item \code{root}: Behaviour at the root (see Maddison et al. 2007,
    FitzJohn et al. 2009).  The possible options are
    \itemize{
      \item \code{ROOT.FLAT}: A flat prior, weighting
      \eqn{D_0}{D0} and \eqn{D_1}{D1} equally.
      \item \code{ROOT.EQUI}: Use the equilibrium distribution
      of the model, as described in Maddison et al. (2007) using
      equation (A6) in Magnuson-Ford and Otto.
      \item \code{ROOT.OBS}: Weight \eqn{D_0}{D0} and
      \eqn{D_1}{D1} by their relative probability of observing the
      data, following FitzJohn et al. 2009:
      \deqn{D = D_0\frac{D_0}{D_0 + D_1} + D_1\frac{D_1}{D_0 + D_1}}{
 	D = D0 * D0/(D0 + D1) + D1 * D1/(D0 + D1)}
      \item \code{ROOT.GIVEN}: Root will be in state 0
      with probability \code{root.p[1]}, and in state 1 with
      probability \code{root.p[2]}.
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
}

\section{Unresolved clade information}{
  % Might be best to comment this out/delete until this is implemented?
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
  \code{\link{make.bisse}} for the model with no state change at nodes.

  \code{\link{tree.bisseness}} for simulating trees under the BiSSE-ness
  model.
  
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bd}} for state-independent birth-death models.

  The help pages for \code{\link{find.mle}} has further examples of ML
  searches on full and constrained BiSSE models.
}

\examples{
## First we simulat a 50 species tree, assuming cladogenetic shifts in 
## the trait (i.e., the trait only changes at speciation).
## Red is state '1', black is state '0', and we let red lineages
## speciate at twice the rate of black lineages.
## The simulation starts in state 0.
set.seed(3)
pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
phy <- tree.bisseness(pars, max.taxa=50, x0=0)
phy$tip.state

h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)

## This builds the likelihood of the data according to BiSSEness:
lik <- make.bisseness(phy, phy$tip.state)
## e.g., the likelihood of the true parameters is:
lik(pars) # -174.7954

## ML search:  First we make hueristic guess at a starting point, based
## on the constant-rate birth-death model assuming anagenesis (uses
## \link{make.bd}).
startp <- starting.point.bisse(phy)

## We then take the total amount of anagenetic change expected across
## the tree and assign half of this change to anagenesis and half to
## cladogenetic change at the nodes as a heuristic starting point:
t <- branching.times(phy)
tryq <- 1/2 * startp[["q01"]] * sum(t)/length(t)
p <- c(startp[1:4], startp[5:6]/2, p0c=tryq, p0a=0.5, p1c=tryq, p1a=0.5)

## Start an ML search from this point.  This takes some time (~12s)
fit <- find.mle(lik, p, method="subplex")
logLik(fit) # -174.0104

## Compare the fit to a constrained model that only allows the trait
## to change along a lineage (anagenesis).  This takes some time (~12s)
lik.no.clado <- constrain(lik, p0c ~ 0, p1c ~ 0)
fit.no.clado <- find.mle(lik.no.clado,p[argnames(lik.no.clado)])
logLik(fit.no.clado) # -174.0577

## This is consistent with what BiSSE finds:
likB <- make.bisse(phy, phy$tip.state)
fitB <- find.mle(likB, startp, method="subplex")
logLik(fitB) # -174.0576

## With only this 50-species tree, there is no statistical support
## for the more complicated BiSSE-ness model that allows cladogenesis:
anova(fit, no.clado=fit.no.clado)
## Note that anova() performs a likelihood ratio test here.

## If the above is repeated with max.taxa=250, BiSSE-ness rejects the
## constrained model in favor of one that allows cladogenetic change.

\dontrun{
## MCMC run: We use the ML estimate from the full model
## as a starting point.
##
## We shift all very small numbers up to 1e-4 to allow the derivatives
## to be calculated.
ml.start.pt <- pmax(coef(fit), 1e-4)

## Make exponential priors for the rate parameters and uniform priors
## for the cladogenetic change probability prarameters.
make.prior.exp_ness <- function(r, min=0, max=1) {
  function(pars) {
    sum(dexp(pars[1:6], rate=r, log=TRUE)) +
      sum(dunif(pars[7:10], min, max, log=TRUE))
  }
}

## Choosing the slice sampling parameter, w (affects speed):
library(numDeriv)
hess <- hessian(lik, ml.start.pt)
vcv <- -solve(hess)
sehess <- sqrt(abs(diag(vcv)))
w <- 2 * pmin(sehess, .2)

## Setting the priors
r <- log(length(phy$tip.label))/max(branching.times(phy))
prior <- make.prior.exp_ness(1/(2*r))
prior(ml.start.pt)

## Running the mcmc chain (only 10 steps are shown for illustration)
steps <- 10
set.seed(1) # For reproducibility
output <- mcmc(lik, ml.start.pt, nsteps=steps, w=w, prior=prior)

## Unresolved tip clade: Here we collapse one clade in the 50 species
## tree (involving sister species sp70 and sp71) and illustrate the use
## of BiSSEness with unresolved tip clades.
slimphy <- drop.tip(phy,c("sp71"))
states <- slimphy$tip.state[slimphy$tip.label]
states["sp70"] <- NA
unresolved <- data.frame(tip.label=c("sp70"), Nc=2, n0=2, n1=0)

## This builds the likelihood of the data according to BiSSEness:
lik.unresolved <- make.bisseness(slimphy, states, unresolved)
## e.g., the likelihood of the true parameters is:
lik.unresolved(pars) # -174.6575

## ML search from the heuristic starting point used above:
fit.unresolved <- find.mle(lik.unresolved, p, method="subplex")
logLik(fit.unresolved) # -173.9136
}
}


\references{
  FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
  trait-dependent speciation and extinction rates from incompletely
  resolved phylogenies. Syst. Biol. 58:595-611.
  
  Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating
  a binary character's effect on speciation and
  extinction. Syst. Biol. 56:701-710.

  Magnuson-Ford, K., and Otto, S.P. 2012. Linking the investigations of
  character evolution and species diversification.  American Naturalist,
  in press.

  Nee S., May R.M., and Harvey P.H. 1994. The reconstructed
  evolutionary process. Philos.  Trans. R. Soc. Lond. B
  Biol. Sci. 344:305-311.
}

\author{Karen Magnuson-Ford}
\keyword{models}

