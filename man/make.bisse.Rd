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
           strict=TRUE, control=list())
starting.point.bisse(tree, q.div=5, yule=FALSE)
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

  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}

  \item{control}{List of control parameters for the ODE solver.  See
    details below.}
  
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
    f(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
      intermediates=FALSE)
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

\section{ODE solver control}{
  The differential equations that define the
  BiSSE model are solved numerically using ODE solvers from the GSL
  library or deSolve's LSODA.  The \code{control} argument to
  \code{make.bisse} controls the behaviour of the integrator.  This is a
  list that may contain elements:
  
  \itemize{
    \item{\code{tol}: Numerical tolerance used for the calculations.
      The default value of \code{1e-8} should be a reasonable trade-off
      between speed and accuracy.  Do not expect too much more than this
      from the abilities of most machines!}

    \item{\code{eps}: A value that when the sum of the D values drops
      below, the integration results will be discarded and the
      integration will be attempted again (the second-chance integration
      will divide a branch in two and try again, recursively until the
      desired accuracy is reached).  The default value of \code{0} will
      only discard integration results when the parameters go negative.
      However, for some problems more restrictive values (on the order
      of \code{control$tol}) will give better stability.}

    \item{\code{backend}: Select the solver.  The three options here are
      \itemize{
	\item{\code{gslode}: (the default).  Use the GSL solvers, by
	  default a Runge Kutta Kash Carp stepper.}
	\item{\code{deSolve}: Use the LSODA solver from the
	  \code{deSolve} package.  This is quite a bit slower at the
	  moment.}
    }}
  }

  \code{deSolve} is the only supported backend on Windows.
}

\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bd}} for state-independent birth-death models.

  The help pages for \code{\link{find.mle}} has further examples of ML
  searches on full and constrained BiSSE models.
}

\examples{
## Due to a change in sample() behaviour in newer R it is necessary to
## use an older algorithm to replicate the previous examples
if (getRversion() >= "3.6.0") {
  RNGkind(sample.kind = "Rounding")
}
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
