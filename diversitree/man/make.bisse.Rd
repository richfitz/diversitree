\name{make.bisse}
\alias{make.bisse}

\title{Binary State Speciation and Extinction}

\description{Prepare to run BiSSE (Binary State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.
}

\usage{
make.bisse(tree, states, unresolved=NULL, sampling.f=NULL, nt.extra=10)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{states}{A vector of character states, each of which must be 0 or
    1.  This vector must have names that correspond to the tip labels in
    the phylogenetic tree (\code{tree$tip.label}).  For tips
    corresponding to unresolved clades, the state should be \code{NA}.}
  \item{unresolved}{Unresolved clade information: see section below for
    structure.}
  \item{sampling.f}{Vector of length 2 with the estimated proportion of
    extant species in state 0 and 1 that are included in the phylogeny.
    A value of \code{c(0.5, 0.75)} means that half of species in state 0
    and three quarters of species in state 1 are included in the
    phylogeny.  By default all species are assumed to be known}
  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}
}

\details{
  \code{make.bisse} returns a function of class \code{bisse}.  This
  function has argument list (and default values)

  \preformatted{
    f(pars, prior=NULL, root=ROOT.OBS,
      condition.surv=TRUE, intermediates=FALSE,
      root.p0=NA, root.p1=NA, fail.value=NULL)
    }
  
    The arguments are interpreted as
    \itemize{
      \item \code{pars} A vector of six parameters, in the order \code{lambda0},
      \code{lambda1}, \code{mu0}, \code{mu1}, \code{q01}, \code{q10}.
      \item \code{prior} A vector of rates for an exponential prior.  The prior
      probability distribution has probability density
      \deqn{r_i e^{-r_i x_i}}{r_i*exp(-r_i*x_i)}
      where the \eqn{i} denotes the \eqn{i}th parameter.  There is no
      ability to specify a form of prior other than an exponential
      function here (yet), and this part of the inteface is subject to
      change!
      \item \code{root}: Behaviour at the root (see Maddison et al. 2007,
      FitzJohn et al. 2009).  The possible options are
      \itemize{
	\item \code{diversitree:::ROOT.FLAT}: A flat prior, weighting
 	\eqn{D_0}{D0} and \eqn{D_1}{D1} equally.
	\item \code{diversitree:::ROOT.EQUI}: Use the equilibrium distribution
 	of the model, as described in Maddison et al. (2007).
	\item \code{diversitree:::ROOT.OBS}: Weight \eqn{D_0}{D0} and
 	\eqn{D_1}{D1} by their relative probability of observing the
 	data, following FitzJohn et al. 2009:
 	\deqn{D = D_0*\frac{D_0}{D_0 + D_1} + D_1*\frac{D_1}{D_0 + D_1}}{
 	  D = D0 * D0/(D0 + D1) + D1 * D1/(D0 + D1)}
	\item \code{diversitree:::ROOT.GIVEN}: Root will be in state 0
	with probability \code{root.p0}, or in state 1 with probability
	\code{root.p1}.  Specify either of these, but just one.
	\item \code{diversitree:::ROOT.BOTH}: Don't do anything at the root,
 	and return both values.  (Note that this will not give you a
 	likelihood!).
      }
      \item \code{intermediates}: Add intermediates to the returned value as
      attributes:
      \itemize{
       \item \code{cache}: Cached tree traversal information.
       \item \code{intermediates}: Mostly branch end information.
       \item \code{vals}: Root \eqn{D} values.
       \item \code{logComp}: Underflow compensation information.
     }
     At this point, you will have to poke about in the source for more
     information on these.
   }
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
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration.
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=30, x0=0)

plot(phy)

f <- make.bisse(phy, phy$tip.state)
f(pars)
}

\author{Richard G. FitzJohn}
\keyword{models}

