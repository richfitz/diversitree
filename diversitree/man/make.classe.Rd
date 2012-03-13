\name{make.classe}
\alias{make.classe}
\alias{starting.point.classe}

\title{Cladogenetic State change Speciation and Extinction Model}

\description{Prepare to run ClaSSE (Cladogenetic State change
  Speciation and Extinction) on a phylogenetic tree and character
  distribution.  This function creates a likelihood function that can be
  used in \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
  make.classe(tree, states, k, sampling.f=NULL, strict=TRUE,
              control=list())
  starting.point.classe(tree, k, eps=0.5)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}

  \item{states}{A vector of character states, each of which must be an
    integer between 1 and \code{k}.  This vector must have names that
    correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).}

  \item{k}{The number of states.  (The maximum now is 31, but that can
    easily be increased if necessary.)}
  
  \item{sampling.f}{Vector of length \code{k} where \code{sampling.f[i]}
    is the proportion of species in state \code{i} that are present in
    the phylogeny.  A value of \code{c(0.5, 0.75, 1)} means that half of
    species in state 1, three quarters of species in state 2, and all
    species in state 3 are included in the phylogeny.  By default all
    species are assumed to be known}

  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are integers on \code{1:k}.  If \code{strict} is
    \code{TRUE} (the default), then the additional check is made that
    \emph{every} state is present.  The likelihood models tend to be
    poorly behaved where states are missing, but there are cases
    (missing intermediate states for meristic characters) where allowing
    such models may be useful.}

  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}

  \item{eps}{Ratio of extinction to speciation rates to be used when
    choosing a starting set of parameters.  The procedure used is based
    on Magallon & Sanderson (2001).}
}

\details{
  The ClaSSE model with \code{k = 2} is equivalent to but a different
  parameterization than the \link[=make.bisseness]{BiSSE-ness} model.
  The \link[=make.geosse]{GeoSSE} model can be constructed from ClaSSE
  with \code{k = 3}; see the example below.

  \code{make.classe} returns a function of class \code{classe}.  The
  arguments and default values for this function are:

  \preformatted{
    f(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
      intermediates=FALSE)
  }

  The arguments are of this function are explained in \link{make.bisse}.
  The speciation rate parameters are lambda_ijk, ordered with k changing
  fastest and insisting on j < k.

  With more than 9 states, lambda_ijk and q_ij can be ambiguous (e.g. is
  q113 1->13 or 11->3?).  To avoid this, the numbers are zero padded (so
  that the above would be q0113 or q1103 for 1->13 and 11->3
  respectively).  It might be easier to rename the arguments in practice
  though.  More human-friendly handling of large speciation rate arrays
  is in the works.

  \code{starting.point.classe} produces a first-guess set of parameters,
    ignoring character states.

  Unresolved clade methods are not available for ClaSSE.

  Tree simulation methods are not yet available for ClaSSE.
}

\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, and \code{\link{mcmc}} for MCMC
  integration.  The help page for \code{\link{find.mle}} has further
  examples of ML searches on full and constrained BiSSE models.  Things
  work similarly for ClaSSE, just with different speciation parameters.
  	
  \code{\link{make.bisse}}, \code{\link{make.bisseness}},
  \code{\link{make.geosse}}, \code{\link{make.musse}} for similar models
  and further relevant examples.
}

\examples{
## GeoSSE equivalence
data("geosse")

lik.g <- make.geosse(geosse.phy, geosse.phy$tip.state)
pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars.g) <- argnames(lik.g)

lik.c <- make.classe(geosse.phy, geosse.phy$tip.state+1, 3)
pars.c <- 0 * starting.point.classe(geosse.phy, 3)
pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
pars.c['lambda123'] <-  pars.g['sAB']
pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
pars.c['q21'] <- pars.g['dA']
pars.c['q31'] <- pars.g['dB']

lik.g(pars.g)   # -309.9812
lik.c(pars.c)   # -309.9812
}

\references{
  FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
  trait-dependent speciation and extinction rates from incompletely
  resolved phylogenies. Syst. Biol. 58:595-611.
  
  Goldberg E.E. and Igic B. Tempo and mode in plant breeding system
  evolution. In review.

  Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary
  character's effect on speciation and extinction. Syst. Biol.
  56:701-710.

  Magallon S. and Sanderson M.J. 2001. Absolute diversification rates in
  angiospem clades. Evol. 55:1762-1780.

  Magnuson-Ford, K., and Otto, S.P. 2012. Linking the investigations of
  character evolution and species diversification.  American Naturalist,
  in press.
}

\author{Emma E. Goldberg}
\keyword{models}
