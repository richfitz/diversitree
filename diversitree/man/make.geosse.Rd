\name{make.geosse}
\alias{make.geosse}
\alias{starting.point.geosse}

\title{Geographic State Speciation and Extinction Model}

\description{Prepare to run GeoSSE (Geographic State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
  make.geosse(tree, states, sampling.f=NULL, strict=TRUE,
    control=list())
  starting.point.geosse(tree, eps=0.5)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be 0
    (in both regions/widespread; AB), 1 or 2 (endemic to one region; A
    or B), or \code{NA} if the state is unknown.  This vector must have
    names that correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).}
  
  \item{sampling.f}{Vector of length 3 with the estimated proportion of
    extant species in states 0, 1 and 2 that are included in the
    phylogeny.  A value of \code{c(0.5, 0.75, 1)} means that half of
    species in state 0, three quarters of species in state 1, and all
    the species in state 2 are included in the phylogeny.  By default
    all species are assumed to be known.}
  
  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0, 1 and 2 only.  If \code{strict} is
    \code{TRUE} (the default), then the additional check is made that
    \emph{every} state is present.  The likelihood models tend to be
    poorly behaved where states are missing.}
  
  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}

  \item{eps}{Ratio of extinction to speciation rates to be used when
    choosing a starting set of parameters.  The procedure used is based
    on Magallon & Sanderson (2001).}
}

\details{
  \code{make.geosse} returns a function of class \code{geosse}.  The
  arguments and default values for this function are:

  \preformatted{
    f(pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
      intermediates=FALSE)
  }

  The arguments of this function are explained in \link{make.bisse}.
  The parameter vector \code{pars} is ordered \code{sA}, \code{sB},
  \code{sAB}, \code{xA}, \code{xB}, \code{dA}, \code{dB}.

  Unresolved clade methods are not available for GeoSSE.  
}


\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  \code{\link{make.bisse}} for further relevant examples.

  The help page for \code{\link{find.mle}} has further examples of ML
  searches on full and constrained BiSSE models.  Things work similarly
  for GeoSSE, just with different parameters.
}

\examples{
## Tree simulation under GeoSSE is currently only available through the
## SimTreeSDD program (\url{http://tigger.uic.edu/~eeg/code/code.html}).
## Here is an example, created with the values in "pars".

pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
data("geosse")

## The likelihood function
lik <- make.geosse(geosse.phy, geosse.phy$tip.state)

## See the data
statecols <- c("AB"="purple", "A"="blue", "B"="red")
plot(geosse.phy, tip.color=statecols[geosse.phy$tip.state+1], cex=0.5)

## With "true" parameter values
names(pars) <- argnames(lik)
lik(pars) # -297.6749

## A guess at a starting point.
p <- starting.point.geosse(geosse.phy)

## Start an ML search from this point (takes several minutes to run).
\dontrun{
fit <- find.mle(lik, p, method="subplex")
logLik(fit) # -295.0275

## Compare with sim values.
rbind(real=pars, estimated=round(coef(fit), 2))

## A model with constraints on the dispersal rates.
lik.d <- constrain(lik, dA ~ dB)
fit.d <- find.mle(lik.d, p[-7])
logLik(fit.d) # -296.5956

## A model with constraints on the speciation rates.
lik.s <- constrain(lik, sA ~ sB, sAB ~ 0)
fit.s <- find.mle(lik.s, p[-c(2,3)])
logLik(fit.s) # -302.5542
}

## "Skeletal tree" sampling is supported.  For example, if your tree
## includes all AB species, half of A species, and a third of B species,
## create the likelihood function like this:
lik.f <- make.geosse(geosse.phy, geosse.phy$tip.state,
                     sampling.f=c(1, 0.5, 1/3))

## If you have external evidence that the base of your tree must have
## been in state 1, say (endemic to region A), you can fix the root 
## when computing the likelihood, like this:
lik(pars, root=ROOT.GIVEN, root.p=c(0,1,0))
}

\references{
  FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
  trait-dependent speciation and extinction rates from incompletely
  resolved phylogenies. Syst. Biol. 58:595-611.
  
  Goldberg E.E., Lancaster L.T., and Ree R.H. 2011. Phylogenetic
  inference of reciprocal effects between geographic range evolution and
  diversification. Syst. Biol. 60:451-465.

  Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating a binary
  character's effect on speciation and extinction. Syst. Biol.
  56:701-710.

  Magallon S. and Sanderson M.J. 2001. Absolute diversification rates
  in angiospem clades. Evol. 55:1762-1780.
}

\author{Emma E. Goldberg}
\keyword{models}
