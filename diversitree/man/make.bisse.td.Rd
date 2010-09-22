\name{make.bisse.td}
\alias{make.bisse.td}

\title{Binary State Speciation and Extinction Model: Time Chunks}

\description{Create a likelihood function for a BiSSE model where
  different chunks of time have different parameters.  This code is
  experimental!}

\usage{
make.bisse.td(tree, states, n.epoch, unresolved=NULL, sampling.f=NULL,
              nt.extra=10, strict=TRUE, safe=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be 0 or
    1, or \code{NA} if the state is unknown.  This vector must have
    names that correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips
    corresponding to unresolved clades, the state should be \code{NA}.}
  
  \item{n.epoch}{Number of epochs.  1 corresponds to plain BiSSE, so
    this will generally be an integer at least 2.}
  
  \item{unresolved}{Unresolved clade information: see
    \code{\link{make.bisse}}. (Currently this is not supported.)}
  
  \item{sampling.f}{Vector of length 2 with the estimated proportion of
    extant species in state 0 and 1 that are included in the phylogeny.
    See \code{\link{make.bisse}}.}
  
  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}

  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0 and 1 only.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present.  The likelihood models tend to be poorly behaved
    where states are missing.}
  
  \item{safe}{Use a safer, slower, integrator?  Leave this alone unless
    crashes are happening.}
}

\details{
  This builds a BiSSE likelihood function where different regions of
  time (epochs) have different parameter sets.  By default, all
  parameters are free to vary between epochs, so some constraining will
  probably be required to get reasonable answers.

  For \code{n} epochs, there are \code{n-1} time points; the first
  \code{n-1} elements of the likelihood's parameter vector are these
  points.  These are measured from the present at time zero, with time
  increasing towards the base of the tree.  The rest of the parameter
  vector are BiSSE parameters; the elements \code{n:(n+6)} are for the
  first epoch (closest to the present), elements \code{(n+7):(n+13)} are
  for the second epoch, and so on.
  
}


\examples{
set.seed(4)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
phy <- tree.bisse(pars, max.t=30, x0=0)

## Suppose we want to see if diversification is different in the most
## recent 3 time units, compared with the rest of the tree (yes, this is
## a totally contrived example!):
plot(phy)
axisPhylo()
abline(v=max(branching.times(phy)) - 3, col="red", lty=3)

## For comparison, make a plain BiSSE likelihood function
lik.b <- make.bisse(phy, phy$tip.state)

## Create the time-dependent likelihood function.  The final argument
## here is the number of 'epochs' that are allowed.  Two epochs is one
## switch point.
lik.t <- make.bisse.td(phy, phy$tip.state, 2)

## The switch point is the first argument.  The remaining 12 parameters
## are the BiSSE parameters, with the first 6 being the most recent
## epoch.
argnames(lik.t)

pars.t <- c(3, pars, pars)
names(pars.t) <- argnames(lik.t)

## Calculations are identical to a reasonable tolerance:
lik.b(pars) - lik.t(pars.t)

## It will often be useful to constrain the time as a fixed quantity.
lik.t2 <- constrain(lik.t, t.1 ~ 3)

## Parameter estimation under maximum likelihood.  This is marked "don't
## run" because the time-dependent fit takes a few minutes.
\dontrun{
## Fit the BiSSE ML model
fit.b <- find.mle(lik.b, pars)

## And fit the BiSSE/td model
fit.t <- find.mle(lik.t2, pars.t[argnames(lik.t2)],
                  control=list(maxit=20000))

## Compare these two fits with a likelihood ratio test (lik.t2 is nested
## within lik.b)
anova(fit.b, td=fit.t)
}
}

\author{Richard G. FitzJohn}
\keyword{models}
