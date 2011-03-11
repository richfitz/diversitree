\name{make.musse.td}
\alias{make.musse.t}
\alias{make.musse.td}

\title{Multiple State Speciation and Extinction Model: Time Dependent
  Models}

\description{Create a likelihood function for a MuSSE model where
  different chunks of time have different parameters.  This code is
  experimental!}

\usage{
make.musse.td(tree, states, k, n.epoch, sampling.f=NULL,
              strict=TRUE, control=list())
make.musse.t(tree, states, k, functions, sampling.f=NULL,
              strict=TRUE, control=list())
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}

  \item{states}{A vector of character states, each of which must be an
    integer between 1 and \code{k}.  This vector must have names that
    correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips corresponding to unresolved
    clades, the state should be \code{NA}.}

  \item{k}{The number of states.}

  \item{n.epoch}{Number of epochs.  1 corresponds to plain MuSSE, so
    this will generally be an integer at least 2.}

  \item{functions}{A named list of functions of time.  See details.}  
  
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
}

\details{
  Please see \code{\link{make.bisse}} for further details.
}

\examples{
## Here we will start with the tree and three-state character set from
## the make.musse example.  This is a poorly contrived example.
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

## Suppose we want to see if diversification is different in the most
## recent 3 time units, compared with the rest of the tree (yes, this is
## a totally contrived example!):
plot(phy)
axisPhylo()
abline(v=max(branching.times(phy)) - 3, col="red", lty=3)

## For comparison, make a plain MuSSE likelihood function
lik.m <- make.musse(phy, phy$tip.state, 3)

## Create the time-dependent likelihood function.  The final argument
## here is the number of 'epochs' that are allowed.  Two epochs is one
## switch point.
lik.t <- make.musse.td(phy, phy$tip.state, 3, 2)

## The switch point is the first argument.  The remaining 12 parameters
## are the MuSSE parameters, with the first 6 being the most recent
## epoch.
argnames(lik.t)

pars.t <- c(3, pars, pars)
names(pars.t) <- argnames(lik.t)

## Calculations are identical to a reasonable tolerance:
lik.m(pars) - lik.t(pars.t)

## It will often be useful to constrain the time as a fixed quantity.
lik.t2 <- constrain(lik.t, t.1 ~ 3)

## Parameter estimation under maximum likelihood.  This is marked "don't
## run" because the time-dependent fit takes a few minutes.
\dontrun{
## Fit the MuSSE ML model
fit.m <- find.mle(lik.m, pars)

## And fit the MuSSE/td model
fit.t <- find.mle(lik.t2, pars.t[argnames(lik.t2)],
                  control=list(maxit=20000))

## Compare these two fits with a likelihood ratio test (lik.t2 is nested
## within lik.m)
anova(fit.m, td=fit.t)
}
}

\author{Richard G. FitzJohn}
\keyword{models}
