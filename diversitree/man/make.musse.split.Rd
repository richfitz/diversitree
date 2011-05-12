\name{make.musse.split}
\alias{make.musse.split}

\title{Multiple State Speciation and Extinction Model: Split Models}

\description{Create a likelihood function for a MuSSE model where the
  tree is partitioned into regions with different parameters.}

\usage{
make.musse.split(tree, states, k, nodes, split.t,
                 sampling.f=NULL, strict=TRUE, control=list())
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  
  \item{states}{A vector of character states, each of which must be an
    integer between 1 and \code{k}.  This vector must have names that
    correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips corresponding to unresolved
    clades, the state should be \code{NA}.}
  
  \item{k}{The number of states.}  
  
  \item{nodes}{Vector of nodes that will be split (see Details).}
  
  \item{split.t}{Vector of split times, same length as \code{nodes} (see
    Details).}
  
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

\details{Branching times can be controlled with the \code{split.t}
argument.  If this is \code{Inf}, split at the base of the branch (as in
MEDUSA).  If \code{0}, split at the top (closest to the present, as in
the new option for MEDUSA).  If \code{0 < split.t < Inf} then we split
at that time on the tree (zero is the present, with time growing
backwards).

TODO: Describe \code{nodes} and \code{split.t} here.
}

\examples{
## This example picks up from the tree used in the ?make.musse example.

## First, simulate the tree:
set.seed(2)
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
phy <- tree.musse(pars, 30, x0=1)

## Here is the phylogeny, with true character history superposed:
h <- history.from.sim.discrete(phy, 1:3)
plot(h, phy, show.node.label=TRUE, font=1, cex=.75, no.margin=TRUE)

## Here is a plain MuSSE function for later comparison:
lik.m <- make.musse(phy, phy$tip.state, 3)
lik.m(pars) # -110.8364

## Split this phylogeny at three points: nd16 and nd25, splitting it
## into three chunks
nodes <- c("nd16", "nd25")
nodelabels(node=match(nodes, phy$node.label) + length(phy$tip.label),
           pch=19, cex=2, col="#FF000099")

## To make a split BiSSE function, pass the node locations and times
## in.  Here, we'll use 'Inf' as the split time to mimick MEDUSA's
## behaviour of placing the split at the base of the branch subtended by
## a node.
lik.s <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf)

## The parameters must be a list of the same length as the number of
## partitions.  Partition '1' is the root partition, and partition i is
## the partition rooted at the node[i-1]:
argnames(lik.s)

## Because we have two nodes, there are three sets of parameters.
## Replicate the original list to get a starting point for the analysis:
pars.s <- rep(pars, 3)
names(pars.s) <- argnames(lik.s)

lik.s(pars.s) # -110.8364

## This is basically identical (to acceptable tolerance) to the plain
## MuSSE version:
lik.s(pars.s) - lik.m(pars)

## The resulting likelihood function can be used in ML analyses with
## find.mle.  However, because of the large number of parameters, this
## may take some time (especially with as few species as there are in
## this tree - getting convergence in a reasonable number of iterations
## is difficult).
\dontrun{
fit.s <- find.mle(lik.s, pars.s, control=list(maxit=20000))
}

## Bayesian analysis also works, using the mcmc function.  Given the
## large number of parameters, priors will be essential, as there will
## be no signal for several parameters.  Here, I am using an exponential
## distribution with a mean of twice the state-independent
## diversification rate.
\dontrun{
prior <- make.prior.exponential(1/(-2*diff(starting.point.bd(phy))))
samples <- mcmc(lik.s, pars.s, 100, prior=prior, w=1, print.every=10)
}
}

\author{Richard G. FitzJohn}
\keyword{models}
