\name{plot.history}
\alias{plot.history}
\title{Plot Character History}

\description{Both stochastic character mapping and simulation may create
  character histories.  This function plots these histories}

\usage{
\method{plot}{history}(x, phy, cols=seq_along(states),
             states=x$states,
             xlim=NULL, ylim=NULL, show.tip.label=TRUE,
             show.node.label=FALSE, show.tip.state=TRUE,
             show.node.state=TRUE, no.margin=FALSE, cex=1, font=3,
             srt=0, adj=0, label.offset=NA, lwd=1, ...)
}

\arguments{
  \item{x}{An object of class \code{history.discrete} containing a
    discrete character history.  This could be made by
    \code{\link{history.from.sim.discrete}}.}
  \item{phy}{The phylogeny used to generate the history.  Few checks are
    made to make sure that this is really correct, and all manner of
    terrible things might happen if these are not compatible.  This may
    change in future.}
  \item{cols}{A vector of colours.}
  \item{states}{The different state types.  Probably best to leave
    alone.}
  \item{xlim}{Plot x-limits (optional).}
  \item{ylim}{Plot y-limits (optional).}
  \item{show.tip.label}{Logical: show the species tip labels?}
  \item{show.node.label}{Logical: show the species node labels?}
  \item{show.tip.state}{Logical: draw a symbol at the tips to indicate
    tip state?}
  \item{show.node.state}{Logical: draw a symbol at the nodes to indicate
    node state?}
  \item{no.margin}{Supress drawing of margins around the plot}
  \item{cex}{Font and symbol scaling factor.}
  \item{font}{Font used for tip and node labels (see
    \code{\link{par}}).}
  \item{srt}{String rotation for tip and node labels.}
  \item{adj}{Label \code{adj}ustment (see \code{\link{par}}).}
  \item{label.offset}{Horizontal offset of tip and node labels, in
    branch length units.}
  \item{lwd}{Line width}
  \item{...}{Additional arguments (currently ignored)}
}

\details{
  This attempts to be as compatible with \code{ape}'s plotting functions
  as possible, but currently implements only right-facing cladegrams.
}

\examples{
## Simulate a tree, but retain extinct species.
pars <- c(.1, .2, .03, .04, 0.05, 0.1) # BiSSE pars
set.seed(2)
phy <- tree.bisse(pars, 20, x0=0, include.extinct=TRUE)

## Create a 'history' from the information produced by the simulation
## and plot this
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, cex=.7)

## Prune the extinct taxa.
phy2 <- prune(phy)

## The history must be recreated for this pruned tree:
h2 <- history.from.sim.discrete(phy2, 0:1)
plot(h2, phy2, cex=.7)
}

\keyword{hplot}
\author{Richard G. FitzJohn}
