\name{plot.history}
\alias{plot.history.discrete}
\title{Plot Character History}

\description{Both stochastic character mapping and simulation may create
  character histories.  This function plots these histories}

\usage{
\method{plot}{history.discrete}(x, phy, cols=c("blue", "red"),
                      xlim=NULL, ylim=NULL, show.tip.label=TRUE,
                      show.node.label=FALSE, show.tip.state=TRUE,
                      show.node.state=TRUE, cex=1, font=3, srt=0,
                      adj=0, label.offset=NULL, ...)
}

\arguments{
  \item{x}{An object of class \code{history.discrete} containing a
    discrete character history}
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
plot(h2, phy2, cex=.7, font=1)
}