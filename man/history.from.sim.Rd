\name{history.from.sim}
\alias{history.from.sim.discrete}

\title{Extract Character Histories From Simulations}

\usage{
history.from.sim.discrete(phy, states)
}

\arguments{
  \item{phy}{A phylogeny produced by \code{\link{tree.bisse}}.}
  \item{states}{Possible states.  For \code{\link{tree.bisse}} this
    should be \code{0:1}.}
}

\description{This function extracts a history object from a simulated
  phylogeny produced by \code{\link{tree.bisse}}.}

\keyword{models}
\author{Richard G. FitzJohn}
