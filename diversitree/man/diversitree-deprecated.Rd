\name{diversitree-deprecated}
\alias{diversitree-deprecated}

\title{Deprecated and Defunct Functions in Diversitree}

\alias{starting.point}
\alias{bisse.starting.point}

\description{These functions have all been replaced by other functions.
  A warning is given when a deprecated function is used, and an error is
  given when a defunct function is used.}

\usage{
% Deprecated in 0.4-1
starting.point(tree, q.div=5)
bisse.starting.point(tree, q.div=5)
}

\details{
  Both \code{starting.point()} and \code{bisse.starting.point()}
  have been replaced by \code{starting.point.bisse()}.  This is to allow
  the other models to (eventually) have starting point generating
  functions, too.
}

\seealso{
  \code{\link{Defunct}}, \code{\link{Deprecated}}
}

\keyword{internal}
