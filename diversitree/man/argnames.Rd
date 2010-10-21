\name{argnames}
\alias{argnames}
\alias{argnames<-}
\alias{argnames.bisse}
\alias{argnames.geosse}
\alias{argnames.constrained}
\alias{argnames<-.bisse}
\alias{argnames<-.geosse}
\alias{argnames<-.constrained}

\title{Argument Names for Vector-Argument Functions}

\description{Functions to get and set \dQuote{argument names} for
  functions that take vectorised arguments.  For example, the likelihood
  function returned by \code{make.bisse} takes a vector of six these
  functions can be used to get the canonical names for these six
  parameters, and also to set them to something more memorable.  These
  names are used by the \code{\link{constrain}} function to specify
  submodels.
}

\usage{
argnames(x, ...)
argnames(x) <- value
\method{argnames}{bisse}(x, ...)
\method{argnames}{bisse}(x) <- value
\method{argnames}{constrained}(x, ...)
\method{argnames}{constrained}(x) <- value
}

\arguments{
  \item{x}{A function taking a vector of parameters as its first
    argument.}
  \item{value}{Vector of names to set the argument names to.}
  \item{...}{Ignored arguments to future methods.}
}

\details{Methods exist for all models: \code{bisse}, \code{geosse}, \code{bd},
  \code{yule}, \code{mk2}, and \code{mkn}.  These are particulary useful
  for \code{mkn} as the number of parameters for the Q matrix can be
  very large.}

% TODO: Need an example showing a three-state mkn.  This can probably be
% copied from the make.mkn file once that is added.
\examples{
## Same example likelihood function as for \link{make.bisse}:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)
f <- make.bisse(phy, phy$tip.state)

argnames(f) # Canonical argument names (set by default)
## Names that might be more informative for a tall/short state
argnames(f) <- c("l.tall", "l.short", "m.tall", "m.short",
                 "q.tall.short", "q.short.tall")
argnames(f)

}

\keyword{programming}
\author{Richard G. FitzJohn}
