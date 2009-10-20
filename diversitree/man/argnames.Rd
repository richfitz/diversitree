\name{argnames}
\alias{argnames}
\alias{argnames<-}
\alias{argnames.bisse}
\alias{argnames.constrained}
\alias{argnames<-.bisse}
\alias{argnames<-.constrained}

\title{Argument Names for Vector-Argument Functions}

\description{Functions to get and set \dQuote{argument names}
  for functions that take vectorised arguments.  The likelihood function
  returned by \code{make.bisse} takes a vector of six parameters.  This
  function can be used to get the canonical names for these six
  parameters, and also to set them to something more memorable.  These
  names are used by the \code{\link{constrain}} function to specify
  submodels.}

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


\examples{
## Same example likelihood function as for \link{make.bisse}:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=30, x0=0)
f <- make.bisse(phy, phy$tip.state)

argnames(f) # Canonical argument names (set by default)
## Names that might be more informative for a tall/short state
argnames(f) <- c("l.tall", "l.short", "m.tall", "m.short",
                 "q.tall.short", "q.short.tall")
argnames(f)
}

\keyword{programming}
