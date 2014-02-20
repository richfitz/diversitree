\name{quasse-common}
\alias{constant.x}
\alias{make.linear.x}
\alias{sigmoid.x}
\alias{sigmoid2.x}
\alias{stepf.x}
\alias{noroptimal.x}
\alias{make.brownian.with.drift}

\title{Support Functions for QuaSSE Models}

\description{Utility functions for working with QuaSSE models.  These
  provide a minimal set of state-varying functions, suitable for use
  with \code{\link{make.quasse}}, and simulation assistance functions
  for use with \code{\link{tree.quasse}}.

  This is currently poorly explained!
}

\usage{
constant.x(x, c)
sigmoid.x(x, y0, y1, xmid, r)
stepf.x(x, y0, y1, xmid)
noroptimal.x(x, y0, y1, xmid, s2)

make.linear.x(x0, x1)

make.brownian.with.drift(drift, diffusion)
}

\arguments{
  \item{x}{Character state}

  \item{c}{Constant.}

  \item{y0}{y value at very small \code{x} (limit as \code{x} tends to
    negative infinity)}
  
  \item{y1}{y value at very large \code{x} (limit as \code{x} tends to
    infinity).  For \code{noroptimal.x}, this is the y value at
    \code{xmid}.}
  
  \item{xmid}{Midpoint (inflection point) of sigmoid or step function}

  \item{r}{Rate at which exponential decay occurs or sigmoid changes -
    higher values are steeper}

  \item{s2}{Variance of the normal distribution specified by
    \code{noroptimal.x}.}
  
  \item{x0}{Lower x limit for the linear function: y will take value at
    x0 for all x smaller than this}
  
  \item{x1}{Upper x limit for the linear function: y will take value at
    x1 for all x greater than this}

  \item{drift}{Rate of drift}

  \item{diffusion}{Rate of diffusion (positive)}
}

\details{
  The linear function returned by  (\code{make.linear.x}) will go to
  zero wherever negative.  This may not always be desired, but is
  required for valid likelihood calculations under QuaSSE.
}

\author{Richard G. FitzJohn}
\keyword{models}
