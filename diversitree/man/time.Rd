\name{time}
\alias{constant.t}
\alias{linear.t}
\alias{exponential.t}
\alias{sigmoid.t}
\alias{stepf.t}
\alias{make.pars.t}

\title{Support Functions for Time-varying Models}

\description{Utility functions for working with time-dependent models.}

\usage{
constant.t(t, c)
linear.t(t, c, m)
exponential.t(t, y0, y1, r)
sigmoid.t(t, y0, y1, tmid, r)
stepf.t(t, y0, y1, tc)

make.pars.t(functions)
}

\arguments{
  \item{t}{Time}

  \item{c}{Constant.  For \code{linear.t} this is the value at t=0.}

  \item{m}{Slope}
  
  \item{y0}{y value at time 0 (\code{exp.t}) or very early times (limit
    as t tends to negative infinity in \code{sigmoid.t})}
    
  \item{y1}{y value approached as t tends to  negative infinity}

  \item{r}{Rate at which exponential decay occurs or sigmoid changes -
    higher values are steeper}

  \item{tmid}{Midpoint (inflection point) of sigmoid}

  \item{tc}{"Critical time" at which y changes from \code{y0} to {y1}.}

  \item{functions}{A list of functions with \code{t} as their first
    argument.}
}

\details{
  The linear function (\code{linear.t}) will go to zero wherever
  negative.  This may not always be desired, but is often useful for
  models of rates.
}

\author{Richard G. FitzJohn}
\keyword{models}
