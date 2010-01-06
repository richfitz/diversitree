\name{protect}
\alias{protect}
\alias{invert}

\title{Protect Function Evaluations}

\description{Ensures non-failure and possibly finite-ness of a function
  evaluation.}

\usage{
protect(f)
invert(f)
}

\arguments{
  \item{f}{A function.}
}

\details{
  \code{protect} returns a function with arguments
  \preformatted{
    g(..., fail.value=NULL, finite=NULL)
  }
  The \code{...} arguments are all passed through to the underlying
  function \code{f}, \code{fail.value} contains the value to return in
  the event of a failure (e.g., an error occuring).  If \code{finite} is
  \code{TRUE}, then \code{fail.value} is also returned where the value
  is \code{NA}, \code{NaN} or infinite.
  
  Some functions, such as \code{\link{optim}} with method
  \code{L-BFGS-B} (and therefore \code{\link{find.mle}}), require that
  every value is finite.  Optimisation with these functions also
  requires that the target functions to not generate errors.
  \code{protect} catches these issues, returning the value of
  \code{fail.value} instead.

  No check is made that \code{f} returns a single value, but it should.
}

\examples{
f <- function(x) log(x)
g <- protect(f)
f(0) # -Inf
g(0, fail.value=-999) # -999

f <- function(x) {
  if ( x < 1 )
    stop("dummmy error")
  x
}
g <- protect(f)
\dontrun{
f(0) # error
}
g(0, fail.value=-999) # -999
}

\keyword{programming}
\keyword{internal}
\author{Richard G. FitzJohn}
