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
  Some functions, such as \code{\link{optim}} with method
  \code{L-BFGS-B} (and therefore \code{\link{find.mle}}), require that
  every value is finite.  Optimisation with these functions also
  requires that the target functions to not generate errors.
  \code{protect} catches these issues, returning the value of
  \code{fail.value} instead.

  When using these, some care might need to be taken to specify a
  value of \code{fail.value} that works well; it should be small enough
  that it is considered a \dQuote{low} probability, but large enough
  that it does not lead to undue bumpyness of the likelihood surface.

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
