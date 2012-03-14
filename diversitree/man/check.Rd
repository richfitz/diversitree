\name{check}
\alias{check.fftC}
\alias{check.cvodes}

\title{Check Capabilities of the Diversitree Install}

\description{These check to see if FFTW or CVODES support was included
  in diversitree.  They rarely need to be called directly.}

\usage{
check.fftC(error=TRUE)
check.cvodes(error=TRUE)
}

\arguments{
  \item{error}{Logical: causes an error if the FFTW or CVODES are not
    available if \code{TRUE}}
}

\author{Richard G. FitzJohn}
\keyword{utilities}
