\name{profiles.plot}
\alias{profiles.plot}
\title{Plot Marginal Distributions from MCMC}
\description{Simple plotting assistance for plotting output from MCMC
  runs}

\usage{
profiles.plot(y, col.line, col.fill, xlim=NULL, n.br=50, ...)
}

\arguments{
  \item{y}{Data frame, columns of which will be plotted as separate
    profiles.}
  \item{col.line}{Vector of colours for the lines.}
  \item{col.fill}{Vector of colours for the fill of the 95\% most
    probable region of the distribution.  If ommited, this will be a
    semi-transparent version of \code{col.line}.}
  \item{xlim}{X-axis limits - calculated automatically if omitted.}
  \item{n.br}{Number of breaks along the range of the data.}
  \item{...}{Additional arguments passed through to \code{\link{plot}}.}
}

\usage{
## For usage, see the example in \code{\link{make.bd}}
}

\keyword{hplot}
\author{Richard G. FitzJohn}
