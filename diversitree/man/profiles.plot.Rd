\name{profiles.plot}
\alias{profiles.plot}
\title{Plot Marginal Distributions from MCMC}
\description{Simple plotting assistance for plotting output from MCMC
  runs}

\usage{
profiles.plot(y, col.line, col.fill, xlim=NULL, ymax=NULL, n.br=50,
              opacity=.5, xlab="Parameter estimate",
              ylab="Probability density", legend.pos=NULL, ...)
}

\arguments{
  \item{y}{Data frame, columns of which will be plotted as separate
    profiles.}

  \item{col.line}{Vector of colours for the lines.}

  \item{col.fill}{Vector of colours for the fill of the 95\% most
    probable region of the distribution.  If ommited, this will be a
    semi-transparent version of \code{col.line}.}

  \item{xlim}{X-axis limits - calculated automatically if omitted.}

  \item{ymax}{Y-axis upper limit - calculated automatically if omitted.}

  \item{n.br}{Number of breaks along the range of the data.}

  \item{opacity}{Opacity of the filled region (0 is transparent, 1 is
    fully opaque).}

  \item{xlab,ylab}{Axis labels for the plot.}

  \item{legend.pos}{String to pass to \code{\link{legend}} to position
    the legend (for automatic legend building based on the names of
    \code{y}).}

  \item{...}{Additional arguments passed through to \code{\link{plot}}.}
}

\examples{
## For usage, see the example in ?make.bd
}

\keyword{hplot}
\author{Richard G. FitzJohn}
