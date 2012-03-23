\name{set.defaults}
\alias{set.defaults}

\title{Set Default Arguments of a Function}

\description{
  Set the default values of formal arguments of a function.
}

\usage{
set.defaults(f, ..., defaults)
}

\arguments{
  \item{f}{A function}

  \item{...}{Named arguments to be set}

  \item{defaults}{A named list of arguments}
}

\details{The repetitive argument lists of many of diversitree's
  likelihood functions are the motivation for this function.

  For example, the likelihood function that \code{make.bisse} produces
  takes arguments \code{condition.surv}, \code{root}, and \code{root.p},
  each with default values.  If you dislike the defaults, you can change
  them by passing in alternative values when computing likelihoods, or
  when doing an ML search.  However, this can get tedious if you are
  using a function a lot, and your code will get cluttered with lots of
  statements like \code{condition.surv=FALSE}, some of which you may
  forget.  See the example below for how to avoid this.
}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)
lik <- make.bisse(phy, phy$tip.state)

## default arguments:
args(lik)

lik.no.cond <- set.defaults(lik, condition.surv=FALSE)
args(lik.no.cond)

## Multiple arguments at once:
lik2 <- set.defaults(lik, root=ROOT.GIVEN, root.p=c(0, 1))
args(lik2)

## Equivalently (using alist, not list -- see ?alist)
defaults <- alist(root=ROOT.GIVEN, root.p=c(0, 1))
lik3 <- set.defaults(lik, defaults=defaults)
identical(lik2, lik3)
}

\author{Richard G. FitzJohn}
\keyword{programming}
