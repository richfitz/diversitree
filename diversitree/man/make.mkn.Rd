\name{make.mkn}
\alias{make.mkn}
\alias{make.mk2}

\title{Mk2 and Mk-n Models of character evolution}

\description{Prepare to run a Mk2/Mk-n model on a phylogenetic tree and
  binary/discrete trait data.  This fits the Pagel 1994 model,
  duplicating the \code{ace} function in ape.  Differences with that
  function include (1) alternative root treatments are possible, (2)
  easier to tweak parameter combinations through
  \code{\link{constrain}}, and (3) run both MCMC and MLE fits to
  parameters.  Rather than exponentiate the Q matrix, this
  implementation solves the ODEs that this matrix defines.  This may or
  may not be robust on trees leading to low probabilities.
}

\usage{
make.mk2(tree, states)
make.mkn(tree, states, k)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{states}{A vector of character states, each of which must be 0 or
    1 for \code{make.mk2} or 1 to \code{k} for \code{make.mkn}.}
  \item{k}{Number of states to model.}
}

\details{
  \code{make.mk2} and \code{make.mkn} return functions of class \code{mkn}.
  These functions have argument list (and default values)
  \preformatted{
    f(pars, pars, prior=NULL, root=ROOT.OBS, root.p=NULL, fail.value=NULL)
  }
  This may change in future versions.
  The arguments are interpreted as
  \itemize{
    \item \code{pars} For \code{make.mk2}, a vector of two parameters,
    in the order \code{q01}, \code{q10}.  For \code{make.mkn}, a
    vector of \code{k(k-1)} parameters, in the order
    \code{q12,q13,...q1k, q21,q23,...,q2k,...qk(k-1)}, corresponding
    to the off-diagonal elements of the \code{Q} matrix in row order.
    \item \code{prior} A vector of rates for an exponential prior.
    The prior probability distribution has probability density
    \deqn{r_i e^{-r_i x_i}}{r_i*exp(-r_i*x_i)}
    where the \eqn{i} denotes the \eqn{i}th parameter.  There is no
    ability to specify a form of prior other than an exponential
    function here (yet), and this part of the inteface is subject to
    change!
    \item \code{root}: Behaviour at the root (see Maddison et al. 2007,
    FitzJohn et al. 2009).  The possible options are
    \itemize{
      \item \code{diversitree:::ROOT.FLAT}: A flat prior, weighting
      \eqn{D_0}{D0} and \eqn{D_1}{D1} equally.
      \item \code{diversitree:::ROOT.EQUI}: Use the equilibrium distribution
      of the model, as described in Maddison et al. (2007).
      \item \code{diversitree:::ROOT.OBS}: Weight \eqn{D_0}{D0} and
      \eqn{D_1}{D1} by their relative probability of observing the
      data, following FitzJohn et al. 2009:
      \deqn{D = D_0*\frac{D_0}{D_0 + D_1} + D_1*\frac{D_1}{D_0 + D_1}}{
 	D = D0 * D0/(D0 + D1) + D1 * D1/(D0 + D1)}
      \item \code{diversitree:::ROOT.GIVEN}: Root will be in state i
      with probability \code{root.p[i]}.
      \item \code{diversitree:::ROOT.BOTH}: Don't do anything at the root,
      and return both values.  (Note that this will not give you a
      likelihood!).
    }
    \item \code{root.p}{Vector of probabilities/weights to use when
      \code{ROOT.GIVEN} is specified.  Must be of length \code{k}.}
    \item \code{fail.value} Optional value to use if the likelihood
    calculation fails.
  }
}

\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration,
  and \code{\link{make.bisse}} for state-dependent birth-death models.
}

\examples{
## Simulate a tree and character distribution and look at the marginal
## likelihood of the transition parameters.
pars <- c(.1, .1, .03, .03, .1, .2)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=25, max.t=Inf, x0=0)[[1]]
p <- c(.1, .1) # initial parameter guess

## Using Mk2
lik.mk2 <- make.mk2(phy, phy$tip.state)
fit.mk2 <- find.mle(lik.mk2, p)
fit.mk2[1:2]

## Using the more general Mk-n (ignore warning)
lik.mkn <- make.mkn(phy, phy$tip.state + 1, 2)
fit.mkn <- find.mle(lik.mkn, p)
fit.mkn[1:2]

## These are the same to a reasonable degree of accuracy
all.equal(fit.mk2[1:2], fit.mkn[1:2], tolerance=1e-7)

## Equivalence to Ape's ace function:
model <- matrix(c(0, 2, 1, 0), 2)
fit.ape <- ace(phy$tip.state, phy, "discrete", model=model, ip=p)

## Rerun the diversitree version with the same root conditions as ape;
## these are the same to a reasonable degree of accuracy, too (the
## matrix exponentiation is slightly less accurate than the ODE
## solving approach.  The make.mk2 version is exact)
fit.mk2 <- find.mle(lik.mk2, p, root=diversitree:::ROOT.GIVEN,
                    root.p=c(1,1))
fit.ape[c("rates", "loglik")]
fit.mk2[1:2]

all.equal(fit.ape[c("rates", "loglik")], fit.mk2[1:2],
          check.attributes=FALSE, tolerance=4e-5)
}

\author{Richard G. FitzJohn}
\keyword{models}
