\name{asr-mkn}
\alias{asr.mkn}
\alias{asr.marginal.mk2}
\alias{asr.marginal.mkn}
\alias{asr.joint.mk2}
\alias{asr.joint.mkn}
\alias{asr.stoch.mk2}
\alias{asr.stoch.mkn}

\title{Ancestral State Reconstruction Under Mk2/Mkn}

\description{Perform ancestral state reconstruction under Mk2 and other
  constant rate Markov models.  Marginal, joint, and stochastic
  reconstructions are supported.  Documentation is still in an early
  stage, and mostly in terms of examples.}

\usage{
\method{asr.marginal}{mkn}(lik, pars, nodes=NULL, ...)
\method{asr.joint}{mkn}(lik, pars, n=1, root.state=NA, simplify=TRUE,
       intermediates=FALSE, ...)
\method{asr.stoch}{mkn}(lik, pars, n=1, root.state=NA, ...)
}

\arguments{
  \item{lik}{A likelihood function, returned by \code{make.mk2} or
    \code{make.mkn}.}
  \item{pars}{A vector of parameters, suitable for \code{lik}.}
  \item{nodes}{For \code{asr.marginal} only; an optional vector of
    nodes to return ancestral states for (using ape's index).  By
    default, all nodes are returned.}
  \item{n}{The number of samples to draw from the joint distribution, or
    number of stochastic reconstructions to make.}
  \item{root.state}{For \code{asr.marginal} and \code{asr.joint}; the
    state to start the root in.  If \code{NA} is given (the default),
    the root state is drawn from the conditional likelihoods at the
    root.}
  \item{simplify}{Should the states be simplified into a matrix or
    left unsimplified as a list?}
  \item{intermediates}{Logical: Return intermediates as an attribute?
    Used internally, and probably not of general interest.}
  \item{...}{Additional arguments passed through to future methods.
    Currently unused.}
}

\details{
  Output will differ slightly when \code{mk2} and \code{mkn} models are
  used as \code{lik}, as \code{mk2} uses states 0/1, while 2-state
  \code{mkn} uses 1/2.

  This is all quite slow.  Faster versions are coming eventually.
}

\examples{
## Start with a simple tree evolved under a constant rates birth-death
## model with asymetric character evolution
pars <- c(.1, .1, .03, .03, .03, .06)
set.seed(1)
phy <- trees(pars, "bisse", max.taxa=50, max.t=Inf, x0=0)[[1]]

## Here is the true history.  The root node appears to be state 1 (red)
## at the root, despite specifying a root of state 0 (x0=0, in statement
## above).  This is because the tree started with a single lineage, but
## had changed state by the time the first speciation event happened.
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, main="True history")

## All of the methods need a likelihood function; build a mk2 function:
lik <- make.mk2(phy, phy$tip.state)

## Using the true parameters, compute the marginal ancestral state
## reconstructions:
st.m <- asr.marginal(lik, pars[5:6])

## There is still not a good stand-along plotting command for nodes.
## For now, use ape's nodelabels().
plot(h, phy, main="Marginal ASR", show.node.state=FALSE)
nodelabels(thermo=t(st.m), piecol=1:2, cex=.5)

## Again, with the true parameters, a sample from the joint
## distribution:
st.j <- asr.joint(lik, pars[5:6])

## Plotting this sample against the true values.
plot(h, phy, main="Joint ASR", show.node.state=FALSE)
nodelabels(pch=19, col=st.j + 1)

## This is just one sample, and is not very accurate in this case!  Make
## 100 such samples and average them:
st.j2 <- asr.joint(lik, pars[5:6], 100)
st.j2.mean <- colMeans(st.j2)

## Plotted, this looks much like the joint sample:
plot(h, phy, main="Joint ASR (averaged)", show.node.state=FALSE)
nodelabels(thermo=1-st.j2.mean, piecol=1:2, cex=.5)

## They are similar, but in this case, the joint method generates
## samples with more state '1' nodes than the marginal method.  Other
## trees would show the opposite pattern.
plot(st.m[2,], st.j2.mean, xlab="Marginal", ylab="Joint", las=1)
abline(0, 1)

## Finally, the stochastic character mapping.  This uses samples from
## the joint distribution at its core.
st.s <- asr.stoch(lik, pars[5:6])
plot(st.s, phy)

## Again, multiple samples can be done at once.  There is a function for
## summarising histories, but it is still in the works.

## Repeating the above with a two-state mkn model:
lik2 <- make.mkn(phy, phy$tip.state + 1, 2, FALSE)

## Everything works:
st2.m <- asr.marginal(lik2, pars[5:6])
st2.j <- asr.joint(lik2, pars[5:6], 100)
st2.s <- asr.stoch(lik2, pars[5:6])

## Marginal likelihoods agree:
all.equal(st.m, st2.m)
## Joint reconstructions are stochastic, so just check with a
## regression:
summary(lm(colMeans(st2.j) - 1 ~ colMeans(st.j2) - 1))

## Integrate parameter uncertainty, and see how far down the tree there
## is any real information on parameter states for this tree:
set.seed(1)
prior <- make.prior.exponential(.5)
samples <- mcmc(lik, pars[5:6], 1000, c(1, 1), 0, prior=prior,
                print.every=100)
st.m.avg <- rowMeans(apply(samples[2:3], 1, asr.joint, lik=lik))

plot(h, phy, main="MCMC Averaged ASR", show.node.state=FALSE)
nodelabels(thermo=1 - st.m.avg, piecol=1:2, cex=.5)
}

\author{Richard G. FitzJohn}
\keyword{model}
