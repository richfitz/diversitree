\name{make.musse.multitrait}
\alias{make.musse.multitrait}
\alias{make.mkn.multitrait}
\alias{musse.multitrait.translate}
\alias{mkn.multitrait.translate}
\alias{starting.point.musse.multitrait}

\title{MuSSE: Multi-State Speciation and Extinction (Multiple Binary Traits
  Version)}

\description{Prepare to run MuSSE or Mkn (Multi-State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.
  
  This is a helper function that wraps the basic MuSSE/Mkn models for
  the case of a combination of several binary traits; its
  parametrisation and argument handling are a little different to the
  other models in diversitree.
}
  
\usage{
make.musse.multitrait(tree, states, sampling.f=NULL,
                      depth=NULL, allow.multistep=FALSE,
                      strict=TRUE, control=list())
make.mkn.multitrait(tree, states,
                    depth=NULL, allow.multistep=FALSE,
                    strict=TRUE, control=list())

musse.multitrait.translate(n.trait, depth=NULL, names=NULL,
                           allow.multistep=FALSE)
mkn.multitrait.translate(n.trait, depth=NULL, names=NULL,
                         allow.multistep=FALSE)

starting.point.musse.multitrait(tree, lik, q.div=5, yule=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}

  \item{states}{A \code{data.frame} of character states, each column of
    which represents a different binary state (with values 0 or 1), and
    each row of which represents a taxon.  The row names of
    \code{states} must be the names that correspond to the tip labels in
    the phylogenetic tree (\code{tree$tip.label}).  The column names
    must be unique and a single character long.  The character "0"
    (zero) is reserved and may not be used. \code{NA} values are allowed
    in one or more columns when one or more traits is unknown for a
    taxon.}

  \item{depth}{A scalar or vector of length 3 indicating the depth of
    interactions to include in the model.  See Details.}

  \item{allow.multistep}{Should transition rates be included that imply
    simultaneous changes in more than one trait?  By default this is not
    allowed, but if set to \code{TRUE} these rates are included at the
    end of the parameter vector.  Warning: treatment of these will
    change in future versions!}
  
  \item{sampling.f}{Scalar with the estimated proportion of extant
    species that are included in the phylogeny.  A value of \code{0.75}
    means that three quarters of extant species are included in the
    phylogeny.  By default all species are assumed to be known.  In the
    future, this will expand to allow state-specific sampling rates.}
  
  \item{strict}{Each column in \code{states} is always checked to make
    sure that the values are 0 or 1.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present.  The likelihood models tend to be poorly behaved
    where states are missing, but there are cases (missing intermediate
    states for meristic characters) where allowing such models may be
    useful.  Note that this model may misbehave even if this check is
    met, due to combinations of traits being absent.}

  \item{control}{List of control parameters for the ODE solver.  See
    details in \code{\link{make.bisse}}.}

  \item{lik}{A likelihood function created by
    \code{make.musse.multitrait}.}

  \item{q.div}{Ratio of diversification rate to character change rate.
    Eventually this will be changed to allow for Mk2 to be used for
    estimating q parameters.}
  
  \item{yule}{Logical: should starting parameters be Yule estimates
    rather than birth-death estimates?}

  \item{n.trait}{Number of binary traits.}

  \item{names}{Vector of names for the traits when using
    musse.multitrait.translate (optional).}
}

\details{
  Suppose that you have two binary traits that may affect speciation and
  extinction.  In previous versions of diversitree, you had to code the
  possible combinations as states 1, 2, 3, 4, which makes the
  interpretation of the speciation rates (\code{lambda1},
  \code{lambda2}, etc) unintuitive.

  Let \code{states} is a data.frame with columns "A" and "B",
  representing the two binary traits.  We can write the speciation rate
  as
  %
  \deqn{\lambda_0 + \lambda_A X_A + \lambda_B X_B + \lambda_{AB}X_AX_B}{
    lambda_0 + lambda_A * X_A + lambda_B * X_B + lambda_AB * X_A * X_B}
  %
  where \eqn{X_A} and \eqn{X_B} are indicator variables that take the
  value of trait A and B respectively (with values 0 or 1).  In this
  form, \eqn{\lambda_0}{lambda_0} is the intercept,
  \eqn{\lambda_A}{lambda_A} and \eqn{\lambda_B}{lambda_B} are "main
  effects" of traits A and B, and \eqn{\lambda_{AB}}{lambda_AB} is the
  "interaction" between these.  We can do a similar trick for the
  extinction rates.

  For character transition rates, we first consider changes only in a
  single trait.  For our two trait case we have four "types" of
  character change allowed (A 0->1, A 1->0, B 0->1, and B 1->0), but the
  rates of change for trait A might depend on the current state of trait
  B (and vice versa).  So we have, for the A0->1 trait change
  \eqn{q_{A01,0} + q_{A01,B} \times X_B}{qA01.0 + qA01.B * X_B}.  Note that one fewer levels of
  interaction are possible for these character changes than for the
  speciation/extinction parameters.

  It may sometimes be desirable to have the multi-trait changes in the
  model.  At present, if \code{allow.multistep} is \code{TRUE}, all the
  multiple change transitions are included at the end of the parameter
  vector.  For the two trait case these are labelled \code{q00.11},
  \code{q10.01}, \code{q01.10}, and \code{q11.00}, where \code{qij.kl}
  represents a change from (A=i, B=j) to (C=k, D=l).  The argument name,
  and treatment, of these may change in future.

  This approach generalises out to more than two traits.  For \code{N}
  traits, interactions are possible up to the \code{N}th order for
  lambda and mu, and up to the \code{N-1}th order for q.  The
  \code{depth} argument controls how many of these are returned.  If
  this is a scalar, then the same level is used for \code{lambda},
  \code{mu} and \code{q}.  If it is a vector of length 3, then different
  depths are used for these three types of parameters.  By default, all
  possible interactions are returned and the model has the same number
  of degrees of freedom as the models returned by \code{make.musse}
  (except for a reduction in the possible q parameters when
  \code{allow.multistep} is \code{FALSE}).  Parameters can then be
  further refined with \code{constrain}.
}

\seealso{\code{\link{make.bisse}} for the basic binary model, and
  \code{\link{make.musse}} for the basic multistate model.}

\examples{
## The translation between these two bases is fairly straightforward; if
## we have a vector of parameters in our new basis 'p' we can convert it
## into the original MuSSE basis ('q') through this matrix:
tr <- musse.multitrait.translate(2)
tr

## Notice that the rows that correspond to transitions in multiple
## traits are all zero by default; this means that these q values will
## be zero regardless of the parameter vector used.
tr["q00.11",]

## And here is the section of the transition matrix corresponding to the
## lambda values; every rate gets a contribution from the intercept term
## (lambda0), lambda10 and lambda11 get a contribution from lambdaA, etc.
tr[1:4,1:4]

## There is currently no nice simulation support for this, so bear with
## an ugly script to generate the tree and traits.
pars <- c(.10, .15, .20, .25, # lambda 00, 10, 01, 11
          .03, .03, .03, .03, # mu 00, 10, 01, 11
          .05, .05, .0,       # q00.10, q00.01, q00.11
          .05, .0,  .05,      # q10.00, q10.01, q10.11
          .05, .0,  .05,      # q01.00, q01.10, q01.11
          .0,  .05, .05)      # q11.00, q11.10, q11.01
set.seed(2)
phy <- tree.musse(pars, 60, x0=1)

states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
rownames(states) <- phy$tip.label

## Here, states has row names corresponding to the different taxa, and
## the states of two traits "A" and "B" are recorded in the columns.
head(states)

## Note that transition from the original MuSSE basis to this basis is
## only possible in general when depth=n.trait and allow.multistep=TRUE
## (as only this generates a square matrix that is invertible).
## However, when it is possible to express the set of parameters in the
## new basis (as it is above), this can be done through a pseudoinverse
## (here, a left inverse).
pars2 <- drop(solve(t(tr) \%*\% tr) \%*\% t(tr) \%*\% pars)

## Going from our new basis to the original MuSSE parameters is always
## straightforward.  This is done automatically in the likelihood
## function.
all.equal(tr \%*\% pars2, pars, check.attributes=FALSE)

## This shows that the two traits act additively on speciation rate
## (lambdaAB is zero), that there is no effect of any trait on
## extinction (the only nonzero mu parameter is mu0) and transition
## rates for one trait are unaffected by other traits (the only nonzero
## q parameters are the qXij.0 parameters; qXij.Y parameters are all
## zero).

## Here is our new MuSSE function parametrised as a multi-trait
## function:
lik <- make.musse.multitrait(phy, states)

## Here are the argument names for the likelihood function.
argnames(lik)

## Basic MuSSE function for comparison
lik.m <- make.musse(phy, phy$tip.state, 4)
argnames(lik.m)

## Rather than fit this complicated model first, let's start with a
## simple model with no state dependent diversification.  This model
## allows the forwards and backwards transition rates to vary, but the
## speciation and extinction rates do not depend on the character
## state:
lik0 <- make.musse.multitrait(phy, states, depth=0)
argnames(lik0)

p <- starting.point.musse.multitrait(phy, lik0)
fit0 <- find.mle(lik0, p)

## Now, allow the speciation rates to vary additively with both
## character states (extinction and character changes are left as in the
## previous model)
lik1 <- make.musse.multitrait(phy, states, depth=c(1, 0, 0))

## Start from the previous ML point:
p <- starting.point.musse.multitrait(phy, lik1)
p[names(coef(fit0))] <- coef(fit0)

fit1 <- find.mle(lik1, p)

## The likelihood improves, but the difference is not statistically
## significant (p = 0.35).
anova(fit1, fit0)

## We can fit an interaction for the speciation rates, too:
lik2 <- make.musse.multitrait(phy, states, depth=c(2, 0, 0))
p <- starting.point.musse.multitrait(phy, lik2)
p[names(coef(fit1))] <- coef(fit1)
fit2 <- find.mle(lik2, p)

## There is next to no support for the interaction term (which is good,
## as the original model did not have any interaction!)
anova(fit2, fit1)

## Constraining also works with these models.  For example, constraining
## the lambdaA parameter to zero:
lik1b <- constrain(lik1, lambdaA ~ 0)
argnames(lik1b)
p <- starting.point.musse.multitrait(phy, lik1b)
p[names(coef(fit0))] <- coef(fit0)
fit1b <- find.mle(lik1b, p)
anova(fit1b, fit0)

## Or constraining both main effects to take the same value:
lik1c <- constrain(lik1, lambdaB ~ lambdaA)
argnames(lik1c)
p <- starting.point.musse.multitrait(phy, lik1c)
p[names(coef(fit0))] <- coef(fit0)
fit1c <- find.mle(lik1c, p)
anova(fit1c, fit0)
}

\author{Richard G. FitzJohn}
\keyword{models}
