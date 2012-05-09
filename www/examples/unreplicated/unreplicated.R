## % Unreplicated correlated evolution

## `r opts_chunk$set(tidy=FALSE)`
##+ echo=FALSE,results=FALSE
options(show.signif.stars=FALSE, continue="  ")

## First, load diversitree and simulate a pure birth tree (with rate
## .1)
library(diversitree)
set.seed(1)
n.spp <- 30
phy <- ladderize(tree.yule(.1, max.taxa=n.spp))

##+ fig.height=5,fig.cap="Simulated pure birth tree"
plot(phy, show.node.label=TRUE, no.margin=TRUE, cex=.8)

## nd14 will be our target node.
target <- "nd14"

## Convert this name to Ape's node index:
idx <- match(target, phy$node.label) + n.spp

## Identify descendants of this node
desc <- diversitree:::descendants(idx, phy$edge)
desc <- sort(desc[desc < n.spp])

## Here is our focal group, highlighted in red
##+ fig.height=5,fig.cap="Focal group at node 'nd14'"
plot(phy, tip.color=ifelse(1:n.spp %in% desc, "#ef2929", "black"),
     no.margin=TRUE)
nodelabels(node=idx, pch=19, cex=2, col="#ef2929")

## Trait 1 is shared across this group, and with nobody else.
t1 <- as.integer(1:n.spp %in% desc)
names(t1) <- phy$tip.label

## Determine rate coefficients for the evolution of this trait,
## assuming a Mk2 model:
lik1 <- make.mk2(phy, t1)
fit1 <- find.mle(lik1, rep(.1, 2))
zapsmall(coef(fit1))

## Now, suppose we have two perfectly codistribted traits "a" and "b";
traits <- data.frame(a=t1, b=t1)

## Making a 2-trait likelihood function, assuming that the two traits
## are independent:
lik2 <- make.mkn.multitrait(phy, traits, depth=0)

## And fitting the ML model
p2 <- rep(coef(fit1), 2)
names(p2) <- NULL
fit2 <- find.mle(lik2, p2)

## We get the same basic coefficients when the traits are assumed to
## evolve independently.
zapsmall(coef(fit1))
zapsmall(coef(fit2))
all.equal(coef(fit2)[1:2], coef(fit1), check.attr=FALSE)

## Now, allow correlated evolution by allowing the $q_{01}$ and
## $q_{10}$ parameters to vary according to the state of the other
## trait:
lik3 <- make.mkn.multitrait(phy, traits, depth=1)
p3 <- rep(0, 8)
names(p3) <- argnames(lik3)
p3[argnames(lik2)] <- coef(fit2)

## Fit the model
##+ fit.lik3,cache=TRUE
fit3 <- find.mle(lik3, rep(.1, 8))

## Basically once trait 'a' has shifted from 0->1, we instantaneously
## shift in state 'b'.
zapsmall(coef(fit3))

## This model is a statistical improvement over the uncorrelated
## model.
anova(fit3, uncorrelated=fit2)

## Likelihood improvement:
dll <- fit3$lnLik - fit2$lnLik
dll

## Simulate 100 traits on the tree with an equal forward/backward
## transition rate of `r`
r <- .02
t2 <- replicate(100,
                sim.character(phy, c(r, r), x=0, model="mk2"),
                simplify=FALSE)

## And fit these the same way:
## This fails to run in sowsear with blank lines below -- fix and push.
f <- function(t2) {
  lik0 <- make.mkn.multitrait(phy, data.frame(a=t1, b=t2), depth=0)
  p0 <- rep(r, 4)
  fit0 <- find.mle(lik0, p0)
  lik1 <- make.mkn.multitrait(phy, data.frame(a=t1, b=t2))
  p1 <- rep(0, 8)
  names(p1) <- argnames(lik1)
  p1[argnames(lik0)] <- coef(fit0)
  fit1 <- find.mle(lik1, p1)
  fit1$lnLik - fit0$lnLik
}

## This takes a minute or so, so run it in parallel:
library(multicore)
##+ fit.random,cache=TRUE
ans <- unlist(mclapply(t2, f))

## Histogram of likelihood improvements, with previous case in red,
## and 5% cut-off in blue:
##+ fig.cap=""
par(mar=c(4.1, 4.1, .5, .5))
hist(ans, xlim=range(ans, dll), main="", las=1,
     xlab="Likelihood impovement")
abline(v=dll, col="red")
abline(v=qchisq(1/20, 4, lower=FALSE)/2, col="blue")

## # BiSSE

## With the tree above, fit a BiSSE model
##+ bisse.ssd,cache=TRUE
lik <- make.bisse(phy, t1, control=list(method="CVODES"))
p <- starting.point.bisse(phy)
fit.ssd <- find.mle(lik, p)

## Fit the state independent model (constrain both speciation and
## extinction rates to not vary with the trait).
##+ bisse.sid,cache=TRUE
lik.sid <- constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0)
fit.sid <- find.mle(lik.sid, p[argnames(lik.sid)])

## There is no significant association between the trait and
## diversification here:
anova(fit.ssd, sid=fit.sid)

## Next, tweak the tree so that the focal clade is apparently
## different.  This is a bit of a pain to do.  First, identify all
## edges descended from that node
desc <- diversitree:::descendants(idx, phy$edge)
i <- phy$edge[,2] %in% desc

##+ fig.height=5,fig.cap="Edges descended from focal node"
plot(phy, no.margin=TRUE,
     edge.color=ifelse(i, "#ef2929", "black"),
     tip.color=ifelse(1:n.spp %in% desc, "#ef2929", "black"))     

## Compress these red edges towards the present, and lengthen the stem
## lineage accordingly.  Scale branches by a factor 'p' (must be less
## than 1)
p <- .3

j <- match(idx, phy$edge[,2])
phy2 <- phy
phy2$edge.length[i] <- phy2$edge.length[i] * p
phy2$edge.length[j] <-
  phy$edge.length[j] + branching.times(phy)[target] * (1 - p)

## The modified tree is still ultrametric
is.ultrametric(phy2)

## Here is the compressed tree:
##+ fig.height=5,fig.cap="Edges descended from focal node"
plot(phy2, no.margin=TRUE,
     edge.color=ifelse(i, "#ef2929", "black"),
     tip.color=ifelse(1:n.spp %in% desc, "#ef2929", "black"))     

## Fit BiSSE models to the tree, first with state-dependent
## diversification
##+ bisse2.sdd,cache=TRUE
lik2 <- make.bisse(phy2, t1, control=list(method="CVODES"))
p <- starting.point.bisse(phy2)
fit2.ssd <- find.mle(lik2, p)

## And without
##+ bisse2.sid,cache=TRUE
lik2.sid <- constrain(lik2, lambda1 ~ lambda0, mu1 ~ mu0)
fit2.sid <- find.mle(lik2.sid, p[argnames(lik2.sid)])

## There is now significant association between the trait and
## diversification!
anova(fit2.ssd, sid=fit2.sid)

## Generalise this a bit to look how the likelihood improvement varies
## with the scaling factor.
g <- function(p) {
  phy2 <- phy
  phy2$edge.length[i] <- phy2$edge.length[i] * p
  phy2$edge.length[j] <-
    phy$edge.length[j] + branching.times(phy)[target] * (1 - p)

  lik2 <- make.bisse(phy2, t1, control=list(method="CVODES"))
  p <- starting.point.bisse(phy2)
  fit2.ssd <- find.mle(lik2, p)

  lik2.sid <- constrain(lik2, lambda1 ~ lambda0, mu1 ~ mu0)
  fit2.sid <- find.mle(lik2.sid, p[argnames(lik2.sid)])
  fit2.ssd$lnLik - fit2.sid$lnLik
}

## Run this with scaling parameters varying from .1 (extremely
## compressed towards the future) and 1 (unchanged from the original
## tree).
##+ bisse.scale,cache=TRUE
p <- seq(1, .1, length=16)
ans.b <- unlist(mclapply(p, g))

## The improvement in the likelihood increases rapidly as the tree
## deviates more and more from the BD model:
##+ fig.cap="Likelihood improvement vs scaling factor, with critical value at the 5% level indicated with a blue line." 
par(mar=c(4.1, 4.1, .5, .5))
plot(ans.b ~ p, las=1, xlab="Scaling factor (1 is unchanged)",
     ylab="Likelihood improvement", type="o", pch=19)
abline(h=qchisq(1/20, 2, lower=FALSE)/2, col="blue")

## However, we can partition the tree and allow rates to vary here.
## First check that this model will fit better without accounting for
## the state-dependent diversification:
lik.bd <- make.bd(phy2)
fit.bd <- find.mle(lik.bd, starting.point.bd(phy2))

## High speciation and extinction rates, compared with the values
## used to simulate the tree (.1):
coef(fit.bd)

lik.bd2 <- make.bd.split(phy2, target)
fit.bd2 <- find.mle(lik.bd2, c(.1, 0, .1, 0), method="subplex")

## The base partition (black branches in the figures) have rates
## similar to the simulated values.  The focal clade (red in the
## figures) has strongly elevated rates.
coef(fit.bd2)

## This improvement is significant;
anova(fit.bd2, unsplit=fit.bd)

## Split BiSSE likelihood function
lik.bs <- make.bisse.split(phy2, t1, target)

## Constrain trait evolution over the whole tree:
lik.bs.sdd <- constrain(lik.bs, q01.2 ~ q01.1, q10.2 ~ q10.1,
                        mu1.1 ~ mu0.1, mu1.2 ~ mu0.2)
## (I also disallowed state-dependent extinction as that causes
## problems trying to explain the pattern of trait evolution)

## And make a version with no state-dependent diversification:
lik.bs.sid <- constrain(lik.bs.sdd,
                        lambda1.1 ~ lambda0.1,
                        lambda1.2 ~ lambda0.2)

## The state-independent version has 6 parameters; two speciation
## rates (one per partition), two extinction rates, and two character
## transition rates.
argnames(lik.bs.sid)

## Starting point, based on the birth-death fit:
p.sid <- rep(0, 6)
names(p.sid) <- argnames(lik.bs.sid)
p.sid[sub("\\.", "0.", names(coef(fit.bd2)))] <- coef(fit.bd2)
p.sid[c("q01.1", "q10.1")] <- .1

##+ bisse.split.sid, cache=TRUE
fit.bs.sid <- find.mle(lik.bs.sid, p.sid)

## The coefficients here look like the birth-death case, with low
## rates of $0\to 1$ transition, and no reverse ($1\to 0$)
## transitions.
zapsmall(coef(fit.bs.sid))

## Expand this point to allow for state-dependent diversification:
p.sdd <- coef(fit.bs.sid, TRUE)[argnames(lik.bs.sdd)]

##+ bisse.split.sdd, cache=TRUE
fit.bs.sdd <- find.mle(lik.bs.sdd, p.sdd)

## This is now not significant (though the improvement in likelihood
## is greater than I would have thought).
anova(fit.bs.sdd, sid=fit.bs.sid)

## I think what is going on here is that currently the stem leading to
## the focal group is included in that group.  However, the (fake)
## diversification shift happens rather closer to the node.  So that
## aspect of the tree still needs explaining.

## ## Document details

## Document compiled on `r as.character(Sys.time(), usetz=TRUE)`
## With R version `r R.version$string`
## and diversitree version `r packageVersion("diversitree")`.
## Original source: [unreplicated.R](unreplicated.R)
