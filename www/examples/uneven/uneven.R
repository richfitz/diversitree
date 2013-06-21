## % Dealing with uneven taxonomic sampling under BiSSE

### Set options
##+ echo=FALSE,results=FALSE
opts_chunk$set(tidy=FALSE, fig.height=5)
options(show.signif.stars=FALSE)

## Not all trees have taxonomically even sampling, and this violates
## an assumption of BiSSE (and other methods) when trying to compute
## likelihoods.  However, we often know about this unevenness; we know
## that there are really X species in one clade, and Y in another.  It
## is possible to use the "split" type models to account for this
## unevenness.

library(diversitree)

## First, simulate a phylogeny with 100 species
set.seed(1)
phy <- ladderize(tree.bisse(c(.1, .1, 0, 0, .03, .03), max.taxa=100))

## This tree breaks into three major clades, subtended by these nodes:
nodes <- c("nd3", "nd4", "nd5")

##+ fig.cap="Simulated tree"
plot(phy, type="fan", no.margin=TRUE)
i <- match(nodes, phy$node.label) + length(phy$tip.label)
nodelabels(node=i, pch=19, cex=4)
nodelabels(nodes, node=i, cex=.7, col="white", bg=NA, frame="none")

## Next, determine the descendents of these three clades:
spp <- lapply(nodes, function(n)
              phy$tip.label[get.descendants(n, phy, TRUE)])

## From each clade, sample 10 species to retain.
spp2 <- lapply(spp, sample, 10)
## To preserve the nodes here, it is necessary to add in a couple of
## species.
spp2[[1]] <- c(spp2[[1]], "sp65")
spp2[[2]] <- c(spp2[[2]], "sp1")

## And from this, determine species to drop.
to.drop <- setdiff(phy$tip.label, unlist(spp2))

## Build a new tree (ape's `drop.tip` tends to shuffle node labels, so
## using a hacked version -- I've been meaning to do a bug report
## around this).
phy2 <- diversitree:::drop.tip.fixed(phy, to.drop)

##+ fig.cap="Simulated tree, unevenly sampled."
plot(phy2, type="fan", no.margin=TRUE)
i2 <- match(nodes, phy2$node.label) + length(phy2$tip.label)
nodelabels(node=i2, pch=19, cex=4)
nodelabels(nodes, node=i2, cex=.7, col="white", bg=NA, frame="none")

## The three clades are sampled unevenly; the smallest clade has 79%
## coverage, while the largest has 17% coverage.
mapply(function(x, y) length(y)/length(x), spp, spp2)

states <- phy$tip.state
mapply(function(x, y) as.numeric(table(states[y]) / table(states[x])),
       spp, spp2, SIMPLIFY=FALSE)

## Right, but how do we use BiSSE with this?

## First, using ML fit models with the correct tree:
p <- starting.point.bisse(phy)
lik <- make.bisse(phy, states, control=list(backend="CVODES"))
fit <- find.mle(lik, p)

## ...and with the incorrect tree but completely ignoring the
## sampling ("partial")
states2 <- states[phy2$tip.label]
lik.p <- make.bisse(phy2, states2, control=list(backend="CVODES"))
fit.p <- find.mle(lik.p, p)

## Note that the speciation rate estimates are reduced in the partial
## tree.
rbind(full=coef(fit), partial=coef(fit.p))[,c("lambda0", "lambda1")]

## Now, correcting for the overall level of sampling, but ignoring the
## pattern of the sampling:
sampling.f <- table(states2) / table(states)

lik.s <- make.bisse(phy2, states2, sampling.f=sampling.f,
                     control=list(backend="CVODES"))
fit.s <- find.mle(lik.s, p)

## Now the coefficients are closer to the ones used in the
## simulations:
rbind(full=coef(fit), partial=coef(fit.p),
      sampled=coef(fit.s))[,c("lambda0", "lambda1")]

## Finally, correcting for sampling.  To do this, we can make a split
## function with three nodes.

## From above, this is the sampling fraction within the three major
## clades:
sampling.f.p <- mapply(function(x, y)
                       as.numeric(table(states[y])/table(states[x])),
                       spp, spp2, SIMPLIFY=FALSE)

## to this, we need to add the background sampling fraction.  Here it
## is complete (there are no branches missing below the nodes
## nd3-nd5), so we specify c(1,1) for the first entry:
sampling.f.p <- c(list(c(1,1)), sampling.f.p)

## In the end, we have a list with four elements (for n nodes, there
## will be n+1 entries).  The ith entry corresponds to the clade
## subtended by the i-1th node.  The first entry always corresponds to
## the background group, and includes the root node.
sampling.f.p

## We then can feed this to the `make.bisse.split` function:
lik.p <- make.bisse.split(phy2, states2,
                          sampling.f=sampling.f.p,
                          nodes=nodes)

## However, this requires many parameters, and we don't actually want
## to fit models that have different parameters in different
## partitions:
argnames(lik.p)

## Here, parameters ending in `.1` correspond to the root partition
## and parameters ending in `.i` correspond to the clade subtended by
## node `i-1`.  We could constrain this with something like
## `constrain(lik, lambda.2 ~ lambda.1, lambda.3 ~ lambda.1)` but this
## quickly becomes tedious.

## Instead using the function `make.bisse.uneven` (new in diversitree
## 0.9-5) simplifies this by allowing varying sampling probabilities,
## but keeping one set of parameters:
lik.u <- make.bisse.uneven(phy2, states2,
                           sampling.f=sampling.f.p,
                           nodes=nodes)

## When all parameters are the same, this gives identical answers to
## repeating the parameters four times:
identical(lik.u(p),
          lik.p(rep(p, 4)))
## (this is because this is *actually what it is doing*!)

## This behaves just as a normal BiSSE likelihood function and can be
## optimised like so:
fit.u <- find.mle(lik.u, p)

## The coefficients here are remarkably close to the true values.
rbind(full=coef(fit), partial=coef(fit.p),
      sampled=coef(fit.s), uneven=coef(fit.u))[,c("lambda0", "lambda1")]

## This model has not been explored in detail statistically, so be
## careful!

## ## Document details

## Document compiled on `r as.character(Sys.time(), usetz=TRUE)`
## With R version `r R.version$string`
## and diversitree version `r packageVersion("diversitree")`.
## Original source: [uneven.R](uneven.R)

## [examples](..), [diversitree home](../..)

