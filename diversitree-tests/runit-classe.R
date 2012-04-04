library(diversitree)
library(RUnit)

tol <- 2e-6

#--------------------------------------------------
# read in test data
#--------------------------------------------------

read.ttn <- function(treefile)
{
    treestr <- readLines(treefile, n=1)
    tree <- read.tree(text=treestr)
    ntips <- Ntip(tree)

    all.states <- read.table(treefile, skip=1)
    tip.states <- all.states[1:ntips, 2]
    names(tip.states) <- all.states[1:ntips, 1]

    return(list(tree = tree, states = tip.states))
}

# sim params = (1.4, 0.6, 0.2, 0.2, 0.5, 0), root state = 0
ttn1 <- read.ttn("bisse-tree.ttn")

# sim params = (1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), root state = 0
ttn2 <- read.ttn("geosse-tree.ttn")

# sim params (lam = 6.2, p = 0.2, q = 0; symmetric):
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#   4.96   1.24      0      0   1.24   4.96      3      0      0      0
ttn3 <- read.ttn("classe-tree.ttn")

#--------------------------------------------------
# prepare likelihood functions
#--------------------------------------------------

lnL1.bisse <- make.bisse(ttn1$tree, ttn1$states)

lnL1.classe <- make.classe(ttn1$tree, ttn1$states+1, 2)
lnL1.c.classe <- make.classe(ttn1$tree, ttn1$states+1, 2, control=list(backend="cvodes"))
lnL1.C.classe <- make.classe(ttn1$tree, ttn1$states+1, 2, control=list(backend="CVODES"))

lnL1.classe2 <- constrain(lnL1.classe, lambda112~0, lambda122~0, lambda211~0, lambda212~0)
lnL1.c.classe2 <- constrain(lnL1.c.classe, lambda112~0, lambda122~0, lambda211~0, lambda212~0)
lnL1.C.classe2 <- constrain(lnL1.C.classe, lambda112~0, lambda122~0, lambda211~0, lambda212~0)

lnL2.geosse <- make.geosse(ttn2$tree, ttn2$states)

lnL2.classe <- make.classe(ttn2$tree, ttn2$states+1, 3)
lnL2.c.classe <- make.classe(ttn2$tree, ttn2$states+1, 3, control=list(backend="cvodes"))
lnL2.C.classe <- make.classe(ttn2$tree, ttn2$states+1, 3, control=list(backend="CVODES"))

lnL2.classe2 <- constrain(lnL2.classe, lambda111~0, lambda122~0,
                            lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                            lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                            lambda313~0, lambda322~0, lambda323~0,
                            lambda112~lambda222, lambda113~lambda333, 
                            mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)
lnL2.c.classe2 <- constrain(lnL2.c.classe, lambda111~0, lambda122~0,
                            lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                            lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                            lambda313~0, lambda322~0, lambda323~0,
                            lambda112~lambda222, lambda113~lambda333, 
                            mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)
lnL2.C.classe2 <- constrain(lnL2.C.classe, lambda111~0, lambda122~0,
                            lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                            lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                            lambda313~0, lambda322~0, lambda323~0,
                            lambda112~lambda222, lambda113~lambda333, 
                            mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)

lnL3.classe <- make.classe(ttn3$tree, ttn3$states+1, 2)
lnL3.c.classe <- make.classe(ttn3$tree, ttn3$states+1, 2, control=list(backend="cvodes"))
lnL3.C.classe <- make.classe(ttn3$tree, ttn3$states+1, 2, control=list(backend="CVODES"))

lnL3.classe2 <- constrain(lnL3.classe, lambda122~0, lambda211~0)
lnL3.c.classe2 <- constrain(lnL3.c.classe, lambda122~0, lambda211~0)
lnL3.C.classe2 <- constrain(lnL3.C.classe, lambda122~0, lambda211~0)

# TODO: add tests against bisseness for ttn3

#--------------------------------------------------
# prepare parameter vectors
#--------------------------------------------------

pars1.bisse <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05)
pars1.classe <- c(rep(0, 6), pars1.bisse[-seq(2)])
names(pars1.classe) <- argnames(lnL1.classe)
pars1.classe[c('lambda111', 'lambda222')] <- pars1.bisse[1:2]

pars2.geosse <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars2.geosse) <- argnames(lnL2.geosse)
pars2.classe <- 0 * starting.point.classe(ttn2$tree, 3)
pars2.classe['lambda222'] <- pars2.classe['lambda112'] <- pars2.geosse['sA']
pars2.classe['lambda333'] <- pars2.classe['lambda113'] <- pars2.geosse['sB']
pars2.classe['lambda123'] <-  pars2.geosse['sAB']
pars2.classe['mu2'] <- pars2.classe['q13'] <- pars2.geosse['xA']
pars2.classe['mu3'] <- pars2.classe['q12'] <- pars2.geosse['xB']
pars2.classe['q21'] <- pars2.geosse['dA']
pars2.classe['q31'] <- pars2.geosse['dB']
pars2.geosse2 <- pars2.geosse[c(3,1,2,4:7)]

pars3.classe <- c(5, 1, 0.1, 0.2, 2, 4, 3, 2.5, 2.1, 2.2)
names(pars3.classe) <- argnames(lnL3.classe)
pars3.classe2 <- pars3.classe[-c(3,4)]

#--------------------------------------------------
# test functions
#--------------------------------------------------

# compare classe with bisse
test.lnL1 <- function()
{
    argvals <- list(condition.surv=T)
    ans <- -284.9406
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)

    argvals <- list(condition.surv=F)
    ans <- -284.9744
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)

    argvals <- list(condition.surv=T, root=ROOT.GIVEN, root.p=c(0.6, 0.4))
    ans <- -285.0867
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.c.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL1.C.classe2, c(list(pars1.bisse), argvals)), ans, tolerance=tol)

    argvals <- list(condition.surv=T, root=ROOT.EQUI)
    ans <- -285.2541
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans, tolerance=tol)
    ## Still has an error...
    checkEquals(do.call(lnL1.classe, c(list(pars1.classe), argvals)), ans, tolerance=tol)
}

# compare classe with geosse
test.lnL2 <- function()
{
    argvals <- list(condition.surv=T)
    ans <- -387.278
    checkEquals(do.call(lnL2.geosse, c(list(pars2.geosse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.c.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.C.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.c.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.C.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)

    argvals <- list(condition.surv=F, root=ROOT.GIVEN, root.p=c(0.5, 0.3, 0.2))
    ans <- -386.9637
    checkEquals(do.call(lnL2.geosse, c(list(pars2.geosse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.c.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.C.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.c.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.C.classe2, c(list(pars2.geosse2), argvals)), ans, tolerance=tol)

    argvals <- list(condition.surv=T, root=ROOT.EQUI)
    ans <- -387.1794
    checkEquals(do.call(lnL2.geosse, c(list(pars2.geosse), argvals)), ans, tolerance=tol)
    checkEquals(do.call(lnL2.classe, c(list(pars2.classe), argvals)), ans, tolerance=tol)
}

test.lnL3 <- function()
{
    ans <- 36.75782
    checkEquals(lnL3.classe(pars3.classe), ans, tolerance=tol)
    checkEquals(lnL3.c.classe(pars3.classe), ans, tolerance=tol)
    checkEquals(lnL3.C.classe(pars3.classe), ans, tolerance=tol)

    ans <- 35.53881
    checkEquals(lnL3.classe2(pars3.classe2), ans, tolerance=tol)
    checkEquals(lnL3.c.classe2(pars3.classe2), ans, tolerance=tol)
    checkEquals(lnL3.C.classe2(pars3.classe2), ans, tolerance=tol)

    checkEquals(lnL3.classe(pars3.classe, root=ROOT.EQUI), 36.75766, tolerance=tol)
    checkEquals(lnL3.classe2(pars3.classe2, root=ROOT.EQUI), 35.53928, tolerance=tol)
}

# no likelihoods, just parameters
test.params <- function()
{
    ans <- c( rep(3.058056, 6), rep(4.587084, 2), rep(4.587084, 2) )
    names(ans) <- argnames(lnL3.classe)
    checkEquals(starting.point.classe(ttn3$tree, 2), ans, tolerance=tol)

    checkEquals(diversitree:::stationary.freq.classe(pars1.classe, 2), 
                diversitree:::stationary.freq.bisse(pars1.bisse), tolerance=tol)

    checkEquals(diversitree:::stationary.freq.classe(pars2.classe, 3), 
                diversitree:::stationary.freq.geosse(pars2.geosse), tolerance=tol)

    pars <- seq_len(27)  # for k = 3 states
    names(pars) <- diversitree:::default.argnames.classe(3)
    parlist <- diversitree:::inflate.pars.classe(pars, 3)
    checkEquals(parlist$lambda[3,1,2], 14, tolerance=tol)
    checkEquals(parlist$mu[2], 20, tolerance=tol)
    checkEquals(parlist$q[1,3], 23, tolerance=tol)
    checkIdentical(pars, diversitree:::flatten.pars.classe(parlist))

    # (not really parameters, but related to stationary.freq.classe)
    ans <- matrix(c(-78, 37, 43, 55, -87, 69, 81, 90, -99), nrow=3)
    checkEquals(diversitree:::projection.matrix.classe(pars, 3), ans)
}
