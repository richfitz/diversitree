library(diversitree)
library(RUnit)

tol <- 1e-6

#--------------------------------------------------
# test data for everything below
#--------------------------------------------------

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")
states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- as.character(seq(Ntip(tree))-1)

lnL.full <- make.geosse(tree, states)
lnL.full.c <- make.geosse(tree, states, control=list(backend="cvodes"))
lnL.full.C <- make.geosse(tree, states, control=list(backend="CVODES"))

rate.names <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
pars0 <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)

#--------------------------------------------------
# test functions
#--------------------------------------------------

test.starting <- function()
{
    ans <- starting.point.geosse(tree)
    ans0 <- c( rep(3.339799, 3), rep(1.669899, 4) )
    names(ans0) <- rate.names
    checkEquals(ans, ans0, tolerance=tol)

    ans <- starting.point.geosse(tree, eps=0)
    ans0 <- c( rep(1.8861312, 3), rep(0, 2), rep(0.1886131, 2) )
    names(ans0) <- rate.names
    checkEquals(ans, ans0, tolerance=tol)
}

test.lnL <- function()
{
    pars <- pars0
    checkEquals(lnL.full(pars), -24.07128, tolerance=tol)
    checkEquals(lnL.full(pars, condition.surv=F), -23.71574, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.EQUI), -24.25509, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.GIVEN, root.p=c(0.6,0.4,0.2)), 
                -24.1432, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.GIVEN, root.p=rep(1,3)/3), 
                lnL.full(pars, root=ROOT.FLAT), tolerance=tol)

    lnL <- constrain(lnL.full, sAB ~ 0)
    checkEquals(lnL(pars[-3]), -23.86826, tolerance=tol)

    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    checkEquals(lnL(pars[-c(3,7)]), -24.05792, tolerance=tol)
}

test.cvodes <- function()
{
    pars <- pars0
    checkEquals(lnL.full(pars), lnL.full.c(pars), tolerance=tol)
    checkEquals(lnL.full(pars), lnL.full.C(pars), tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.EQUI), lnL.full.c(pars, root=ROOT.EQUI), tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.EQUI), lnL.full.C(pars, root=ROOT.EQUI), tolerance=tol)

    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    lnL.c <- constrain(lnL.full.c, dB ~ dA, sAB ~ 0)
    lnL.C <- constrain(lnL.full.C, dB ~ dA, sAB ~ 0)
    checkEquals(lnL(pars[-c(3,7)]), lnL.c(pars[-c(3,7)]), tolerance=tol)
    checkEquals(lnL(pars[-c(3,7)]), lnL.C(pars[-c(3,7)]), tolerance=tol)
}

test.mle <- function()
{
    pars <- pars0
    names(pars) <- argnames(lnL.full)
    ans <- find.mle(lnL.full, pars)
    ans0.par <- c(1.483118e+00, 3.653205e-01, 2.478823e-07, 2.010653e-05, 
                  3.944644e-06, 1.275225e+00, 1.206178e+00)
    names(ans0.par) <- rate.names
    checkEquals(coef(ans), ans0.par, tolerance=tol)
    checkEquals(logLik(ans)[1], -19.32895, tolerance=tol)

    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    pars <- pars0[-c(3,7)]
    names(pars) <- argnames(lnL)
    ans <- find.mle(lnL, pars)
    ans0.par <- c(1.481470e+00, 3.698169e-01, 9.211896e-09, 1.526980e-08, 1.270324e+00)
    names(ans0.par) <- rate.names[-c(3,7)]
    checkEquals(coef(ans), ans0.par, tolerance=tol)
    checkEquals(logLik(ans)[1], -19.32908, tolerance=tol)
}

test.mcmc <- function()
{
    set.seed(1)
    pars <- pars0
    ans <- mcmc(lnL.full, pars, nsteps=3, lower=0, upper=3, w=3/5, 
                prior=make.prior.exponential(seq(0.8, by=0.1, length.out=7)),
                print.every=0) 
    checkEquals(ans$p[3], -30.56511, tolerance=tol)

    set.seed(1)
    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    pars <- pars0[-c(3,7)]
    ans <- mcmc(lnL, pars, nsteps=3, lower=0, upper=3, w=seq(0.8, by=0.1,
                length.out=5), prior=make.prior.exponential(1), print.every=0)
    checkEquals(ans$p[3], -28.18477, tolerance=tol)
}

test.split <- function()
{
    # with 0.8-4, only split.t = Inf works (for any model)
    lnL.split <- make.geosse.split(tree, states, c(27, 29), split.t = Inf)
    checkEquals(lnL.full(pars0), lnL.split(rep(pars0, 3)), tolerance=tol)
    checkEquals(lnL.split(c(pars0, pars0*1.5, pars0*0.5)), -23.82277, tolerance=tol)
}

# simulator
test.sim <- function()
{
    pars.g <- pars0
    names(pars.g) <- diversitree:::default.argnames.geosse()
    pars.c <- diversitree:::pars.ge.to.cl(pars.g)

    set.seed(1)
    phy <- trees(pars.g, type="geosse", n=2, max.t=4, x0=0)
    checkEquals(lapply(phy, Ntip), list(20, 132))
    lnL <- make.geosse(phy[[2]], phy[[2]]$tip.state)
    checkEquals(lnL(pars.g), -252.7173, tolerance=tol)

    set.seed(3)
    phy <- tree.geosse(pars.g, max.t=5)
    checkEquals(as.numeric(table(phy$tip.state)), c(31, 54, 33))
    lnL.g <- make.geosse(phy, phy$tip.state)

    set.seed(3)
    phy2 <- tree.classe(pars.c, max.t=5)
    lnL.c <- make.classe(phy2, phy2$tip.state, 3)
    checkEquals(lnL.c(pars.c), lnL.g(pars.g), tolerance=tol)
}
