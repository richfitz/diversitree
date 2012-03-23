## OU motion model from Luke's geiger package.
make.all.branches.ou.vcv <- function(cache, control) {
  n.tip <- cache$n.tip
  states <- cache$states
  states.sd <- cache$states.sd
  vcv <- cache$vcv
  one <- rep(1, n.tip)

  function(pars, intermediates) {
    s2 <- pars[1]
    alpha <- pars[2]
    vv <- s2 * ouMatrix(vcv, alpha)
    diag(vv) <- diag(vv) + states.sd^2
    mu <- phylogMean(vv, states)
    dmvnorm2(states, rep(mu, n.tip), vv, solve(vv), log=TRUE)
  }
}

## From geiger:
phylogMean<-function(vcv, data) {
  o <- rep(1, length(data))
  ci <- solve(vcv)
  ## By my calculations t(1) %*% VI %*% 1 = sum(VI)
  m1 <- solve(t(o) %*% ci %*% o) # == 1/sum(ci)
  ## t(one) %*% VI %*% states = sum(colSums(VI) * states)  
  m2 <- t(o) %*% ci %*% data
  ## sum(colSums(VI) * states) / sum(VI)
  m1 %*% m2
}

ouMatrix <- function(vcv, alpha) {
  ## follows Hansen 1997; does not assume ultrametricity (AH 12 dec
  ## 07) vectorized by LJH
  vcv.d <- diag(vcv)
  diagi <- matrix(vcv.d, length(vcv.d), length(vcv.d))
  diagj <- matrix(vcv.d, length(vcv.d), length(vcv.d), byrow=TRUE)
  Tij <- diagi + diagj - (2 * vcv)
  (1 / (2 * alpha)) * exp(-alpha * Tij) * (1 - exp(-2 * alpha * vcv))
}
