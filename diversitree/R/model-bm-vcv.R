## Brownian motion model from Luke's geiger package, so that I can use
## this without loading the entire package (and also have some more
## fun with different models in the future).

## This is actually slightly slower than Luke's version at the
## moment, which might be due to the shift from a log basis to a
## linear basis (the log basis is probably "easier", linearising the
## effects and helping optimisation).
make.all.branches.bm.vcv <- function(cache, control) {
  n.tip <- cache$n.tip
  states <- cache$states
  states.sd <- cache$states.sd
  vcv <- cache$vcv
  one <- rep(1, n.tip)

  if ( all(states.sd == 0) ) {
    VI.tmp <- solve(vcv)
    function(x, intermediates) {
      VI <- VI.tmp / x
      ## By my calculations t(1) %*% VI %*% 1 = sum(VI)
      ## t(one) %*% VI %*% states = sum(colSums(VI) * states)
      ## mu <- solve(t(one) %*% VI %*% one) %*% (t(one) %*% VI %*% states)
      mu <- sum(colSums(VI) * states) / sum(VI)
      dmvnorm2(states, rep(mu, n.tip), vcv * x, VI, log=TRUE)
    }
  } else {
    function(x, intermediates) {
      vv <- x*vcv
      diag(vv) <- diag(vv) + states.sd^2
      VI <- solve(vv)
      mu <- sum(colSums(VI) * states) / sum(VI)
      dmvnorm2(states, rep(mu, n.tip), vv, VI, log=TRUE)
    }
  }
}

rootfunc.bm.vcv <- function(res, pars, root, root.x, intermediates) {
  if ( root != ROOT.MAX )
    stop('root cannot be modified -- use method="direct"')
  if ( intermediates )
    stop('intermediates cannot be produced -- use method="direct"')
  res
}
