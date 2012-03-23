make.all.branches.ou.direct <- function(cache, control) {
  function(pars, intermediates, preset=NULL)
    all.branches.matrix(pars, cache,
                        initial.conditions.bm.direct,
                        branches.ou, preset)
}

## In ou, for a branch that ends with mean m (y[1]), the mean moves to
##   exp(len * alpha) * (m - theta) + theta
## And the variance becomes
##   (exp(2 * len * alpha) - 1) / (2 * alpha) + exp(2 * len * alpha) * v
## With normalising constant
##   exp(t * alpha)

## The last line of the second option comes from the limit of
##   (exp(2*len*alpha) - 1) * sigma2 / (2*alpha)
## as alpha -> 0 being len * sigma2
##   if ( alpha > 0 )
##     c(len * alpha + z,
##       exp(len * alpha) * (m - theta) + theta,
##       (exp(2*len*alpha) - 1) * sigma2 / (2*alpha) + exp(2*len*alpha) * v)
##   else
##     c(len * alpha + z,
##       exp(len * alpha) * (m - theta) + theta,
##       len * sigma2 + exp(2*len*alpha) * v)
branches.ou <- function(y, len, pars, t0, idx) {
  m <- y[1]
  v <- y[2]
  z <- y[3]

  sigma2 <- pars[1]
  alpha  <- pars[2]
  theta  <- pars[3]

  if ( alpha > 0 )
    list(len * alpha + z,
         c(exp(len * alpha) * (m - theta) + theta,
           (exp(2*len*alpha) - 1) * sigma2 / (2*alpha) +
           exp(2*len*alpha) * v,
           0))
  else
    list(len * alpha + z,
         c(exp(len * alpha) * (m - theta) + theta,
           len * sigma2 + v,
           0))
}
