## Time chunks.

## The make.branches.td function converts a normal "branches" function
## into a time-dependent function.

## The returned function takes an argument list of the format:

## This makes a time-dependent *anything* by wrapping up the branches
## function.

## Much of the difficulty here is in trying to This is somewhat complicated by trying to do the tips with multiple
## different times at once.

## This function may change in the future to accomodate QuaSSE (QuaSSE
## is difficult because the parameters are a list rather than a
## vector, and the returned variables can be of different lengths
## depending on the time).
make.branches.td <- function(branches) {
  function(y, len, pars, t0) {
    t <- pars[,1]
    bpars <- pars[,-1]

    t1 <- t0 + len
    lq <- 0

    i <- which(t > t0)[1]
    j <- which(t > t1[length(t1)])[1]

    done <- rep(FALSE, length(len))
    out <- vector("list", length(t))

    for ( epoch in i:j ) {
      not.last <- epoch < j
      k <- t1 < t[epoch] & !done
      times <- c(t1[k], if (not.last) t[epoch])
      n <- length(times)

      z <- branches(y, times-t0, bpars[epoch,], t0)[,,drop=FALSE]
      z[,1] <- z[,1] + lq

      if ( any(k) ) {
        out[[epoch]] <- if ( not.last ) z[-n,] else z
        done[k] <- TRUE
      }
      
      if ( not.last ) {
        y <- z[n,-1]
        lq <- z[n,1]
        t0 <- t[epoch]
      }
    }

    do.call(rbind, out)
  }
}

make.initial.conditions.td <- function(initial.conditions)
  function(init, pars, t, is.root=FALSE)
  initial.conditions(init, get.par.td(pars, t), t, is.root)

get.par.td <- function(pars, t)
  pars[which(pars[,1] > t)[1],-1]

## This is identical to the version in diversitree-branches.R, except
## for the root parameter treatment.
ll.xxsse.td <- function(pars, cache, initial.conditions,
                        branches, condition.surv, root, root.p,
                        intermediates) {
  pars.root <- get.par.td(pars, cache$depth[cache$root])  
  ans <- all.branches(pars, cache, initial.conditions, branches)
  vals <- ans$init[[cache$root]]
  root.p <- root.p.xxsse(vals, pars.root, root, root.p)
  loglik <- root.xxsse(vals, pars.root, ans$lq, condition.surv,
                       root.p)
  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}
