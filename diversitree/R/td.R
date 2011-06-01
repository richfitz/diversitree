## Time chunks.

## The make.branches.td function converts a normal "branches" function
## into a time-dependent function.
make.branches.td <- function(branches) {
  function(y, len, pars, t0) {
    t <- pars[,1]
    bpars <- pars[,-1]

    t1 <- t0 + len
    lq.carry <- 0

    i <- which(t > t0)[1]
    j <- which(t > t1[length(t1)])[1]

    done <- rep(FALSE, length(len))
    out <- vector("list", length(t))

    base <- matrix(NA, length(y), length(len))
    lq <- numeric(length(len))

    for ( epoch in i:j ) {
      not.last <- epoch < j
      ## Identify branches that terminate in this epoch:
      k <- t1 < t[epoch] & !done
      
      times <- c(t1[k], if (not.last) t[epoch])
      n <- length(times)

      res <- branches(y, times-t0, bpars[epoch,], t0)
      res[[1]] <- res[[1]] + lq.carry

      if ( any(k) ) {
        base[,k] <- if ( not.last ) res[[2]][,-n] else res[[2]]
        lq[k]    <- if ( not.last ) res[[1]][-n]  else res[[1]]
        done[k] <- TRUE
      }
      
      if ( not.last ) {
        y <- res[[2]][,n]
        lq.carry <- res[[1]][n]
        t0 <- t[epoch]
      }
    }

    list(lq, base)
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
  
  ans <- all.branches.matrix(pars, cache, initial.conditions, branches)
  vals <- ans$init[,cache$root]
  root.p <- root.p.xxsse(vals, pars.root, root, root.p)
  loglik <- root.xxsse(vals, pars.root, ans$lq, condition.surv,
                       root.p)

  if ( intermediates ) {
    ans$root.p <- root.p
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- vals
  }

  loglik  
}
