## It might be the case that auxilliary variables need computing
## at splits.  For example, in BiSSE (and friends), when switching
## parameter space, rather than take the E values from the
## daughter lineage, I need to recompute them to use the new
## parameters at this point (a lineage arising just below the
## split would have a totally different probability of extinction
## than one after the split).  In this case I have to run a new
## branch down with the parent set and parent sampling.f.  The 'D'
## values from doing this would be ignored.  For mk2 (and friends)
## this is not an issue.

## If parameters are shared for some partitions across calls, there is
## no need to recompute things.  However, going beyond the simple
## implementation below is a lot of book-keeping.  The reason for
## this, is that if we want to remember multiple points previously
## (which would be ideal given the way that MCMC and MLE searches
## work, generally modifying a single partition at once) then we have
## to remember parent/offspring parameter relationships.
## Suppose a tree has nested sements A-B-C (A tipward, B rootward).
## and suppose that the following parameter sets have been run:
##    A1 B1 C1
##    A2 B2 C2
## and now we want to run
##    A2 B2 C1
## we would mark the A and B partitions as case=0, and mark the C
## partition as case=2, starting from the C2 initial conditions.
## If we then evaluate
##    A1 B2 C2
## we would mark things 0, 2, 2 (the C=2 would happen automatically)
## and then
##    A1 B2 C1 -> 0, 0, 2
## However, this would require that we are storing a large number of
## intermediates, which might not be feasble anyway for QuaSSE.

## This is surprisingly slow for QuaSSE (can take about 0.1s).
## recycle.get <- function(cache, pars) {
##   e <- cache$prev
##   desc <- cache$desc.parts

##   n <- length(pars)
##   run <- rep(1L, n)
##   vals <- vector("list", n)
##   take <- integer(n)

##   for ( i in seq_len(n) ) {
##     tmp <- e$res[[i]]
##     for ( j in seq_along(tmp) ) {
##       x <- tmp[[j]]
##       if ( identical(x$pars, pars[[i]]) ) {
##         run[i] <- if ( identical(x$extra, pars[desc[[i]]]) )
##           0L else 2L
##         vals[[i]] <- x
##         take[i] <- j
##         break
##       }
##     }
##   }

##   list(run=run, vals=vals, take=take)
## }

## recycle.set <- function(cache, pars, vals, run, desc) {
##   e <- cache$prev
##   desc <- cache$desc.parts

##   n <- length(pars)
##   for ( i in seq_len(n) ) {
##     if ( run[i] == 1 ) {
##       if ( !is.null(keep <- cache$cache[[i]]$recycle.keep) ) {
##         tmp <- vals[[i]]$intermediates[c("base", "lq")]
##         tmp$base[-keep] <- list(NULL)
##         vals[[i]]$intermediates <- tmp
##       } else {
##         ## This was <- NULL to blank the intermediates.
##         vals[[i]]$intermediates <- vals[[i]]$intermediates["lq"]
##       }

##       x <- list(pars=pars[[i]], vals=vals[[i]],
##                 extra=pars[desc[[i]]])
      
##       e$res[[i]] <- push.stack(e$res[[i]], x)
##     }
##   }
## }

## ## This is not really a stack, but it will do for now.
## make.stack <- function(n) {
##   structure(list(),
##             length=n,
##             class="stack")
## }

## push.stack <- function(x, val) {
##   n <- attr(x, "length")
##   ret <- c(list(val), if (length(x) >= n) x[seq_len(n-1)] else x)
##   attr(ret, "length") <- n
##   class(ret) <- "stack"
##   ret
## }

