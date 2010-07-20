## This is surprisingly slow for QuaSSE (can take about 0.1s).
recycle.get <- function(cache, pars) {
  e <- cache$prev
  desc <- cache$desc.parts

  n <- length(pars)
  run <- rep(1L, n)
  vals <- vector("list", n)
  take <- integer(n)

  for ( i in seq_len(n) ) {
    tmp <- e$res[[i]]
    for ( j in seq_along(tmp) ) {
      x <- tmp[[j]]
      if ( identical(x$pars, pars[[i]]) ) {
        run[i] <- if ( identical(x$extra, pars[desc[[i]]]) )
          0L else 2L
        vals[[i]] <- x
        take[i] <- j
        break
      }
    }
  }

  list(run=run, vals=vals, take=take)
}

recycle.set <- function(cache, pars, vals, run, desc) {
  e <- cache$prev
  desc <- cache$desc.parts

  n <- length(pars)
  for ( i in seq_len(n) ) {
    if ( run[i] == 1 ) {
      if ( !is.null(keep <- cache$cache[[i]]$recycle.keep) ) {
        tmp <- vals[[i]]$intermediates[c("base", "lq")]
        tmp$base[-keep] <- list(NULL)
        vals[[i]]$intermediates <- tmp
      } else {
        ## This was <- NULL to blank the intermediates.
        vals[[i]]$intermediates <- vals[[i]]$intermediates["lq"]
      }

      x <- list(pars=pars[[i]], vals=vals[[i]],
                extra=pars[desc[[i]]])
      
      e$res[[i]] <- push.stack(e$res[[i]], x)
    }
  }
}

## This is not really a stack, but it will do for now.
make.stack <- function(n) {
  structure(list(),
            length=n,
            class="stack")
}

push.stack <- function(x, val) {
  n <- attr(x, "length")
  ret <- c(list(val), if (length(x) >= n) x[seq_len(n-1)] else x)
  attr(ret, "length") <- n
  class(ret) <- "stack"
  ret
}

