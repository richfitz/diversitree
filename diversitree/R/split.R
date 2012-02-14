## Check that the splits are OK.
check.split <- function(phy, nodes, split.t) {
  if ( length(split.t) == 1 && length(nodes) > 1 ) {
    if ( split.t == 0 || split.t == Inf )
      split.t <- rep(split.t, length(nodes))
    else
      stop("If split.t is length 1, it must be '0' or 'Inf'")
  } else if ( length(split.t) != length(nodes) ) {
    stop("'nodes' and 'split.t' must be the same length")
  }
    
  ## Check that all nodes are ok
  n.tip <- length(phy$tip.label)
  if ( is.character(nodes) )
    nodes <- match(nodes, phy$node.label) + n.tip
  if ( any(is.na(nodes) | nodes <= n.tip+1 | nodes > n.tip + phy$Nnode) )
    stop("Invalid node specification")
  
  ## Check that split times are OK
  partial <- !(split.t == Inf | split.t == 0)
  if ( any(partial) )
    stop("Partial split times not yet allowed")

  list(nodes=c(NA, nodes),
       split.with.edge=c(NA, split.t == Inf))
}

## Determine the group membership of every edge.
## Augment the usual cache vector with some extra information:
make.cache.split <- function(phy, cache, nodes, split.t) {
  tmp <- split.group(phy, nodes, split.t)
  nodes <- tmp$node

  cache$nodes <- tmp$nodes
  cache$group.nodes <- tmp$group.nodes
  cache$group.branches <- tmp$group.branches
  cache$split.with.edge <- tmp$split.with.edge

  cache$n.part <- length(nodes)
  cache$aux.use <- rep(FALSE, length(cache$group.branches))
  cache$aux.use[nodes[-1]] <- TRUE

  cache
}

split.group <- function(phy, nodes, split.t) {
  tmp <- check.split(phy, nodes, split.t)
  nodes <- tmp$nodes
  split.with.edge <- tmp$split.with.edge
  n.tip <- length(phy$tip.label)
  edge <- phy$edge

  descendants.idx <- function(node)
    which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
  
  f <- function(node, phy, group) {
    if ( is.na(node) )
      return(list(group=group, base=NA))
    i <- which(edge[,1] == node | edge[,2] %in% descendants(node, edge))
    base <- group[edge[,2] == node]
    group[i[group[i] == base]] <- max(group) + 1
    list(group=group, base=base)
  }

  group <- rep(1, nrow(edge))

  base <- integer(length(nodes))
  for ( i in seq_along(nodes) ) {
    tmp <- f(nodes[i], phy, group)
    group <- tmp$group
    base[i] <- tmp$base
  }

  group <- group[match(seq_len(max(edge)), edge[,2])]
  group[n.tip + 1] <- 1

  group.branches <- group
  i <- which(!split.with.edge)
  if ( length(i) > 0 )
    group.branches[which(i)] <- base[which(i)]

  list(nodes=c(n.tip+1, nodes[-1]),
       group.nodes=group, group.branches=group.branches,
       split.with.edge=split.with.edge)
}

make.branches.split <- function(cache, branches, branches.aux) {
  group <- cache$group.branches
  split.with.edge <- cache$split.with.edge
  aux.use <- cache$aux.use
  aux.i <- cache$aux.i

  branches.c     <- make.caching.branches(cache, branches)
  branches.aux.c <- make.caching.branches.aux(cache, branches.aux)

  function(y, len, pars, t0, idx) {
    g <- group[idx]

    if ( length(len) == 1 ) {
      ## cat(sprintf("%d: %d, len=%2.5f, t0=%2.5f\n", idx, g, len, t0))
      ## 1: This might be an internal edge, so it is possible that
      ## there are auxilliary variables that need computing.
      ## (this could also be a terminal edge where there is only a
      ## single group)
      p <- pars[[g]]

      if ( !aux.use[idx] ) {
        ## Normal case: just run the branch
        branches.c(y, len, p, t0, idx)
      } else if ( split.with.edge[g] ) {
        ## This is a join.  In this one here, the edge also has the
        ## derived character state:
        ans <- branches.c(y, len, p, t0, idx)

        g.parent <- group[cache$parent[idx]]
        aux <- branches.aux.c(g.parent, t0+len, pars[[g.parent]], idx)
        if ( is.list(ans) )
          ans[[2]][aux.i] <- aux
        else
          ans[-1][aux.i] <- aux

        ans # return
      } else {
        ## Also a join, but now the edge is in the ancestral character
        ## state:
        y[aux.i] <- branches.aux.c(g, t0, p, idx)
        branches.c(y, len, p, t0, idx)
      }
    } else if ( length(unique(g)) == 1 ) {
      ## 2: This must be a set of tips, but it is only of one type, so
      ## let's do it simply.
      branches.c(y, len, pars[[g[1]]], t0, idx)
    } else {
      ## 3: This is a set of tips that have more than one group type
      ## in it, so loop over the different group types.
      grp <- sort.default(unique.default(group[idx]))
      i <- split.default(seq_along(idx), group[idx])

      if ( is.null(cache$vars.in.list) ) {
        lq <- numeric(length(len))
        ans <- matrix(NA, cache$ny, length(len))

        for ( g.i in seq_along(grp) ) {
          j <- i[[g.i]]
          tmp <- branches.c(y, len[i[[g.i]]], pars[[grp[g.i]]], t0, idx[j])
          lq[j] <- tmp[[1]]
          ans[,j] <- tmp[[2]]
        }
        list(lq, ans)
      } else {
        ans <- vector("list", length(len))
        for ( g.i in seq_along(grp) ) {
          j <- i[[g]]
          ans[j] <- branches.c(y, len[i[[g]]], pars[[grp[g]]], t0, idx[j])
        }
        ans
      }
    }
  }
}

make.initial.conditions.split <- function(cache, initial.conditions) {
  group <- cache$group.nodes
  function(init, pars, t, idx)
    initial.conditions(init, pars[[group[idx]]], t, idx)
}

ll.xxsse.split <- function(pars, cache, initial.conditions,
                           branches, condition.surv, root, root.p,
                           intermediates) {
  ## always pars[[1]], but being safe...
  pars.root <- pars[[cache$group.nodes[cache$root]]]

  if ( cache$control$caching.branches )
    caching.branches.set.pars(pars, branches)

  ans <- all.branches.matrix(pars, cache, initial.conditions, branches)
  vals <- ans$init[,cache$root]
  root.p <- root.p.xxsse(vals, pars.root, root, root.p)
  loglik <- root.xxsse(vals, pars.root, ans$lq, condition.surv, root.p)

  if ( intermediates ) {
    ans$root.p <- root.p
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- vals
  }

  loglik  
}

make.initial.tip.xxsse.split <- function(cache) {
  k <- cache$k
  n.part <- cache$n.part
  grp <- cache$group.branches  
  idx.D <- (k+1):(2*k)

  out <- vector("list", length(cache$y))
  for ( i in seq_along(cache$y)) {
    y.i <- cache$y[[i]]
    out[[i]] <- vector("list", n.part)
    
    for ( group in seq_len(n.part) ) {
      keep <- grp[y.i$target] == group
      if ( any(keep) ) {
        s.f <- cache$sampling.f[[group]]
        out[[i]][[group]] <- 
          list(y     = c(1-s.f, s.f * y.i$y[idx.D]),
               y.i   = y.i$y.i,
               target= y.i$target[keep],
               t     = y.i$t[keep])
      }
    }
  }
  unlist(out, FALSE)
}

## The first of these that I'll do is fully passive.  It checks to
## see if the parameters and input are identical.  This should be
## the most basic, and most foolproof.  It may not be the fastest
## though, as some of these checks add overhead.  The use of
## identical should keep things fairly fast though.

## TODO: We only really need to bother checking the variables for
## cases where the group has a descendent.

## A caching all.branches:
all.branches.caching <- function(pars, cache, initial.conditions,
                                 branches, branches.aux,
                                 type="matrix") {
  if ( type == "matrix" )
    res <- all.branches.matrix(pars, cache, initial.conditions, branches)
  else
    res <- all.branches.list(pars, cache, initial.conditions, branches)

  ## Now, update the environment of the branches and branches.aux
  ## functions with the parameters used and the results:
  e <- 
  e.aux <- environment(branches.aux)

  environment(branches)$prev.pars <-
    environment(branches.aux)$prev.pars <- pars

  res
}


## I want to make caching versions of the branches function.
caching.branches.set.pars <- function(p, branches, branches.aux) {
  e <- environment(branches)
  e.b <- environment(e$branches.c)
  e.a <- environment(e$branches.aux.c)

  e.a$pars.same <- e.b$pars.same <- i <-
    mapply(identical, p, e.b$prev.pars)
  e.a$prev.pars <- e.b$prev.pars <- p

  i
}

make.caching.branches <- function(cache, branches) {
  if ( !cache$control$caching.branches )
    return(branches)

  cat("Using experimental caching branches!\n")
  
  n <- length(cache$len)
  n.part <- cache$n.part
  group <- cache$group.branches

  pars.same <- logical(n.part)
  prev.pars <- vector("list", n.part)
  prev.lq <- numeric(n)
  
  if ( is.null(cache$vars.in.list) ) {
    prev.vars <- prev.base <- matrix(NA, cache$ny, n)
    branches.caching <- function(y, len, p, t0, idx) {
      g <- group[idx[1]]
      if ( pars.same[g] && identical(prev.vars[,idx[1]], y) ) {
        list(prev.lq[idx],
             prev.base[,idx,drop=FALSE])
      } else {
        prev.vars[,idx] <<- y
        ret <- branches(y, len, p, t0, idx)
        prev.lq[idx]    <<- ret[[1]]
        prev.base[,idx] <<- ret[[2]]
        ret
      }
    }
  } else { # Just QuaSSE, really.
    prev.vars <- prev.base <- vector("list", n)
    branches.caching <- function(y, len, p, t0, idx) {
      g <- group[idx[1]]
      if ( length(idx) > 1 ) # Will take some work, but not used atm.
        stop("vector idx not yet allowed")
      if ( pars.same[g] && identical(prev.vars[[idx]], y) ) {
        c(prev.lq[idx],
          prev.base[[idx]])
      } else {
        prev.vars[[idx]] <<- y
        ret <- branches(y, len, p, t0, idx)
        prev.lq[idx]   <<- ret[1]
        prev.base[[idx]] <<- ret[-1]
        ret
      }
    }
  }

  branches.caching
}

## Note that for now, the returned function here includes an
## additional argument to the input function:
##   normal aux: (g, t, p)
##   returned:   (g, t, p, idx)
make.caching.branches.aux <- function(cache, branches.aux) {
  if ( !cache$control$caching.branches )
    return(function(g, t, p, idx) branches.aux(g, t, p))
  
  n <- length(cache$len)
  prev.aux <- vector("list", n)
  pars.same <- logical(cache$n.part)

  branches.aux.caching <- function(g, t, p, idx) {
    if ( pars.same[g] ) {
      prev.aux[[idx]]
    } else {
      prev.aux[[idx]] <<- ret <- branches.aux(g, t, p)
      ret
    }
  }

  branches.aux.caching
}
