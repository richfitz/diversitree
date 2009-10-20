## I need to tidy this up a little bit to allow for different tree
## types.  I also cannot use functions beginning with simulate(), as
## this is a standard R generic function.
##
## It might be useful to have the simulations use somthing like the
## equilibrium distribution for characters at the bottom of the tree.
##
## It is also worth noting that Luke Harmon has a birthdeath.tree
## function that simulates a tree under a birth death process in
## geiger.

## There is also some missing logic with single branch trees, that I
## need to work through.

## Main interface.  In the hope that I will make this generic over a
## 'model' object, I will design the calling structure in a way that
## is similar to S3 generics/methods.
trees <- function(pars, type=c("bisse", "birthdeath"), n=1,
                  max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
                  ...) {
  if ( is.infinite(max.taxa) && is.infinite(max.t) )
    stop("At least one of max.taxa and max.t must be finite")
  type <- match.arg(type)
  f <- switch(type,
              bisse=tree.bisse,
              birthdeath=tree.birthdeath)
  trees <- vector("list", n)  
  i <- 1

  while ( i <= n ) {
    trees[[i]] <- phy <-
      f(pars, max.taxa, max.t, include.extinct, ...)
    if ( include.extinct || !is.null(phy) )
      i <- i + 1
  }

  trees
}

tree.bisse <- function(pars, max.taxa=Inf, max.t=Inf,
                       include.extinct=FALSE, x0=NA, ...) {
  if ( is.na(x0) )
    x0 <- as.integer(runif(1) > stationary.freq(pars))
  else if ( length(x0) != 1 || !(x0 == 0 || x0 == 1) )
    stop("Invalid root state")
  stopifnot(length(pars) == 6 && all(pars >= 0))
  
  info <- make.tree.bisse(pars, max.taxa, max.t, x0)
  phy <- me.to.ape.bisse(info[-1,])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

tree.birthdeath <- function(pars, max.taxa=Inf, max.t=Inf,
                            include.extinct=FALSE, ...) {
  info <- make.tree.birthdeath(pars, max.taxa, max.t)
  phy <- me.to.ape.birthdeath(info[-1,])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

## This is written in an odd style because I want to shift this into C
## soon after getting it working.
## TODO: I think the single.lineage bit is broken, so I am not
## exposing it in the next function up (the two-lineage case is
## broken).  This will interact particularly with the
##   me.to.ape.bisse(info[-1])
## line in tree.bisse() above.
make.tree.bisse <- function(pars, max.taxa=Inf, max.t=Inf, x0,
                            single.lineage=TRUE) {
  pars <- matrix(pars, 2, 3)

  extinct <- FALSE
  split   <- FALSE
  parent <- 0

  n.i <- c(0, 0)
  r.i <- rowSums(pars)
  len <- 0
  t <- 0
  
  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0+1] <- 1
  } else {
    states <- rep(x0, 2)
    n.taxa <- lineages <- n.i[x0+1] <- 2
  }

  while ( n.taxa < max.taxa && n.taxa > 0) {
    ## When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Proceed.  What state does an event happen to?
    state <- as.integer(runif(1) > r.n[1]/r.tot)
    state.i <- state + 1

    ## Pick a lineage for that state:
    j <- sample(n.i[state.i], 1)
    lineage <- lineages[states[lineages] == state][j]

    type <- sample(3, 1, FALSE, pars[state.i,])

    if ( type == 1 ) {
      ## Speciating:
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      states[new.i] <- state
      parent[new.i] <- lineage
      len[new.i] <- 0

      n.i[state.i] <- n.i[state.i] + 1
      n.taxa <- n.taxa + 1

      lineages <- which(!split & !extinct)

    } else if ( type == 2 ) {
      ## Extinct
      extinct[lineage] <- TRUE

      lineages <- which(!split & !extinct)

      n.i[state.i] <- n.i[state.i] - 1
      n.taxa <- n.taxa - 1
    } else {
      ## Character switch:
      n.i <- n.i + if ( state == 0 ) c(-1,1) else c(1,-1)
      states[lineage] <- 1 - state
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     state=states, extinct=extinct, split=split)
  attr(info, "t") <- t
  info
}

me.to.ape.bisse <- function(x) {
  if ( nrow(x) == 0 )
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)

  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[ x$split] <- order(x$idx[x$split]) + n.tips + 1

  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  x$name <- NA
  x$name[!x$split] <- tip.label

  tip.state <- x$state[match(1:n.tips, x$idx2)]
  names(tip.state) <- tip.label
  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                               Nnode=Nnode,
                               tip.label=tip.label,
                               tip.state=tip.state,
                               node.label=node.label,
                               edge.length=x$len,
                               orig=x),
                          class="phylo"))

  phy$edge.state <- x$state[match(phy$edge[,2], x$idx2)]
  phy
}

make.tree.birthdeath <- function(pars, max.taxa=Inf, max.t=Inf) {
  extinct <- FALSE
  split   <- FALSE
  parent <- 0

  lambda <- pars[1]
  mu <- pars[2]
  r <- lambda + mu
  len <- 0
  t <- 0
  n.taxa <- 1
  lineages <- 1

  pr.speciation <- lambda/(lambda + mu)

  while ( n.taxa < max.taxa && n.taxa > 0) {
    ## When does an event happen?
    r.n <- r * n.taxa
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
      break
    }

    len[lineages] <- len[lineages] + dt

    ## Pick a lineage for the event to happen to:
    lineage.i <- sample(n.taxa, 1)
    lineage <- lineages[lineage.i]

    if ( runif(1) < pr.speciation ) {
      ## Speciating:
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE
      parent[new.i] <- lineage
      len[new.i] <- 0

      n.taxa <- n.taxa + 1

      ## lineages <- which(!split & !extinct)
      lineages <- c(lineages[-lineage.i], new.i)
    } else {
      ## Extinct
      extinct[lineage] <- TRUE
      ## lineages <- which(!split & !extinct)
      lineages <- lineages[-lineage.i]
      n.taxa <- n.taxa - 1
    }

  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     extinct=extinct, split=split)
  attr(info, "t") <- t
  info
}

me.to.ape.birthdeath <- function(x) {
  if ( nrow(x) == 0 )
    return(NULL)
  Nnode <- sum(!x$split) - 1
  n.tips <- sum(!x$split)

  x$idx2 <- NA
  x$idx2[!x$split] <- 1:n.tips
  x$idx2[ x$split] <- order(x$idx[x$split]) + n.tips + 1

  i <- match(x$parent, x$idx)
  x$parent2 <- x$idx2[i]
  x$parent2[is.na(x$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(x, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  x$name <- NA
  x$name[!x$split] <- tip.label

  phy <- reorder(structure(list(edge=cbind(x$parent2, x$idx2),
                                Nnode=Nnode,
                                tip.label=tip.label,
                                node.label=node.label,
                                edge.length=x$len,
                                orig=x),
                           class="phylo"))
  phy
}


## Adapted from prune.extinct.taxa() in geiger
prune <- function(phy) {
  to.drop <- subset(phy$orig, !split)$extinct
  if ( all(to.drop) ) {
    NULL
  } else if ( any(to.drop) ) {
    phy2 <- drop.tip(phy, phy$tip.label[to.drop])
    ## phy2$orig <- subset(phy2$orig, !extinct) # Check NOTE
    phy2$orig <- phy2$orig[!phy2$orig$extinct,]
    phy2$tip.state <- phy2$tip.state[!to.drop]
    phy2
  } else {
    phy
  }
}

