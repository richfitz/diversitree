sim.character <- function(tree, pars, x0=0, model="bm", br=NULL) {
  if ( is.null(br) ) {
    model <- match.arg(model, c("bm", "ou", "bbm", "mk", "meristic"))
    br <- match.fun(sprintf("make.br.%s", model))(pars)
  }
  ## TODO: Add checking for br().

  edge <- tree$edge
  idx <- seq_len(max(edge))
  n.tip <- length(tree$tip.label)
  root <- n.tip + 1
  is.tip <- idx <= n.tip
  children <- lapply(idx[!is.tip], function(x)
                     edge[edge[,1] == x, 2])
  children <- rbind(matrix(NA, n.tip, 2),
                    t(matrix(unlist(children), 2)))
  order <- rev(get.ordering(children, is.tip, root))
  len <- tree$edge.length[match(idx, edge[, 2])]

  y <- rep.int(NA, length(len))
  y[root] <- x0

  for ( i in order ) {
    j <- children[i,]
    y[j] <- br(y[i], len[j])
  }

  ret <- y[is.tip]
  names(ret) <- tree$tip.label
  ret
}

## BM: one parameter -- s2
make.br.bm <- function(pars) {
  if ( pars <= 0 )
    stop("s2 must be positive")
  function(x0, t)
    rnorm(length(t), x0, sqrt(pars*t))
}

## OU: three parameters -- s2, alpha, theta
make.br.bm <- function(pars) {
  s2    <- pars[1]
  alpha <- pars[2]
  theta <- pars[3]
  if ( s2 <= 0 )
    stop("s2 must be positive")
  if ( alpha <= 0 )
    stop("alpha must be positive")

  function(x0, t)  
    rnorm(length(t),
          theta + (x0 - theta)*exp(-alpha*t),
          sqrt(s2/(2*alpha) * (1 - exp(-2*alpha*t))))
}

## BBM: three parameters (drift not yet allowed!)
make.br.bbm <- function(pars) {
  s2 <- pars[1]
  c <- pars[2]
  d <- pars[3]

  if ( s2 <= 0 )
    stop("s2 must be positive")
  if ( c > d )
    stop("c must be smaller than d")

  ## TODO: x0 within [c,d], rewrite to be scalar.
  function(x0, t) {
    x <- rnorm(length(t), x0, sqrt(s2*t))
    while ( any(x < c | x > d) ) {
      if ( any(i <- x < c) )
        x[i] <- 2*c - x[i]
      if ( any(i <- x > d) )
        x[i] <- 2*d - x[i]
    }
    x
  }
}

## Three parameters 
make.br.meristic <- function(pars) {
  k <- pars[1]
  up <- pars[2]
  down <- pars[3]
  
  Q <- matrix(0, k, k)
  idx <- cbind(1:(k-1), 2:k)
  Q[idx] <- up
  Q[idx[,2:1]] <- down
  diag(Q) <- -rowSums(Q)

  make.br.mk(Q)
}

make.br.mk <- function(pars) {
  if ( is.matrix(pars) ) {
    Q <- pars
    k <- nrow(Q)
    if ( k == 0 || k != ncol(Q) )
      stop("'pars' must be a square matrix")
  }
  if ( is.vector(pars) ) {
    k <- sqrt(length(pars))
    Q <- matrix(pars, k, k)
  }
  r <- -diag(Q)
  function(x0, t) {
    if ( length(t) != 1 )
      stop("t must be a scalar.")
    while ( r[x0] > 0 && (t <- t - rexp(1, r[x0])) > 0 )
      x0 <- sample((1:k)[-x0], 1, prob=Q[x0,-x0])
    x0
  }
}
