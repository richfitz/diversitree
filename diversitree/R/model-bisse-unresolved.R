## Compute the probabilities of generating every state, given that we
## started with one species in state 0 or state 1.  For lt = length(t)
## times, returns a n(space) x lt x 2 array.
##   'scal': a multiplier to get use of extra precision (the vector
##           will be scaled so that the total probability sums to
##           scal.  Underflow is possible where scal=1, though).
##   'tol': tolerance, used as exit condition (when changes in the
##          probability space do not exceed tol, the step size is good
##          enough).  Lower values _may_ be faster and _may_ have
##          lower accuracy (but with large numbers of clades being
##          calculated, this effect can be small).  It is not clear to
##          me if tol should be resacled by scal before being passed
##          in.
##   'm': parameter affecting internal arguments.  Small numbers often
##        faster, sometimes slower.  Play around.
bucexp <- function(nt, mua, mub, laa, lab, qba, qab, t, scal=1,
                   tol=1e-10, m=15) {
  stopifnot(all(c(mua, mub, laa, lab, qba, qab) >= 0))
  lt <- length(t)
  stopifnot(lt > 0)
  stopifnot(scal >= 1)
  stopifnot(m > 1)
  n <- (nt*(nt+1)/2+1)
  res <- .Fortran("bucexp",
                  nt   = as.integer(nt),
                  mua  = as.numeric(mua),
                  mub  = as.numeric(mub),
                  laa  = as.numeric(laa),
                  lab  = as.numeric(lab),
                  qba  = as.numeric(qba),
                  qab  = as.numeric(qab),
                  t    = as.numeric(t),
                  lt   = as.integer(lt),
                  scal = as.numeric(scal),
                  tol  = as.numeric(tol),
                  m    = as.integer(m),
                  w    = numeric(2*n*lt),
                  flag = integer(1))
  if ( res$flag < 0 )
    stop("Failure in the Fortran code")
  else if ( res$flag > 0 )
    cat(sprintf("Flag: %d\n", res$flag))
  array(res$w, c(n, lt, 2))
}

## This is more useful; given model parameters, and vectors of times
## (t), number of species in each clade (Nc), numbers of species in
## known states (nsc) and numbers of species in state b (k), calculate
## the probabilities of generating the clade, and of extinction.
## Returns a length(t) x 4 matrix, where the first two columns are the
## probabilities of generating the clade, and the second two are the
## probabilities of extinction.
bucexpl <- function(nt, mua, mub, laa, lab, qba, qab, t, Nc, nsc, k,
                    scal=1, tol=1e-10, m=15) {
  tt <- sort(unique(t))
  ti <- as.integer(factor(t))

  lt <- length(tt)
  lc <- length(Nc)
  stopifnot(lc > 0, lt > 0, max(Nc) < nt,
            length(Nc) == lc, length(nsc) == lc, length(k) == lc,
            all(nsc <= Nc), all(nsc >= 0), all(k <= nsc),
            scal >= 1, m > 1)

  res <- .Fortran("bucexpl",
                  nt   = as.integer(nt),
                  mua  = as.numeric(mua),
                  mub  = as.numeric(mub),
                  laa  = as.numeric(laa),
                  lab  = as.numeric(lab),
                  qba  = as.numeric(qba),
                  qab  = as.numeric(qab),
                  t    = as.numeric(tt),
                  lt   = as.integer(lt),
                  ti   = as.integer(ti),
                  Nc   = as.integer(Nc),
                  nsc  = as.integer(nsc),
                  k    = as.integer(k),
                  lc   = as.integer(lc),
                  scal = as.numeric(scal),
                  tol  = as.numeric(tol),
                  m   = as.integer(m),
                  ans  = numeric(4*lc),
                  flag = integer(1))
  if ( res$flag < 0 )
    stop("Failure in the Fortran code")
  else if ( res$flag > 0 )
    cat(sprintf("Flag: %d\n", res$flag))
  matrix(res$ans, ncol=4)
}

## Compute the PDF of the hypergeometric distribution, using the same
## parameters that Wikipedia uses, and I used in the Fortran
## function.
hyperg <- function(N, m, n, k) dhyper(k, m, N-m, n)

## The function 'bucexp.n' creates a data.frame with the number of
## species in state a, b and total for a bucexp state-space absorbing
## at n species.
bucexp.n <- function(n) {
  z <- sapply(0:(n-1), seq, from=0)
  nb <- unlist(z)
  na <- unlist(lapply(z, rev))
  nt <- na + nb
  rbind(data.frame(na, nb, nt), c(NA, NA, n))
}

## Pack a nt x nt matrix with probabilities returned by bucexp().
## Cases where na + nb > nt are given zero probabilities (or change,
## with the 'default' argument; e.g., default=NA will set them to be
## NA).
repack <- function(p, default=0) {
  n <- (sqrt(8*length(p) - 7) - 1)/2
  m <- matrix(default, n, n)
  idx <- bucexp.n(n)
  m[with(idx[-nrow(idx),]+1, cbind(na, nb))] <- p[-length(p)]
  m
}

## Given two of na, nb and nt, calculate the position in the state
## vector.
index <- function(na, nb, n=na + nb) {
  if ( missing(nb) ) nb <- n - na
  if ( missing(na) ) na <- n - nb
  stopifnot(na + nb == n && all(na >= 0) && all(nb >= 0))
  n*(n + 1)/2 + 1 + nb
}

## Construct the transition matrix.  Again, non-R style, as this was
## used as a template for constructing the same in Fortran.
make.matrix <- function(nt, mu.a, mu.b, lambda.a, lambda.b, q.ba,
                        q.ab) {
  ## Diagonals (are linear sums of these: -(na*pa + nb*pb)):
  pa <- lambda.a + mu.a + q.ba
  pb <- lambda.b + mu.b + q.ab

  ## A few useful numbers:
  n1 <- nt*(nt-1)/2
  nm <- n1 + nt + 1

  ## Create some useful vectors
  j1 <- integer(n1)
  na <- integer(n1)
  nb <- integer(n1)
  k <- 1
  for ( i in 1:(nt-1) ) {
    for ( j in 1:i ) {
      j1[k] <- i + k
      na[k] <- i - j + 1
      nb[k] <- j
      k <- k + 1
    }
  }

  ## Character state transitions and extinctions
  m <- matrix(0, nm, nm)
  for ( i in 1:n1 ) {
    m[j1[i],   i]       <- na[i]*mu.a
    m[j1[i]+1, i]       <- nb[i]*mu.b
    m[j1[i],   j1[i]+1] <- na[i]*q.ba
    m[j1[i]+1, j1[i]]   <- nb[i]*q.ab
  }

  ## Speciation and diagonals:
  for ( i in 2:n1 ) {
    if ( na[i] > 1 )
      m[i, j1[i]]   <- (na[i]-1)*lambda.a
    if ( nb[i] > 1 )
    m[i, j1[i]+1] <- (nb[i]-1)*lambda.b
    m[i,i]        <- - (na[i]-1)*pa - (nb[i]-1)*pb
  }

  ## Speciation in the special final column, diagonals for the last
  ## class.
  k <- nt*(nt-1)/2
  for ( i in 1:nt ) {
    m[k+i, nm]  <-  (nt-i)*lambda.a + (i-1)*lambda.b
    m[k+i, k+i] <- -(nt-i)*pa       - (i-1)*pb
  }

  m
}
