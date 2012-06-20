## Internal time machine tests.  At some point I may directly expose
## the time machine, when this will become more important.  However,
## it is used directly in a number of places.

test.time.machine.bd <- function() {
  make.time.machine <- diversitree:::make.time.machine
  functions <- c(lambda="linear.t", mu="constant.t")
  t.max <- 10

  tm <- make.time.machine(functions, c(0, t.max))

  p1.0 <- c(.1, 0, .05)
  tm$set(p1.0)
  checkIdentical(tm$get(0),  p1.0[c(1,3)])
  checkIdentical(tm$get(10), p1.0[c(1,3)])

  ## Now, add a slope
  p1.1 <- c(.1, .05, .05)
  tm$set(p1.1)
  checkIdentical(tm$get(0),  p1.1[c(1,3)])
  checkIdentical(tm$get(10), c(p1.1[1] + p1.1[2]*10, p1.1[3]))

  ## Try a spline:  
  functions <- c(lambda="spline.t", mu="constant.t")
  x <- seq(0, t.max, length.out=101)
  y <- sin(x/t.max*2*pi)
  spline.data <- list(t=x, y=sin(x/t.max*2*pi))

  tm <- make.time.machine(functions, c(0, t.max),
                          spline.data=spline.data)
  p2.0 <- c(.1, .1, .05)
  p2.1 <- c(.1, .2, .05)
  tm$set(p2.0)

  checkIdentical(tm$get(0),  p2.0[c(1,3)])
  checkIdentical(tm$get(10), p2.0[c(1,3)])

  tm$set(p2.1)
  tt <- seq(0, 10, length=101)
  tmp <- tm$getv(tt)
  checkTrue(all(sapply(tmp[,2], identical, p2.0[3])))
  checkEquals(tmp[,1],
              p2.1[1] + diff(p2.1[1:2])*(sin(tt/t.max*2*pi)+1)/2)

  functions[1] <- "spline.linear.t"
  tm <- make.time.machine(functions, c(0, t.max),
                          spline.data=spline.data)
  p2.0 <- c(.1, .1, 0, .05)
  p2.1 <- c(.1, .2, .02, .05)
  tm$set(p2.0)

  tt <- seq(0, 10, length=101)
  tmp <- tm$getv(tt)

  checkIdentical(tm$get(0),  p2.0[c(1,4)])
  checkIdentical(tm$get(10), p2.0[c(1,4)])

  tm$set(p2.1)
  tt <- seq(0, 10, length=101)
  tmp <- tm$getv(tt)
  checkTrue(all(sapply(tmp[,2], identical, p2.1[4])))
  checkEquals(tmp[,1],
              p2.1[1] + diff(p2.1[1:2])*(sin(tt/t.max*2*pi)+1)/2
              + p2.1[3] * tt)
}

## MuSSE is more complicated, as there are Q matrices to deal with:
test.time.machine.musse <- function() {
  make.time.machine <- diversitree:::make.time.machine
  t.max <- 10  
  ## MuSSE version, to test out the Q matrix.  K of 3.
  k <- 3
  functions <- c(rep("linear.t", k),
                 rep("constant.t", k),
                 rep("constant.t", k*(k-1)))
  names(functions) <- diversitree:::default.argnames.musse(3)

  tm <- make.time.machine(functions, c(0, t.max), k=3)
  
  f <- diversitree:::make.pars.musse(k)

  p3.0 <- c(.1, 0, .2, 0, .3, 0,
            .01,   .02,   .03,
            .4, .5, .6, .7, .8, .9)
  tm$set(p3.0)

  checkIdentical(tm$get(0),  f(p3.0[-c(2,4,6)]))
  checkIdentical(tm$get(10), f(p3.0[-c(2,4,6)]))

  p3.1 <- c(.1, .05, .2, .05, .3, .05,
            .01,     .02,     .03,
            .4, .5, .6, .7, .8, .9)
  ## tm$set(p3.1)

  ##checkIdentical(tm$get(0),  f(p3.1[-c(2,4,6)]))
  ##checkIdentical(tm$get(10),
  ##               f(c(p3.1[c(1,3,5)] + p3.1[c(2,4,6)]*10, p3.1[-(1:6)])))

  ## Now, time-variable Q (this was not possible previously)
  functions <- c(rep("constant.t", k),
                 rep("constant.t", k),
                 rep("linear.t", k*(k-1)))
  names(functions) <- diversitree:::default.argnames.musse(3)
  tm <- make.time.machine(functions, c(0, t.max), k=3)

  p4.0 <- c(.1,  .2,  .3,
            .01, .02, .03,
            .4, 0, .5, 0, .6, 0, .7, 0, .8, 0, .9, 0)
  tm$set(p4.0)

  i <- seq(8, length(p4.0), by=2)
  checkIdentical(tm$get(0),  f(p4.0[-i]))
  checkIdentical(tm$get(10), f(p4.0[-i]))

  p4.1 <- c(.1,  .2,  .3,
            .01, .02, .03,
            .4, 0.02, .5, 0.02, .6, 0.02, .7, 0.02, .8, 0.02, .9, 0.02)
  tm$set(p4.1)

  checkIdentical(tm$get(0),  f(p4.1[-i]))
  checkIdentical(tm$get(10),
                 f(c(p4.1[c(1:6)], p4.1[i-1] + 10*p4.1[i])))

  ## Test negative parameter check:
  p4.2 <- p4.1
  p4.2[7] <- -.1
  checkException(tm$set(p4.2), silent=TRUE)
  p4.2 <- p4.1
  p4.2[1] <- -.1
  checkException(tm$set(p4.2), silent=TRUE)  
  p4.2 <- p4.1
  p4.2[8] <- -.1
  checkException(tm$set(p4.2), silent=TRUE)  
}
