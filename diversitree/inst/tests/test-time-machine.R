## OK, so there is an old verison of time machine that hopefully
## excercises most of what we want.

library(diversitree)
library(testthat)

context("Time machine")

make.time.machine <- diversitree:::make.time.machine
functions <- c(lambda="linear.t", mu="constant.t")
t.max <- 10

tm <- make.time.machine(functions, c(0, t.max))

p1.0 <- c(.1, 0, .05)
tm$set(p1.0)
expect_that(tm$get( 0), is_identical_to(p1.0[c(1,3)]))
expect_that(tm$get(10), is_identical_to(p1.0[c(1,3)]))

tm2 <- diversitree:::make.time.machine2(functions, c(0, t.max))

tm2$set(p1.0)
expect_that(tm2$get( 0), is_identical_to(p1.0[c(1,3)]))
expect_that(tm2$get(10), is_identical_to(p1.0[c(1,3)]))

## Now, add a slope
p1.1 <- c(.1, .05, .05)
tm$set(p1.1)

expect_that(tm$get(0),  is_identical_to(p1.1[c(1,3)]))
expect_that(tm$get(10), is_identical_to(c(p1.1[1] + p1.1[2]*10, p1.1[3])))

tm2$set(p1.1)
expect_that(tm2$get(0),  is_identical_to(p1.1[c(1,3)]))
expect_that(tm2$get(10), is_identical_to(c(p1.1[1] + p1.1[2]*10, p1.1[3])))

## Here I should consider splines...

## MuSSE.
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
## Indices of the slope and intercept parameters here for use in
## testing code.
i.c <- c(1, 3, 5)
i.m <- c(2, 4, 6)

tm$set(p3.0)

expect_that(tm$get( 0), is_identical_to(f(p3.0[-i.m])))
expect_that(tm$get(10), is_identical_to(f(p3.0[-i.m])))

tm2 <- diversitree:::make.time.machine2(functions, c(0, t.max), k=3)
tm2$set(p3.0)

expect_that(tm2$get( 0), is_identical_to(f(p3.0[-i.m])))
expect_that(tm2$get(10), is_identical_to(f(p3.0[-i.m])))

## With a slope
p3.1 <- c(.1, .05, .2, .05, .3, .05,
          .01,     .02,     .03,
          .4, .5, .6, .7, .8, .9)

ff <- function(t)
  f(c(p3.1[i.c] + p3.1[i.m] * t,
      p3.1[-c(i.c, i.m)]))

tm$set(p3.1)
tm2$set(p3.1)

expect_that(tm$get(0), is_identical_to(f(p3.0[-i.m])))
expect_that(tm$get(0), is_identical_to(ff(0)))
expect_that(tm$get(10), is_identical_to(ff(10)))

expect_that(tm2$get(0), is_identical_to(f(p3.0[-i.m])))
expect_that(tm2$get(0), is_identical_to(ff(0)))
expect_that(tm2$get(10), is_identical_to(ff(10)))

## Now, time-variable Q
functions <- c(rep("constant.t", k),
               rep("constant.t", k),
               rep("linear.t", k*(k-1)))
names(functions) <- diversitree:::default.argnames.musse(3)
tm <- make.time.machine(functions, c(0, t.max), k=3)
tm2 <- diversitree:::make.time.machine2(functions, c(0, t.max), k=3)

p4.0 <- c(.1,  .2,  .3,
          .01, .02, .03,
          .4, 0, .5, 0, .6, 0, .7, 0, .8, 0, .9, 0)
tm$set(p4.0)
tm2$set(p4.0)

i <- seq(8, length(p4.0), by=2)
expect_that(tm$get( 0), is_identical_to(f(p4.0[-i])))
expect_that(tm$get(10), is_identical_to(f(p4.0[-i])))

expect_that(tm2$get( 0), is_identical_to(f(p4.0[-i])))
expect_that(tm2$get(10), is_identical_to(f(p4.0[-i])))

p4.1 <- c(.1,  .2,  .3,
          .01, .02, .03,
          .4, 0.02, .5, 0.02, .6, 0.02, .7, 0.02, .8, 0.02, .9, 0.02)
tm$set(p4.1)
tm2$set(p4.1)

expect_that(tm$get( 0), is_identical_to(f(p4.1[-i])))
expect_that(tm$get(10),
            is_identical_to(f(c(p4.1[c(1:6)], p4.1[i-1] + 10*p4.1[i]))))

expect_that(tm2$get(0), is_identical_to(f(p4.1[-i])))
expect_that(tm2$get(10),
            is_identical_to(f(c(p4.1[c(1:6)], p4.1[i-1] + 10*p4.1[i]))))

## Negative parameter check -- constant parameter:
p4.2 <- p4.1
p4.2[1] <- -.1
expect_that(tm$set(p4.2), throws_error())
expect_that(tm2$set(p4.2), throws_error())

## Negative parameter check -- variable parameter:
p4.2 <- p4.1
p4.2[7] <- -.1
expect_that(tm$set(p4.2), throws_error())

## Throws later here.
tm2$set(p4.2)
expect_that(tm2$get(0), throws_error())

## Then, do same with the flattened version
tm2 <- diversitree:::make.time.machine2(functions, c(0, t.max), k=3,
                                        truncate=TRUE)

## Constant parameter:
p4.2 <- p4.1
p4.2[1] <- -.1
tm2$set(p4.2)
## truncated parameter:
expect_that(tm2$get(0)[1], is_identical_to(0))

## Variable parameter (input 7 ends up in position 10)
p4.2 <- p4.1
p4.2[7] <- -.1
tm2$set(p4.2)
expect_that(tm2$get(0)[10], is_identical_to(0))

## And the negative-OK version
tm2 <- diversitree:::make.time.machine2(functions, c(0, t.max), k=3,
                                        nonnegative=FALSE)

## Constant parameter:
p4.2 <- p4.1
p4.2[1] <- -.1
tm2$set(p4.2)
## truncated parameter:
expect_that(tm2$get(0)[1], is_identical_to(-.1))

## Variable parameter (input 7 ends up in position 10)
p4.2 <- p4.1
p4.2[7] <- -.1
tm2$set(p4.2)
expect_that(tm2$get(0)[10], is_identical_to(-.1))
