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
