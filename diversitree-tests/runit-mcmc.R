ignore <- function() {
  library(diversitree)
  library(RUnit)
  library(lubridate)

  library(diversitree)
  make.mvn <- function(mean, vcv) {
    logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
    tmp <- length(mean) * log(2 * pi) + logdet
    vcv.i <- solve(vcv)
    
    function(x) {
      dx <- x - mean
      -(tmp + rowSums((dx %*% vcv.i) * dx))/2
    }
  }

  vcv <- matrix(c(1, .25, .25, .75), 2, 2)
  lik <- make.mvn(c(0, 0), vcv)
  set.seed(1)
  samples1 <- mcmc(lik, c(0, 0), 100, 1, print.every=10)

  samples2 <- mcmc(lik, NULL, 200, 1, print.every=10,
                   previous=samples1)

  ## samples2 <- mcmc(lik, NULL, 200, 1, print.every=10,
  ##                  previous=samples1)

  set.seed(1)
  samples3 <- mcmc(lik, c(0, 0), 200, 1, print.every=10)

  set.seed(1)
  samples3 <- mcmc(lik, c(0, 0), 200, 1, print.every=10)

  identical(samples2, samples3)

  ## Then, with saving:
  samples1 <- mcmc(lik, c(0, 0), 10000, 1, print.every=10,
                   save.every=10, save.file="test.rds")

  library(lubridate)
  samples1 <- mcmc(lik, c(0, 0), 10000, 1, print.every=10,
                   save.every.dt=seconds(4), save.file="test.rds")
}
