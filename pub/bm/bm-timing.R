## Run this file by sourcing it, ot doing
##   Rscript bm-timing.R
## 
## This caches the timing results (which takes about ) in some .Rdata files.
library(diversitree)

source("bm-timing-fun.R")

n.taxa <- 2^(4:12)

obj <- run.cached("timings.Rdata",
                  timings(n.taxa, 10))
obj.err <- run.cached("timings-err.Rdata",
                  timings(n.taxa, 10, states.sd=1e-7))

times <- sapply(obj, sapply, function(x) mean(x[,"t.each"]))
times.err <- sapply(obj.err, sapply, function(x) mean(x[,"t.each"]))

to.pdf("bm-timing.pdf", 5, 5,
       fig.timing(n.taxa, times, times.err))

if ( FALSE ) { # don't run by default
  ## Here is the estimated exponent of qudratic growth:
  i <- (length(n.taxa)-1):(length(n.taxa))
  diff(log(times[i,]))/diff(log(n.taxa[i]))
  diff(log(times.err[i,]))/diff(log(n.taxa[i]))

  ## Timing ratios:
  times[9,1]     / times[9,]
  times.err[9,1] / times.err[9,]
}
