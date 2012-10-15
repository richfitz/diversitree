library(diversitree)
library(multicore)
source("musse-simulations-fun.R")

dir.create("trees", FALSE)
dir.create("output", FALSE)

n.taxa <- c(50, 100, 200, 400) # logarithmic series
n.trait <- c(1, 2, 3, 4)       # 1..4 should be sufficient.
nn <- expand.grid(n.taxa=n.taxa, n.trait=n.trait)

## There are some integration errors here, which are probably driven
## by the non-strict state treatment.  This doesn't take too long to
## run, and can easily be run on one (multicore) computer.
## fits <- mapply(run.ml, nn$n.trait, nn$n.taxa, cores=6,
##                SIMPLIFY=FALSE)
## However, I don't explore this here.

if ( file.exists("output/obj.Rdata") ) {
  load("output/obj.Rdata")
} else {
  ## This takes considerably longer to run.  I ran these on our cluster,
  ## with different runs on different nodes
  ##   ans <- run.mcmc(1,  50)
  ##   ans <- run.mcmc(1, 100)
  ##   ...
  ##   ans <- run.mcmc(4, 400)
  ## and copied the results (about 1GB of them) back to my computer.
  ## This will load the returned results though.  Loading the results
  ## takes a bit of time, as they are quite large.
  ## 
  ## Loading this takes a while (and a lot of memory) as there are about
  ## 1GB of files to process.
  samples <- mapply(run.mcmc, nn$n.trait, nn$n.taxa, SIMPLIFY=FALSE)

  ## Process the samples to get means and 95% credibility intervals for
  ## each speciation main effect.
  f <- function(x)
    sapply(x[grep("lambda[A-Z]", names(x))],
           function(x) c(mean(x), hdr.uniroot(x)))

  obj <- lapply(samples, lapply, f)
  save(obj, file="output/obj.Rdata")
  rm(samples)
  gc()# force garbage collection before next step.
}

if ( file.exists("output/obj.inc.Rdata") ) {
  load("output/obj.inc.Rdata")
} else {
  ## The "incorrect" cases, dropping state 'A' were run in a similar way
  ## to the mcmc chains above:
  ##   ans <- run.mcmc(2,  50, incorrect=TRUE)
  ##   ans <- run.mcmc(2, 100, incorrect=TRUE)
  ##   ...
  ##   ans <- run.mcmc(4, 400, incorrect=TRUE)
  ##
  ## Loading also takes a while, for about 400MB of files.
  nn.inc <- subset(nn, n.trait > 1)
  samples.inc <- mapply(run.mcmc, nn.inc$n.trait, nn.inc$n.taxa,
                        incorrect=TRUE, SIMPLIFY=FALSE)

  ## The first lot are mislablled, but computed correctly:
  for ( i in which(nn.inc$n.trait == 2) )
    for ( j in seq_along(samples.inc[[i]]) )
      names(samples.inc[[c(i,j)]]) <-
        sub("A", "B", names(samples.inc[[c(i,j)]]))

  obj.inc <- lapply(samples.inc, lapply, f)
  save(obj.inc, file="output/obj.inc.Rdata")  
  
  rm(samples.inc, i, j)
  gc()
}

## fig.musse.mcmc(nn, obj)
## fig.power(nn, obj, obj.inc)

to.pdf("musse-multitrait-ci.pdf", 6, 6,
       fig.musse.mcmc(nn, obj))
to.pdf("musse-multitrait-power.pdf", 6, 6,
       fig.power(nn, obj, obj.inc))
