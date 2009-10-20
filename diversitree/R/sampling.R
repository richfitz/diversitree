sample.phylo <- function(phy, type=c("random", "clade"), p) {
}

## Random sample of a tree - this one is easy!
sample.phylo.random <- function(phy, p) {
  n.taxa <- length(phy$tip.label)
  n.drop <- round((1-p) * n.taxa)
  to.drop <- sample(n.taxa, n.drop)
  phy2 <- drop.tip(phy, to.drop)
  phy2$orig <- phy
  phy2
}

