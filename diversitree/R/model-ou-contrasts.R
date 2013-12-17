## OK; I think that this one is going to be the general case (CF
## special case for BM and possibly also for lambda)
##
##   1. rescale the tree
##   2. compute contrasts
##   3. compute likelihood
##
## The actual likelihood calculation part can probably be factored
## out, I suspect.  It will be shared between the bm and ou cases
## only.

## Most of what is done below should be done in the make.cache.ou,
## really.  But it's here so that I can get this basically working.
make.all.branches.ou.contrasts <- function(cache, control) {
  if (any(cache$states.sd > 0))
    stop("Cannot (yet) do contrasts based bm models with state error")

  ## This is annoying; we'll have to do this in make.bm.  If we don't,
  ## then we will reorder every time and that's going to hurt,
  ## timewise.  If I rewrite the pic() calculations to use my
  ## ordering, we can skip this step.  Given that eventually I think I
  ## want access to more of the pic components, that seems like a good
  ## idea.
  tree <- reorder(cache$info$phy, "pruningwise")
  states <- cache$states[tree$tip.label] # does reorder change this?
  rescale <- make.rescale.phylo.ou(tree)
  n <- length(tree$tip.label)

  function(pars, intermediates, preset=NULL) {
    s2 <- pars[[1]]
    alpha <- pars[[2]]
    tree <- rescale(alpha)

    ## Copied from model-bm.R:make.cache.bm
    ## There is duplication here with make.cache.pgls; perhaps merge
    ## these?  That might help with some of the root treatment
    ## things.
    pics <- pic(cache$states, tree, var.contrasts=TRUE)
    u <- pics[,"contrasts"]
    V <- pics[,"variance"]
    V0 <- pgls.root.var.bm(tree)
    ## This step is brutal.  We could rewrite the branch lengths there
    ## once the lookups are done properly.  This is particularly bad
    ## because it will involve doing another tree reordering.
    ## Costly!
    root.x <- pgls.root.mean.bm(tree, cache$states)

    ## This is the bit that is shared with all.branches.bm
    ll <- -(n * log(2 * pi * s2) +
            sum(log(V)) +
            log(V0) +
            sum(u * u) / s2) * 0.5
    list(loglik=ll,
         root.x=root.x,
         root.v=s2 * V0,
         # not sure if these are needed...
         V=V, V0=V0)
  }
}
