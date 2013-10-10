\name{trait.plot}
\alias{trait.plot}

\title{Plot a Phylogeny and Traits}

\description{Plot a phylogeny and label the tips with traits.
  This function is experimental, and may change soon.  Currently it can
  handle discrete-valued traits and two basic tree shapes.}

\usage{
trait.plot(tree, dat, cols, lab=names(cols), str=NULL, 
           class=NULL, type="f", w=1/50, 
           legend=length(cols) > 1, cex.lab=.5, 
           font.lab=3, cex.legend=.75, margin=1/4, 
           check=TRUE, quiet=FALSE, ...)
}

\arguments{
  \item{tree}{Phylogenetic tree, in ape format.}

  \item{dat}{A \code{data.frame} of trait values.  The row names must be
    the same names as the tree (\code{tree$tip.label}), and each column
    contains the states (0, 1, etc., or \code{NA}).  The column names
    must give the trait names.}

  \item{cols}{A list with colors.  Each element corresponds to a trait
    and must be named so that all names appear in \code{names(dat)}.
    Each of these elements is a vector of colors, with length matching
    the number of states for that trait.  Traits will be plotted in the
    order given by \code{cols}.}

  \item{lab}{Alternative names for the legend (perhaps longer or more
    informative).  Must be in the same order as \code{cols}.}

  \item{str}{Strings used for the states in the legend.  If \code{NULL}
    (the default), the values in \code{dat} are used.}

  \item{class}{A vector along \code{phy$tip.label} giving a higher level
    classification (e.g., genus or family).  No checking is done to
    ensure that such classifications are not polyphyletic.}

  \item{type}{Plot type (same as \code{type} in \code{?plot.phylo}).
    Currently only \code{f} (fan) and \code{p} (rightwards phylogram)
    are implemented.}

  \item{w}{Width of the trait plot, as a fraction of the tree depth.}

  \item{legend}{Logical: should a legend be plotted?}

  \item{cex.lab, font.lab}{Font size and type for the tip labels.}

  \item{cex.legend}{Font size for the legend.}

  \item{margin}{How much space, relative to the total tree depth, should
    be reserved when plotting a higher level classification.}

  \item{check}{When TRUE (by default), this will check that the classes
    specified by \code{class} are monophyletic.  If not, classes will be 
    concatenated and a warning raised.}

  \item{quiet}{When TRUE (FALSE by default), this suppresses the warning
    caused by \code{check=TRUE}.}

  \item{...}{Additional arguments passed through to phylogeny plotting
    code (similar to \code{ape}'s \code{plot.phylo}).}
}


\examples{
## These are the parameters: they are a single speciation and extinction
## rate, then 0->1 (trait A), 1->0 (A), 0->1 (B) and 1->0 (B).
colnames(musse.multitrait.translate(2, depth=0))

## Simulate a tree where trait A changes slowly and B changes rapidly.
set.seed(1)
phy <- tree.musse.multitrait(c(.1, 0, .01, .01, .05, .05),
                             n.trait=2, depth=0, max.taxa=100,
                             x0=c(0,0))
## Here is the matrix of tip states (each row is a species, each column
## is a trait).
head(phy$tip.state)

trait.plot(phy, phy$tip.state,
           cols=list(A=c("pink", "red"), B=c("lightblue", "blue")))

nodes <- c("nd5", "nd4", "nd7", "nd11", "nd10", "nd8")
grp <- lapply(nodes, get.descendants, phy, tips.only=TRUE)
class <- rep(NA, 100)
for ( i in seq_along(grp) )
  class[grp[[i]]] <- paste("group", LETTERS[i])

## Now, 'class' is a vector along phy$tip.label indicating which of six
## groups each species belongs.

## Plotting the phylogeny with these groups:
trait.plot(phy, phy$tip.state,
           cols=list(A=c("pink", "red"), B=c("lightblue", "blue")),
           class=class, font=1, cex.lab=1, cex.legend=1)

## Add another state, showing values 1:3, names, and trait ordering.
tmp <- sim.character(phy, c(-.1, .05, .05, .05, -.1, .05, .05, 0.05, -.1), 
                     model="mkn", x0=1)
phy$tip.state <- data.frame(phy$tip.state, C=tmp)
trait.plot(phy, phy$tip.state,
           cols=list(C=c("palegreen", "green3", "darkgreen"),
                     A=c("pink", "red"), B=c("lightblue", "blue")),
           lab=c("Animal", "Vegetable", "Mineral"),
           str=list(c("crane", "toad", "snail"), c("kale", "carrot"), 
                    c("calcite", "beryl")))

## Rectangular/phylogram plot with groups.
trait.plot(ladderize(phy, right=FALSE), phy$tip.state, type="p",
           cols=list(A=c("pink", "red"), B=c("lightblue", "blue"),
                     C=c("palegreen", "green3", "darkgreen")),
           class=class, font=1, cex.lab=1)
}

\keyword{hplot}
\author{Richard G. FitzJohn}
