diversitree: comparative phylogenetic analyses of diversification
=================================================================

This repository contains "diversitree".  It is my experimental sources
though, and may not compile or work for you.  You may prefer the
released version from CRAN:

    > install.packages("diversitree")


The interesting directories are:

* diversitree: the R package
* diversitree-tests: testing functions that exercise most of the
package's main features (run `Rscript zzz.R` in this directory to run
all tests).
* doc: Vignettes, and their required data files.

Installing diversitree
----------------------

Clone the repository with

    git clone git://github.com/richfitz/diversitree.git


You need to generate the configure script.  In the
`diversitree/diversitree` directory, run

    autoheader
    autoconf

The package should then be installable the usual way.

Branches
--------

The "master" branch contains the bleeding edge version of diversitree.
It may or may not work for you.  The "release" branch contains stable
releases, which you can usually get from the diversitree website.

