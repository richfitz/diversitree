# diversitree: comparative phylogenetic analyses of diversification

[![Build Status](https://travis-ci.org/richfitz/diversitree.png?branch=master)](https://travis-ci.org/richfitz/diversitree)

This repository contains "diversitree".  It is my experimental sources
though, and may not compile or work for you.  You may prefer the
released version from CRAN:

    > install.packages("diversitree")


The interesting directories are:

* inst/tests: testing functions that exercise most of the
package's main features.  Running `./run_tests.R` will run the tests.  These take too long to run on CRAN (over a minute), so are not set up in the usual way.
* doc: Vignettes, and their required data files.

## Installing diversitree

Clone the repository with

    git clone git://github.com/richfitz/diversitree.git

The package should then be installable the usual way.  You'll need a C, C++ and Fortran compiler.

To install and specify the location of the fftw library in a
non-standard place, a line like this is required:
  R CMD INSTALL diversitree --configure-args='--with-fftw=/path/to/fftw'
where the path will be the path so that the files
  /path/to/fftw/include/fftw3.h
  /path/to/fftw/lib/lib/fftw3.so
exist.

On Windows, set the environment variable LIB_FFTW to point to the
directory that contains include/fftw.h, and install the package.

If fftw is not found, installation will continue, but the (relatively)
fast C based QuaSSE integration will not be available.  The R based
fft integrator and the method-of-lines integrator will be available.

## Branches

The "master" branch contains the bleeding edge version of diversitree.
It may or may not work for you.  The "release" branch contains stable
releases.
