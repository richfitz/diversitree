#!/bin/sh
## cp -p logo/logo.* logo/beetles.eps www
make -C www/bib all install
rsync --delete -e ssh -rvzL www/ zoology.ubc.ca:www/diversitree
