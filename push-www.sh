#!/bin/sh
cp -p logo/logo.* logo/beetles.eps www
cd www/bib; make all install
cd -
rsync --delete --exclude .hg -e ssh -rvzL www/ zoology:www/diversitree
