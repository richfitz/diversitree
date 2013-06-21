#!/bin/sh
## cp -p logo/logo.* logo/beetles.eps www
make -C www/bib all install
rsync -aLvz -e ssh --delete www/ zoology.ubc.ca:www/diversitree \
    --exclude 'examples/*/cache'
## Avoid future confusion about which file to add into!
rm -f www/studies.bib
