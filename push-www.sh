#!/bin/sh
rsync --delete --exclude .hg -e ssh -rvzL www/ zoology:www/diversitree
