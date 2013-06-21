#!/bin/sh

## TODO: This should add some templating; the end-of-file part should
## be automatically added, perhaps.  A diversitree.footer() R function
## with 'asis' or something.  The pandoc call should include the
## google analytics fragment so I can see what people find
## interesting.

BASE=$1

cd $BASE
sowsear.sh $1.R

pandoc $BASE.md -o index.html --standalone --highlight-style=tango \
    -c ../stylesheet.css
