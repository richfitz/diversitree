#!/bin/sh

BASE=$1

cd $BASE
sowsear.sh $1.R

pandoc $BASE.md -o index.html --standalone --highlight-style=tango \
    -c ../stylesheet.css
