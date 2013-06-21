#!/bin/sh

BASE=$1

cd $BASE
rm -f $BASE.Rmd $BASE.md
rm -rf cache
rm -rf figure
