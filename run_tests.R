#!/usr/bin/env Rscript
library(testthat, quietly=TRUE)
suppressMessages(library(geiger, quietly=TRUE))
suppressMessages(library(diversitree, quietly=TRUE))
library(parallel)
library(ggplot2)
test_dir("inst/tests")
