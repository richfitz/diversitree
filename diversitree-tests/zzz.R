#!/usr/bin/env Rscript
library(RUnit)
library(diversitree)

options(warn=1)

testSuite <- defineTestSuite(name="diversitree",
                             dirs=".",
                             testFileRegexp="runit.+\\.R$",
                             testFuncRegexp="^test.+",
                             rngKind="default",
                             rngNormalKind="default")

## testSuite <- defineTestSuite(name="diversitree",
##                              dirs=".",
##                              testFileRegexp="runit-quasse.R$",
##                              testFuncRegexp="^test.+",
##                              rngKind="default",
##                              rngNormalKind="default")

res <- runTestSuite(testSuite)
printTextProtocol(res)
