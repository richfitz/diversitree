# usage: R --slave < run_tests.R

library(RUnit)
library(diversitree)

# Defaults to requiring "test." in function names and "runit" in file names
# The full set of tests is run through zzz.R.  This file is just handy for
# testing one model at a time.

#test.suite <- defineTestSuite("geosse_classe", dirs=getwd(), rngKind=RNGkind()[1])
test.suite <- defineTestSuite("geosse_classe", dirs=".", testFileRegexp="runit-classe.R$", rngKind=RNGkind()[1])

# Careful! The RUnit default random number generator is different than the R
# default.  This way, runs by hand and within runit give the same results.
 
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
