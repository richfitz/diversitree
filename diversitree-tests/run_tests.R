# usage: R --slave < run_tests.R

library(RUnit)
library(diversitree)

test.suite <- defineTestSuite("geosse_classe", dirs=getwd(), rngKind=RNGkind()[1])
# Defaults to requiring "test." in function names and "runit" in file names
# Careful! The RUnit default random number generator is different than the R default.  
#   This way, runs by hand and within runit give the same results.
 
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result)
