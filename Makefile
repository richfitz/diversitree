install:
	R CMD INSTALL .

clean:
	rm -f src/*.o src/diversitree.so

test:
	./run_tests.R

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr diversitree*gz | tail -n1`
	@rm -f `ls -1tr diversitree*gz | tail -n1`
	@rm -rf forest.Rcheck

.PHONY: all build install clean test
