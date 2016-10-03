PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')

all: compile_dll

compile_dll:
	Rscript -e "devtools::compile_dll()"

clean:
	rm -f src/*.o src/*.so src/*.dll

attributes:
	Rscript -e "Rcpp::compileAttributes()"

staticdocs:
	@mkdir -p inst/staticdocs
	Rscript -e "library(methods); staticdocs::build_site()"

publish_pages:
	cd inst && ./update-gh-pages.sh

install:
	R CMD INSTALL .

build:
	R CMD build .

check: build
	R CMD check --as-cran --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

test:
	./run_tests.R

run_examples: install
	make -C inst/examples

.PHONY: attributes document roxygen staticdocs publish_pages install clean build check
