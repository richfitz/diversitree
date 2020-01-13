all: compile_dll

compile_dll:
	Rscript -e "pkgbuild::compile_dll()"

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

check:
	_R_CHECK_CRAN_INCOMING_=FALSE make check_all

check_all:
	Rscript -e "rcmdcheck::rcmdcheck(args = c('--as-cran', '--no-manual'))"

test:
	./run_tests.R

run_examples: install
	make -C inst/examples

.PHONY: attributes document roxygen staticdocs publish_pages install clean build check
