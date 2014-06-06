SOURCES = $(shell ls *.cpp *.c *.f)

TARGETS = diversitree.so

# Have to set this up manually, otherwise we can't find Rcpp or BH's
# headers.
RCPP_FLAG = $(shell Rscript -e "Rcpp:::CxxFlags()")
export PKG_CPPFLAGS = $(RCPP_FLAG)

all: $(TARGETS)

diversitree.so: $(SOURCES)
	R CMD SHLIB -o $@ $(SOURCES)

clean:
	rm -f *.o *.so

.PHONY: all clean
