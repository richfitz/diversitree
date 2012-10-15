all: primates.pdf primates-incl.tex

primates.pdf: primates.tex
	pdflatex primates
	bibtex primates
	pdflatex primates
	pdflatex primates

primates.tex: primates.Rnw
	mkdir -p cache
	Sweave.sh -t --cache=cache primates.Rnw
	rm -f Rplots.pdf

primates-incl.tex: primates.tex
	awk /BEGIN/,/END/ primates.tex > primates-incl.tex

tidy:
	rm -f *.aux *.bbl *.blg *.log *.log *.out *.toc Rplots.pdf

clean: tidy
	rm -f cache/primates-*pdf
	rm -f primates.pdf primates.tex

.PHONY: all tidy clean
