all: diversitree-tutorial.pdf

diversitree-tutorial.pdf: diversitree-tutorial.tex
	pdflatex diversitree-tutorial
	bibtex diversitree-tutorial
	pdflatex diversitree-tutorial
	pdflatex diversitree-tutorial

diversitree-tutorial.tex: diversitree-tutorial.Rnw
	R CMD Sweave diversitree-tutorial.Rnw
	rm -f Rplots.pdf

tidy:
	rm -f *.aux *.bbl *.blg *.log *.log *.out *.toc Rplots.pdf

clean: tidy
	rm -f cache/diversitree-tutorial-*pdf
	rm -f diversitree-tutorial.pdf diversitree-tutorial.tex
