all: diversitree-tutorial.pdf diversitree-ext.pdf

diversitree-tutorial.tex: diversitree-tutorial.Rnw
	Sweave.sh -t --cache=cache diversitree-tutorial.Rnw
	rm -f Rplots.pdf

diversitree-tutorial.pdf: diversitree-tutorial.tex
	pdflatex diversitree-tutorial
	bibtex diversitree-tutorial
	pdflatex diversitree-tutorial
	pdflatex diversitree-tutorial

tidy:
	rm -f *.aux *.bbl *.blg *.log *.log *.out *.toc Rplots.pdf
	rm -f Makevars
	rm -f *.so *.o

clean: tidy
	rm -f cache/diversitree-tutorial-*pdf
	rm -f diversitree-tutorial.pdf diversitree-tutorial.tex
	rm -f cache/diversitree-ext-*pdf
	rm -f diversitree-ext.pdf diversitree-ext.tex

diversitree-ext.tex: diversitree-ext.Rnw
	Sweave.sh -t --cache=cache diversitree-ext.Rnw
	rm -f Rplots.pdf

diversitree-ext.pdf: diversitree-ext.tex
	pdflatex diversitree-ext
	pdflatex diversitree-ext

www: 
	tar -zcvf data.tar.gz data
	cp -p data.tar.gz \
		diversitree-tutorial.Rnw diversitree-tutorial.pdf \
		diversitree-ext.Rnw diversitree-ext.pdf \
		../www/doc
	cp -p *.c refs.bib refstyle.bst ../www/doc/files