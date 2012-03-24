all: diversitree-tutorial.pdf

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

diversitree-ext.tex: diversitree-ext.Rnw
	Sweave.sh -t --cache=cache diversitree-ext.Rnw
	rm -f Rplots.pdf

diversitree-ext.pdf: diversitree-ext.tex
	pdflatex diversitree-ext
	pdflatex diversitree-ext



# all: diversitree-tutorial.pdf

# diversitree-tutorial.pdf: diversitree-tutorial.tex
# 	pdflatex diversitree-tutorial
# 	bibtex diversitree-tutorial
# 	pdflatex diversitree-tutorial
# 	pdflatex diversitree-tutorial

# diversitree-tutorial.tex: diversitree-tutorial.Rnw
# 	Sweave.sh -c diversitree-tutorial.Rnw
# 	rm -f Rplots.pdf

# diversitree-ext.pdf: diversitree-ext.tex
# 	pdflatex diversitree-ext
# #	bibtex diversitree-ext
# #	pdflatex diversitree-ext
# 	pdflatex diversitree-ext

# diversitree-ext.tex: diversitree-ext.Rnw
# 	Sweave.sh -c diversitree-ext.Rnw
# 	rm -f Rplots.pdf

# diversitree-ext2.pdf: diversitree-ext2.tex
# 	pdflatex diversitree-ext2
# 	pdflatex diversitree-ext2

# diversitree-ext2.tex: diversitree-ext2.Rnw
# 	Sweave.sh -c diversitree-ext2.Rnw
# 	rm -f Rplots.pdf

# tmp.tex: tmp.Rnw
# 	Sweave.sh -c tmp.Rnw
# tmp.pdf: tmp.tex
# 	pdflatex tmp

# tidy:
# 	rm -f *.aux *.bbl *.blg *.log *.log *.out *.toc Rplots.pdf

# clean: tidy
# 	rm -f cache/diversitree-tutorial-*pdf
# 	rm -f diversitree-tutorial.pdf diversitree-tutorial.tex
