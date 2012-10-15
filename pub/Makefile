all: diversitree.pdf diversitree-suppl.pdf

diversitree.pdf: diversitree.tex example/primates-incl.tex
	pdflatex diversitree
	bibtex diversitree
	pdflatex diversitree
	pdflatex diversitree

diversitree-suppl.pdf: diversitree-suppl.tex
	pdflatex diversitree-suppl
	bibtex diversitree-suppl
	pdflatex diversitree-suppl
	pdflatex diversitree-suppl

bm/bm-timing.pdf:
	cd bm && Rscript bm-timing.R

example/primates-incl.tex:
	cd example && make

simulations/musse-multitrait-power.pdf:
	cd simulations && Rscript musse-simulations.R


