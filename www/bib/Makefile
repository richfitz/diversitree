export TMPDIR=.

all: studies.html

install: ../studies.html

../studies.html: studies.html
	cp -p studies.html studies.bib ..

studies.html: studies_tmp.html template.html
	python process.py

studies_tmp.html: studies.bib
	bibtex2html --style sysbio3 --no-keys -a --dl \
		 -nobibsource -nodoc  -nofooter -o studies_tmp studies.bib 

clean:
	rm -f studies_tmp.html

.PHONY: all clean install
