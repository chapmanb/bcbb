all: clean chapman_bcbio.bib
	pdflatex chapman_bcbio.tex
	bibtex chapman_bcbio
	pdflatex chapman_bcbio.tex
	pdflatex chapman_bcbio.tex

chapman_bcbio.bib:
	python fix_bibfile.py bcbio-nextgen.bib chapman_bcbio.bib

clean:
	rm -f chapman_bcbio.aux
	rm -f chapman_bcbio.bbl
	rm -f chapman_bcbio.blg
	rm -f chapman_bcbio.log

