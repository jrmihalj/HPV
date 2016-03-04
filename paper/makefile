# Makefile for the HPV article

all: main.pdf

main.pdf: main.tex refs.bib
	# make the pdf
	latex main
	bibtex main
	latex main
	latex main
	dvipdf main

	# clean up some intermediate files
	rm *.dvi
	rm *.blg
	rm *.aux
	
