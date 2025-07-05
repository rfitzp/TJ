# Makefile for Documentation

include ../Machines.in

# List of main .tex files (without .tex extension)
TEXFILES := $(wildcard *.tex)
TARGETS  := $(basename $(TEXFILES)) 

# Default target: build all PDFs
all: $(TARGETS:=.pdf)

# Rule to build a single PDF from a .tex file
%.pdf: %.tex
	$(L_COMPILER) $<
	$(L_COMPILER) $<
	$(L_COMPILER) $<

# Clean auxiliary files
clean:
	rm -f *.aux *.bib *.log *.bbl *.blg *.out *.toc *.lof *.lot *.synctex.gz

vclean:
	rm -f *.pdf *.aux *.bib *.log *.bbl *.blg *.out *.toc *.lof *.lot *.synctex.gz



