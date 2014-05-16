knit: impatiens_knit.R impatiens_phylogeography.Rnw
	Rscript impatiens_knit.R

impatiens_phylogeography.aux: impatiens_phylogeography.tex
	cp ~/Library/impatiens_phylogeography.bib .
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	bibtex impatiens_phylogeography.aux

all: knit impatiens_phylogeography.aux impatiens_phylogeography.tex impatiens_phylogeography.bib
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	make clean-partial

clean-partial:
	-rm impatiens_phylogeography.tex
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *.out

clean: clean-partial
	-rm *.pd
