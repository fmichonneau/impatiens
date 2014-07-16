all: impatiens_phylogeography.aux impatiens_phylogeography.tex
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	make clean-partial

impatiens_phylogeography.tex: impatiens_phylogeography.Rnw code/impatiens_analysis.R code/impatiens_map.R
	Rscript -e "library(knitr); knit('impatiens_phylogeography.Rnw')"

impatiens_phylogeography_nourl.bib: impatiens_phylogeography.tex impatiens_phylogeography.bib ~/Library/impatiens_phylogeography.bib
	cp ~/Library/impatiens_phylogeography.bib .
	Rscript parseURLs.R

impatiens_phylogeography.aux: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	bibtex impatiens_phylogeography.aux

clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~

clean: clean-partial
	-rm *.pdf
	-rm impatiens_phylogeography.tex
