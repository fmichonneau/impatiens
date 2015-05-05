all: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib clean-partial
	-xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	-bibtex impatiens_phylogeography
	-xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex

impatiens_phylogeography_nourl.bib: impatiens_phylogeography.tex ~/Library/impatiens_phylogeography.bib
	cp ~/Library/impatiens_phylogeography.bib .
	Rscript parseURLs.R


clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~

clean: clean-partial
	-rm impatiens_phylogeography.pdf
	-rm impatiens_phylogeography.tex
