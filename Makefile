all: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib clean-partial
	-xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	-bibtex impatiens_phylogeography
	-xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex
	xelatex  -interaction=nonstopmode "\input" impatiens_phylogeography.tex

impatiens_phylogeography.tex: impatiens_phylogeography.Rnw code/impatiens_analysis.R
	Rscript -e "library(remake); make('$@');"

impatiens_phylogeography_nourl.bib: impatiens_phylogeography.tex ~/Library/impatiens_phylogeography.bib
	-cp ~/Library/impatiens_phylogeography.bib .
	Rscript parseURLs.R

impatiens_phylogeography.docx: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib
	# pandoc -t docx -o $@ --csl systematic-biology.csl --bibliography impatiens_phylogeography_nourl.bib $^
	pandoc -t docx -o $@ $^


clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~

clean: clean-partial
	-rm impatiens_phylogeography.pdf
	-rm impatiens_phylogeography.tex
