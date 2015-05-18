dev_pdf = pdf
dev_tikz = tikz

ifeq ($(TRAVIS),true)
  dev=$(dev_pdf)
else
  dev=$(dev_tikz)
endif

all: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib clean-partial
	-xelatex  -interaction=nonstopmode "\input" $<
	-bibtex impatiens_phylogeography
	-xelatex  -interaction=nonstopmode "\input" $<
	xelatex  -interaction=nonstopmode "\input" $<

impatiens_phylogeography.tex: impatiens_phylogeography.Rnw code/impatiens_analysis.R
	Rscript -e "library(remake); Sys.setenv(DEV_TYPE = '$(dev)'); make('$@');"

impatiens_phylogeography_nourl.bib: impatiens_phylogeography.bib
	-cp ~/Library/impatiens_phylogeography.bib .
	Rscript parseURLs.R

impatiens_phylogeography.docx: impatiens_phylogeography.tex impatiens_phylogeography_nourl.bib
	pandoc -t docx -o $@ --csl systematic-biology.csl --bibliography impatiens_phylogeography_nourl.bib $<


clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~
	-rm impatiens_phylogeography_nourl.bib

clean: clean-partial
	-rm impatiens_phylogeography.pdf
	-rm impatiens_phylogeography.tex
