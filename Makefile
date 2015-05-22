dev_pdf = pdf
dev_tikz = tikz

ifeq ($(TRAVIS),true)
  dev=$(dev_pdf)
else
  dev=$(dev_tikz)
endif

all: impatiens.tex impatiens_nourl.bib clean-partial
	-xelatex  -interaction=nonstopmode "\input" $<
	-bibtex impatiens
	-xelatex  -interaction=nonstopmode "\input" $<
	xelatex  -interaction=nonstopmode "\input" $<

impatiens.tex: impatiens.Rnw remake.yml
	-rm impatiens_nourl.bib
	Rscript -e "library(remake); Sys.setenv(DEV_TYPE = '$(dev)'); make('all');"

impatiens_nourl.bib:
	-cp ~/Library/impatiens.bib .
	Rscript parseURLs.R

impatiens.docx: impatiens.tex impatiens_nourl.bib
	pandoc -t docx -o $@ --csl systematic-biology.csl --bibliography impatiens_nourl.bib $<


clean-partial:
	-rm *.bbl
	-rm *.blg
	-rm *.aux
	-rm *.log
	-rm *~

clean: clean-partial
	-rm impatiens.pdf
	-rm impatiens.tex
