all: BiGGR.Rnw
	$(R_HOME)/bin/R   CMD Sweave BiGGR.Rnw	
	$(R_HOME)/bin/R   CMD pdflatex BiGGR.tex
	$(R_HOME)/bin/R   CMD bibtex BiGGR	
	$(R_HOME)/bin/R   CMD pdflatex BiGGR.tex
	$(R_HOME)/bin/R   CMD pdflatex BiGGR.tex
	make clean

clean:
	rm -fv *.tex *.aux *.dvi *.log *.toc *.bak *~ *.blg *.bbl *.lot *.lof
	rm -fv *.nav *.snm *.out \#*\# _region_* _tmp.* *.vrb *.lim
	rm -fv Rplots.pdf *.RData