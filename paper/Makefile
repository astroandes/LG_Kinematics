LATEX = pdflatex
MAIN  = ms


main: 	$(MAIN).tex  
	$(LATEX) $(MAIN).tex 
	bibtex $(MAIN)
	$(LATEX) $(MAIN).tex 
	$(LATEX) $(MAIN).tex 
#	dvips $(MAIN).dvi -o $(MAIN).ps



clean:
	rm -f $(MAIN)~
	rm -f $(MAIN).aux
	rm -f $(MAIN).bbl
	rm -f $(MAIN).ps
	rm -f $(MAIN).pdf
	rm -f $(MAIN).dvi
	rm -f $(MAIN).blg
	rm -f $(MAIN).log
	rm -f *~
	rm -f $(MAIN).tex.flc

view: 
	open $(MAIN).pdf &

