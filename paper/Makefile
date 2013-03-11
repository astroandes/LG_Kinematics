LATEX = latex
MAIN  = ms


main: 	$(MAIN).tex  
	$(LATEX) $(MAIN).tex 
	bibtex $(MAIN)
	$(LATEX) $(MAIN).tex 
	$(LATEX) $(MAIN).tex 
	dvips $(MAIN).dvi -o $(MAIN).ps
	ps2pdf $(MAIN).ps


clean:
	rm -f $(MAIN).{aux,bbl,ps,pdf,div,blg,log}
	rm -f *~

view: 
	open $(MAIN).pdf &

