TEXFILE= BentancorClarke_AdolescentFertility

# all runs All R-Scripts (which have subfunctions) and then compiles TeX
all :
	../Source/Births/BirthEstimates.Rout: ../Source/Births/BirthEstimates.R
		R CMD BATCH $<
	../Source/Deaths/DeathEstimates.Rout: ../Source/Deaths/DeathEstimates.R
		R CMD BATCH $<
	../Source/Summary/sumStats.Rout:      ../Source/Summary/sumStats.R
		R CMD BATCH $<
	../Source/placebo/placeboTest.Rout:   ../Source/placebo/placeboTest.R
		R CMD BATCH $<

	$(TEXFILE).pdf: $(TEXFILE).tex
		pdflatex -shell-escape $(TEXFILE)
		bibtex $(TEXFILE).aux
		pdflatex -shell-escape $(TEXFILE)
		pdflatex -shell-escape $(TEXFILE)


clean :
	rm -fv *.aux *.log *.blg *.bbl
	rm -fv $(TEXFILE).pdf

veryclean :
	rm -fv *.aux *.log *.blg *.bbl
	rm -fv $(TEXFILE).pdf
	rm -fv .Rout