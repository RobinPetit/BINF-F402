NB_PLOTS=20
LAST_PLOT_IDX=$(shell expr $(NB_PLOTS) - 1)
INDICES=$(shell seq 0 $(LAST_PLOT_IDX))
DENSITY=500

plots: Rplots.pdf
	for idx in $(INDICES); do \
	    echo Extracting $$idx; \
	    convert -density $(DENSITY) Rplots.pdf[$$idx] plots/plot$$idx.eps; \
		done
	    #convert -density $(DENSITY) plots/plot$$idx.eps plots/plot$$idx.png; done

Rplots.pdf: assignment1.R
	./assignment1.R
