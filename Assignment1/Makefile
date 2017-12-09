NB_PLOTS=20
LAST_PLOT_IDX=$(shell expr $(NB_PLOTS) - 1)
INDICES=$(shell seq 0 $(LAST_PLOT_IDX))

plots: Rplots.pdf
	for idx in $(INDICES); do \
	    echo Extracting $$idx; \
	    convert -density 500 Rplots.pdf[$$idx] plots/plot$$idx.pdf; \
	    convert -density 500 plots/plot$$idx.pdf plots/plot$$idx.png; \
	done

Rplots.pdf: assignment1.R
	./assignment.R
