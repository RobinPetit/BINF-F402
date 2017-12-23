all: report.pdf
	evince $<

report.pdf: report.tex eps
	pdflatex $<
	pdflatex $<

eps:
	./make_all.sh

clean:
	rm -f report.aux report.log report.pdf *.paf figs/*

.PHONY: all eps clean
