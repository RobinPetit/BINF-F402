all: human.eps x_human.eps

hg38.fa: compile_hg38.sh download_hg38.sh
	rm -f hg38.fa
	./download_hg38.sh
	./compile_hg38.sh 5

x_human.eps: x_human.paf
	minidot -f15 -w1800 $< > $@

human.eps: human.paf
	minidot -f15 -w1800 $< > $@

human.paf: hg38.fa
	minimap2 -xasm5    $< $< > $@
	
x_human.paf: hg38.fa
	minimap2 -xasm5 -X $< $< > $@
	
