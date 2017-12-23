#!/bin/bash

download() {
	if [ ! -f chr$1.fa ]; then
		wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr$1.fa.gz
		gunzip chr$1.fa.gz
		rm -f chr$1.fa.gz
	fi
}

for i in $(echo {1..22} X Y); do
	download $i &
done
