#!/bin/bash

HG_FILE=hg38.fa

chr=($(echo {1..22} X Y))

if [ $# -eq 2 ]; then
	output=$2
else
	output=${HG_FILE}
fi

rm -f $output

for idx in ${!chr[@]}; do
	if [ $idx -eq $1 ]; then
		break
	fi
	cat chr${chr[$idx]}.fa >> $output
done
