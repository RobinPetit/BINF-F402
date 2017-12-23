#!/bin/bash

fasta_files=($(ls *.fasta))
array_size=$(echo "${#fasta_files[@]}-1" | bc)

for i in $(seq 0 $array_size); do
	for j in $(seq 0 $array_size); do
		paf_file=${fasta_files[$i]%.fasta}_${fasta_files[$j]%.fasta}.paf
		eps_file=figs/${paf_file%.paf}.eps
		if [ ! -f $paf_file ]; then
			if [ $i -eq $j ]; then
				minimap2 -xasm5 -X ${fasta_files[$i]} ${fasta_files[$j]} > $paf_file
			else
				minimap2 -xasm5    ${fasta_files[$i]} ${fasta_files[$j]} > $paf_file
			fi
			echo -e "\n"  # makes logs more readble
		fi
		if [ ! -f $eps_file ]; then
			minidot -f20 -w1200 -d $paf_file > $eps_file
		fi
	done
done
