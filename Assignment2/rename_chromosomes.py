#!/usr/bin/python3

from sys import argv

from os.path import exists
from os import system

from roman import toRoman

try:
    header = argv[1]
    path = argv[2]
    assert exists(argv[2])
    nb_chr = int(argv[3])
except:
    print('Usage: ./rename_chromosomes.py <header> <path> <# chromosomes>.\n\te.g.: ./rename_chromosomes.py "Kyo_" Kyokai7.fasta 16')
    exit()

for i in range(1, nb_chr+1):
    roman = toRoman(i)
    system('sed -i \'s/{0}{1}$/{0}{2}/\' {3}'.format(header, roman, i, path))
