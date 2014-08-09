#!/bin/bash

usage="fasta2oneline.sh input_fasta \n
Removes newlines from the sequences in a fasta file."

header_lines=$(grep -n '>' $1 | cut -d: -f1)

for line in header_lines; do
    sed -n "$(line){p;q}" $1
    

#echo $header
#tail -n +2 < $1 | tr -d '\n'
#echo
