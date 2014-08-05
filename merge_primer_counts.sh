#!/bin/bash
#
# Joins two primer files specified by the first two arguments. The
# first file is the foreground primers file, the second is the
# background primers file. 
# 
# The output (to stdout) are the common primers, followed by the
# background binding count and then the foreground binding count.
# 
# Erik Clarke - ecl@mail.med.upenn.edu

usage="merge_primer_counts.sh foreground_primers background_primers"

if [[ -f "$1" && -f "$2" ]]; then
    awk 'FNR==NR{a[$1]=$2 FS $3;next}{OFS = "\t"; print $0, a[$1]}' $1 $2 | awk 'NF > 2' | awk '{$1=$1}1'
else 
    echo $usage
    echo "Error: no vaid foreground/background primers specified."
    exit 1
fi




