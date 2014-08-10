#!/bin/bash

if [ $# -eq 0 ]; then
    echo "usage: fasta_flattener.sh fasta_file"
    exit 1
fi

sed 's/>.*/>/' $1 | tr -d '\n' 