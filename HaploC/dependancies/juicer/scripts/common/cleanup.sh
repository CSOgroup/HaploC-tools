#!/bin/bash
#########
# Script to clean up big repetitive files and zip fastqs. Run after you are 
# sure the pipeline ran successfully.  Run from top directory (HIC001 e.g.).
# Juicer version 1.5
total=`ls -l aligned/merged_sort.txt | awk '{print $5}'`
total2=`ls -l aligned/merged_nodups.txt aligned/dups.txt aligned/opt_dups.txt | awk '{sum = sum + $5}END{print sum}'`
if [ $total -eq $total2 ] 
then 
    rm aligned/merged_sort.txt 
    rm -r splits 
    for i in fastq/*.fastq
    do
        gzip $i
    done
    gzip aligned/merged_nodups.txt
    gzip aligned/dups.txt
    gzip aligned/opt_dups.txt
    gzip aligned/abnormal.sam
    gzip aligned/collisions.txt
    gzip aligned/unmapped.sam
else 
    echo "Problem: The sum of merged_nodups and the dups files is not the same size as merged_sort.txt"
    echo "Did NOT clean up";
fi
