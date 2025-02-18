#!/bin/bash

cd /home/janek/Desktop/Intership/work/data

dataNL=original_data/all5gt.nolen.csv
metadata_original=original_data/sample.tissue.correct.cluadj.txt

#1 remove subpar quality samples from data (column3, keep the 0)
m1=metadata_1.txt

echo "Removing inferior quality samples from $metadata_original"

awk '$3 == 0' $metadata_original > $m1 #creates metadata_1 (m1)

echo "Writing results to $m1" #sample.tissue.correct.cluadj.txt no more use
echo ""

#2 match metadata and data samples, remove unneeded from metadata (metdata has pollen data present, which is nonexistent in data)
m2=metadata_2.txt

echo "Matching samples in $dataNL and $m1"

dataNL_new=data_nolen.csv

Rscript ../scripts/MetaFixer.R $dataNL $m1 $dataNL_new $m2 #create metadata_2 (m2)

echo ""

#3 tissue abreviation mapping
m3=metadata_3.txt

echo "Clarifying the tissue code column in $m2"

python3 ../scripts/tissue_abv.py $m2 $m3 #create m3

echo ""

#Last: Split table by species(follow metadata species guideline)

Rscript ../scripts/speciator.R $dataNL_new $m3 ./ ./ #destinations of data and metadata tables are in /home/janek/Desktop/Intership/work/data/

#mclapply USE ^^

#POst processing (in R result): join by tissue replicate? (1 for cob, 1 for etc)
