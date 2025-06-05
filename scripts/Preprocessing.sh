#!/bin/bash

cd /home/jj/Desktop/Bioinformatics/3rd_year/Internship/HelmholtzZentrum_Munchen/Network_WGCNA/gitSent/data

option="${1:-nolen}" #unless input is given, default option will be nolen


data=original_data/all5gt.$option.csv
metadata_original=original_data/sample.tissue.correct.cluadj.txt

echo "Working on $data"

#1 remove subpar quality samples from data (column3, keep the 0)
m1=metadata_1.txt

echo "Removing inferior quality samples from $metadata_original"

awk '$3 == 0' $metadata_original > $m1 #creates metadata_1 (m1)

echo "Writing results to $m1" #sample.tissue.correct.cluadj.txt no more use
echo ""

#2 match metadata and data samples, remove unneeded from metadata (metdata has pollen data present, which is nonexistent in data)
m2=metadata_2.txt

echo "Matching samples in $data and $m1"

data_new=data_$option.csv

Rscript ../scripts/MetaFixer.R $data $m1 $data_new $m2 #create metadata_2 (m2)

echo ""

#3 tissue abreviation mapping
m3=metadata_3.txt

echo "Clarifying the tissue code column in $m2"

python3 ../scripts/tissue_abv.py $m2 $m3 #create m3

echo ""

#Last: Split table by lines(follow metadata lines guideline)

Rscript ../scripts/speciator.R $data_new $m3 ./ ./

echo ""

#Rename metadata, remove old metadata
mv $m3 metadata.out #.out is so that its not recognized in the nex step
rm metadata_*.txt -v

#Move out all species files to their corresponding length option folder
mkdir $option -v #create nolen or wlen folder

mv *.csv $option -v
mv *.txt $option -v #shift all species data into the new folder

mv metadata.out metadata.txt
