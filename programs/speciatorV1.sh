#!/bin/bash

cd /home/janek/Desktop/Intership/work/data

file=$1
metadata=sample.tissue.correct.cluadj.txt

species=$(cut -d$'\t' -f2 $metadata |sort |uniq)

#all genes have Zm (R all(grepl("Zm", rownames(data))) done in transposed data
R --vanilla --slave <<EOF
data<-read.delim("$1", row.names=1, stringsAsFactors=TRUE)
transposed <- t(data)
write.table(transposed, "T_$1", row.names=TRUE, sep="\t", eol="\n", col.names = NA)
EOF

for sp in $species; do
	echo $sp
	grep "Zm" T_$1 > $sp.csv 
	grep $sp T_$1>> $sp.csv
	echo "Data done"
	grep $sp $metadata>> "$sp"_m.csv
	echo "Metadata done"
	
done

#cut -d',' -f1 --complement B73.csv > h.csv
