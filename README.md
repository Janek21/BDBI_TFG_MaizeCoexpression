# PGSB_CoexpressionNW

_From the earth I rise and to the earth I shall return, alas while I still draw breath I will be eating corn_
---

⚠️Repository in developement⚠️
Not all scripts are complete
---
This repository contains the files and data for executing a network co-expression analysis on maize lines B73, DK105, EP1, F7 and PE75.
It has programs for both of the the analysis done in my [internship](./programs/IndividualNetworks.Rmd) and in my bachelor thesis().

The original [expression data](./data/original_data) has been kept, but there is also the processed and split data (containing the gene [length](./data/wlen/) and without the gene [length](./data/nolen)), the mercator annotation data can also be found, as well as the data belonging to _Zea Mays_ genes and the [protein](./data/annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa) and [canonical](./data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical.cds.fa) sequences of this same species.

---
## Data preprocessing

The programs for preprocessing the data can be found in the [`scripts`](./scripts) directory and are all executed through a [bash](./scripts/Preprocessing.sh) file, a R markdown [file](./scripts/PreliminaryAnalysis.Rmd) containing the program for an inital data analysis is also found in this directory.

## Data analysis

## Annotation

etc
