⚠️Repository in developement⚠️
# PGSB_CoexpressionNW

_From the earth I rise and to the earth I shall return, alas while I still draw breath I will be eating corn_
---
This repository contains the files and data for executing a network co-expression analysis on maize lines B73, DK105, EP1, F7 and PE75.
It has programs for both of the the analysis done in my [internship](./programs/IndividualNetworks.Rmd) and in my bachelor thesis().

The original [expression data](./data/original_data) has been kept, but there is also the processed and split data (containing the gene [length](./data/wlen/) and without the gene [length](./data/nolen)), the mercator annotation data can also be found, as well as the data belonging to _Zea Mays_ genes and the [protein](./data/annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa) and [canonical](./data/annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical.cds.fa) sequences of this same species.

---
## [`Scripts`](./scripts) directory

[`Preliminary analysis`](./scripts/PreliminaryAnalysis.Rmd): R markdown for initial viewing of the data

[`Bash execution file`](./scripts/Preprocessing.sh): Executes all the steps for the initial data cleaning

- [`MetaFixer.R`](./scripts/MetaFixer.R): Matches the amount of samples to the lowest one among data and metadata
- [`tissue_abv.py`](./scripts/tissue_abv.py): Assigns to each metadata tissue code a tissue abreviation for a better comprehension.
- [`speciator.R`](./scripts/speciator.R): Splits the data and metadata tables by lines, necessary only for the internship project.

## [`Programs`](./programs) directory

## Annotation

etc
