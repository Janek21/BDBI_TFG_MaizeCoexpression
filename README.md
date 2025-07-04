
# Maize Co-Expression Network Analysis

---
This repository contains the files and data for executing a network co-expression analysis on maize lines B73, DK105, EP1, F7 and PE75.
It has programs for both: the the analysis for [individual lines](./coexpression_code/IndividualNetworks.Rmd) and for [multiple lines](./coexpression_code/WholeLineNetworks.Rmd) simultaneously.

---
## [`Data`](./data) directory:

[`original_data`](./data/original_data): Original expression data

[`Length data`](./data/wlen): Data tables that contain the length(used for RPKM), one table with joint data for all the lines and many of data split by the lines ([scripts](#Scripts-directory))

[`Regular data`](./data/nolen): Contains a table with the data for all lines and tables where the data is split([scripts](#Scripts-directory)), one table per line.

[`Annotation`](./data/annotation): Contains files with gene annotation for _Zea Mays_.

## [`Scripts`](./scripts) directory

Contains programs for the initial data analysis and pre-processing

[`Preliminary analysis`](./scripts/PreliminaryAnalysis.Rmd): R markdown for initial viewing of the data

[`Bash execution file`](./scripts/Preprocessing.sh): Executes all the steps for the initial data cleaning

* [`MetaFixer.R`](./scripts/MetaFixer.R): Matches the amount of samples to the lowest one among data and metadata
* [`tissue_abv.py`](./scripts/tissue_abv.py): Assigns to each metadata tissue code a tissue abreviation for a better comprehension.
* [`speciator.R`](./scripts/speciator.R): Splits the data and metadata tables by lines, necessary only for separate analysis of the lines.

## [`Co-expression`](./coexpression_code) directory

Contains programs and results belonging to the co-expression network creation and analysis

[`IndividualNetworks.Rmd`](./coexpression_code/IndividualNetworks.Rmd): Co-expression network creation and correaltion analysis for individual lines.

[`WholeLineNetworks.Rmd`](./coexpression_code/WholeLineNetworks.Rmd): Co-expression network creation and correaltion analysis simultaneous for all lines.

[`geneModule.txt`](./coexpression_code/geneModule.txt): Table containing each gene and the module to which it belongs

[`Distribution plots`](./coexpression_code/DistrPlots): Contains the data distributions for the raw, normalized(through CPM and RPKM), filtered and unfiltered data.

[`Expression plots`](./coexpression_code/ExprPlots): Where PCA and UMAP data visualizations are saved

[`Correlation plots`](./coexpression_code/CorPlots): Contains the plots and tables where the module co-expression across tissues (and lines) is reflected.


## [`Annotation`](./annotation_code) directory

Contains programs and results which belong to functional annotation of both, modules and genes.

### [`Module functions`](./annotation_code/modules): Contains the code and results for obtaining the main functions of modules

* [`FileRead_class.py`](./annotation_code/modules/FileRead_class.py): Classes for parsing files (mercator annotation [file](./data/annotationb73.mercator.v4.7.txt) and P-value [tables](./annotation_code/modules/Pres)).
* [`FunctionPval_tables.py`](./annotation_code/modules/FunctionPval_tables.py): Calculates the significance for functions in each module, saves them in [tables](./annotation_code/modules/Pres).
* [`SignificantPerModule_plots.R`](./annotation_code/modules/SignificantPerModule_plots.R): Plots the significance per function in each module and the top significant functions, saves them in a [directory](./annotation_code/modules/SigPlots).
* [`TopSignificantFunctions.py`](./annotation_code/modules/TopSignificantFunctions.py): Writes a [table](./annotation_code/modules/database.txt) of most significant functions in each module for each annotation type.
* [`Pres`](./annotation_code/modules/Pres): Tables containing function, significance and module for each annotation type, created with [FileRead_class.py](./annotation_code/modules/FileRead_class.py).
* [`SigPlots`](./annotation_code/modules/SigPlots): Contains plots and tables of significance per function in each module and of the top significant functions per annotation type, which are created using [SignificantPerModule_plots.R](./annotation_code/modules/SignificantPerModule_plots.R).


### [`Gene functions`](./annotation_code/genes): Contains the code and results for obtaining the gene communities and annotations.

  [`Annotation tables`](./annotation_code/genes/community_tables): Tables containing gene annotations for most of the genes in each module, as well as a [script](./annotation_code/genes/community_tables/ComRatio.sh) that counts known and unknown genes for each community in a module.

  [`Graph representations`](./annotation_code/genes/graphFigs): Image representations of module correlations and of communities within modules.

* [`maize.B73.AGPv4.aggregate.gaf`](./annotation_code/genes/maize.B73.AGPv4.aggregate.gaf): File containing GO term, taxon, evidence code and more data belonging to _Zea Mays_ genes.
* [`go-basic.obo`](./annotation_code/genes/go-basic.obo): Contains GO term information.
* [`geneClustering.Rmd`](./annotation_code/genes/geneClustering.Rmd): Performs a K-means and DBscan clustering for a given module, also calculates and saves the [correlation table](./annotation_code/correlation_tables) belonging to this module.
* [`GOparser_class.py`](./annotation_code/genes/GOparser_class.py): Class for parsing the maize [GAF file](./annotation_code/genes/maize.B73.AGPv4.aggregate.gaf), also contains a function to translate from v4 IDs to v5 IDs based on a [translation file](./data/annotation/genes_all.txt)
* [`GOannotation.py`](./annotation_code/genes/GOannotation.py): Finds unknown genes and n closely related genes, the GO terms and functions of the known genes and extrapolates the functions of the unknown genes from these ones.
* [`comm_annotation.py`](./annotation_code/genes/comm_annotation.py): Finds Louvain communities from a given [correlation table](./annotation_code/correlation_tables) and annotates known and unknown genes in this communities using [GOannotation.py](./annotation_code/genes/GOannotation.py). Then creates a [table](./annotation_code/genes/community_tables) which contains the annotations.

## [`Figures`](./all_figure) directory

Contains all the folders that hold figures or tables produced as a result by a program

## Contact
For any questions or further information, contact me at: jan.izquierdo@alum.esci.upf.edu
