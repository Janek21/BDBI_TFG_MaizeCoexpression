#!/usr/bin/env Rscript

#data matcher function
dataMatcher<-function(data, metadata){
  options(warn=-1)
  
  cat("Data begins with:" , dim(data))
  cat("\nMetadata begins with:", dim(metadata))
  #Match data to metadata
  data <- data[, order(colnames(data))]
  metadata <- metadata[order(rownames(metadata)), ]
  
  cat("\nColumns data = Rows of metatdata?", all(rownames(metadata) == colnames(data)))
  #If TRUE, columns of data and rows of metadata are matched
  
  cat("\nRemove the excess from data")
  data<-data[,colnames(data) %in% rownames(metadata)] #remove data not present in metadata
  cat("\nData end with:" , dim(data))
  
  cat("\nRemove the excess from metadata")
  metadata<-metadata[rownames(metadata) %in% colnames(data),] #remove metadata not present in data
  cat("\nMetadata end with:" , dim(metadata), "\n")
  
  options(warn=0)
  return(list(data, metadata))
}

#get args from bash
args = commandArgs(trailingOnly=TRUE)

dataPath<- args[1]
metadataPath<- args[2]
resultDataPath<- args[3]
resultMetaPath<- args[4]

data<-read.delim(dataPath, row.names=1, stringsAsFactors=TRUE)
metadata<-read.delim(metadataPath, header=FALSE, row.names=1, stringsAsFactors=TRUE)

#data<-read.delim("../data/original_data/all5gt.nolen.csv", row.names=1, stringsAsFactors=TRUE)
#metadata<-read.delim("../data/original_data/sample.tissue.correct.cluadj.txt", header=F, row.names=1, stringsAsFactors=TRUE)


#If dataNL is data with Length: create a row in metadata named Length as well, so that the 2 datasets can be joined, later we will remove it
if ("Length" %in% colnames(data)){
  metadata<-as.data.frame(t(metadata))
  metadata$Length<-rep(0, nrow(metadata))
  metadata<-as.data.frame(t(metadata))
  
}

jointData<-dataMatcher(data, metadata)

data_fixed<-jointData[[1]]
metadata_fixed<-jointData[[2]]

cat("\nWriting resulting data table to:", resultDataPath)
write.table(data_fixed, resultDataPath, row.names=TRUE, sep="\t", eol="\n", col.names = NA)

cat("\nWriting resulting metadata table to:", resultMetaPath, "\n")
write.table(metadata_fixed, resultMetaPath, row.names=TRUE, sep="\t", eol="\n", col.names = NA)