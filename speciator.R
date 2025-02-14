#!/usr/bin/env Rscript

#### Functions

specieCreator<-function(data_joint, spName){ #selects rows of a dataset based on specie column, then splits between data and metadata
  
  cat("\nSeparating by ", spName)
  data_specie<-data_joint[data_joint$specie == spName, ] #select rows according to metadata$specie column (now in joint),
  #data_specie is now a single specie dataset  
  
  #Separate back to data and metadata
  cat("\nSplitting data and metadata")
  geneCols<-grepl("Zm", names(data_specie)) #all gene columns(data) have Zm, so we search for them
  
  data<-data_specie[, geneCols] #data will be gene columns
  metadata<-data_specie[, !geneCols] #metadata will be all that is not gene columns
  
  cat("\n",spName, "specie table created\n")
  return(list(data, metadata))
}


listSplitter<-function(specieName, ALLspecieList){ #takes specie name and list of specie lists and makes data and metadata object for each specie
  
  cat("\nSplitting list for ", specieName,"\n")
  dataName<-paste0(specieName, "_data") #B73_data
  mdataName<-paste0(specieName, "_metadata") #B73_metadata
  specieData<-ALLspecieList[[specieName]][[1]] #Takes B73 list(produced by specieCreator) from MainList, then takes first item in B73 list (belongs to data)
  specieMetadata<-ALLspecieList[[specieName]][[2]] #same but 2nd slot in inner list, so metadata
  
  assign(paste0(dataName), specieData, envir = .GlobalEnv) #assigns 1rst inner list to B73_data name created
  assign(paste0(mdataName), specieMetadata, envir = .GlobalEnv) #same but 2nd and metadata
  
  cat("Data and metadata dataframes for ", specieName, " are done")
  return(NULL)
}


tableWriter<-function(specieName, 
                      data,
                      metadata, 
                      resultDataPath=paste0(specieName,".csv"), 
                      resultMetaPath=paste0(specieName,"_m.txt")){ #write the data and metadata tables
  
  View(metadata)
  cat("\nWriting", data,"data table to:", resultDataPath)
  write.table(data, resultDataPath, row.names=TRUE, sep="\t", eol="\n", col.names = NA)
  
  cat("\nWriting", metadata,"metadata table to:", resultMetaPath, "\n")
  write.table(metadata, resultMetaPath, row.names=TRUE, sep="\t", eol="\n", col.names = NA)
  
  return (NULL)
}

################################################################################


#### Dataset collection and preparation

#dataPath<- args[1]
#metadataPath<- args[2]
dataPath<-"../data/data_nolen.csv"
metadataPath<-"../data/metadata_3.txt"

dataNL<-read.delim(dataPath, row.names=1, stringsAsFactors=TRUE)
metadata<-read.delim(metadataPath, header=T, row.names=1, stringsAsFactors=TRUE)
colnames(metadata)<-c("specie", "quality", "tissue_abv", "rep", "location")

dataNL<-as.data.frame(t(dataNL))

datALL<-merge(dataNL, metadata, by="row.names")
rownames(datALL)<-datALL$Row.names
datALL$Row.names<-NULL

################################################################################


#### Program execution

specieList<-levels(metadata$specie) #get all specie names
#result is B73, a2, etc...
cat("The species are:\n", specieList)

specieResults <- lapply(specieList, function(specieName) specieCreator(datALL, specieName)) #creates a Main list containing each specie list produced by specieCreator
names(specieResults) <- specieList #names the lists according to specie (if not it would be 1,2,3,etc)

tmp<-lapply(specieList, function(specieName) listSplitter(specieName, specieResults)) #takes dataframes out of nested list and assigns specie name as variable

listRows<-sapply(specieResults, function(specieSublist) nrow(specieSublist[[1]])) #get number of rows of each specie, and put them on a vector

cat("\nNumber of rows of original dataframe:", nrow(datALL))
cat("\nTotal rows of all resulting databases combined:", sum(listRows))

View(paste0(specieName[[1]], "_data"))
################################################################################

#### Writing the datasets

#tmp<-lapply(names(specieResults), function(specieName) tableWriter(specieName, data=paste0(specieName, "_data"), metadata=paste0(specieName, "_metadata")))
