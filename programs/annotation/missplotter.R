library(ggplot2)
library(stringr)
library(dplyr)

#read data
database<-read.delim("~/Desktop/Intership/work/programs/annotation/Pres/database.txt")
pvdata<-read.table("./Pres/d_1_sheet.txt", sep='\t', header=TRUE)
genemodule<-read.table("../geneModule.txt", header=TRUE)
gFile<-read.delim("~/Desktop/Intership/work/data/annotation/b73.mercator.v4.7.txt", fill=TRUE, na.strings=c("''", "NA", "#N/A"))

#remove '
gFile[] <- lapply(gFile, function(x) gsub("'", "", x))

#Clean Unnanotated and empty rows
gFile<-gFile[!is.na(gFile$IDENTIFIER),]
geneFunction<-gFile[!grepl("not assigned", gFile$NAME),]  #["not assigned.not annotated"!=gFile$NAME,]
geneFunction<-geneFunction[!grepl("Enzyme", geneFunction$NAME),]

#Clean Gene IDs for later match with MEs
geneFunction$IDENTIFIER<-gsub("_[^_]*$", "", geneFunction$IDENTIFIER)
geneFunction$IDENTIFIER<-gsub("z", "Z", geneFunction$IDENTIFIER)

#Get functions at depth1
geneFunction[c('D1Function', 'Function')] <- str_split_fixed(geneFunction$NAME, '\\.', 2)

#associate function to module (function>gene>module)
genemodule$IDENTIFIER<-rownames(genemodule)
temp<-left_join(geneFunction, genemodule, by="IDENTIFIER")
geneFmodule<-data.frame("genes"=temp$IDENTIFIER, "D1Function"=temp$D1Function, "Function"=temp$Function, "Module"=temp$modules)

#match gene and function to pvalue (only 1 function/module combination > find pvalue)
genPval_table<-left_join(geneFmodule, pvdata, by = c("Module", "D1Function"="Function"))

#Set data tyoes
genPval_table$genes<-as.factor(genPval_table$genes)
genPval_table$D1Function<-as.factor(genPval_table$D1Function)
genPval_table$Function<-as.factor(genPval_table$Function)
genPval_table$Module<-as.factor(genPval_table$Module)
genPval_table$Relevancy<-as.logical(genPval_table$Relevancy)

#Get relevant functions only, remove NA
PvalTrue<-genPval_table[genPval_table$Relevancy,]
PvalTrue<-PvalTrue[!is.na(PvalTrue$D1Function),]

ggplot(PvalTrue, aes(x=genes, y=D1Function, alpha=Pvalue, color=Pvalue))+
  geom_point()+
  labs(x="Genes", y="Functions")+
  facet_wrap(PvalTrue$Module)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(size=10, color = "black"),
        legend.position="right",
        legend.key.width=unit(1, "cm"), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=16),
        plot.background = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = "grey"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.5,),
        strip.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.ticks=element_blank())

#Focus on a ME
ME<-"blue"
PvalFocus<-PvalTrue[PvalTrue$Module==ME,]

ggplot(PvalFocus, aes(x=genes, y=D1Function, color=Pvalue, size=Pvalue))+
  geom_point()+
  labs(x="Genes", y="Functions", title=ME)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y=element_text(size=10, color = "black"),
        title=element_text(size=20, face="bold"),
        legend.position="right",
        legend.key.width=unit(1, "cm"), 
        legend.text=element_text(size=14), 
        legend.title=element_text(size=16),
        plot.background = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = "grey"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        panel.grid.major = element_line(colour = "lightgrey", linetype = "solid", size = 0.5,),
        strip.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.ticks=element_blank())



##Use tissue instead of function? pvaluue of each function per tissue
dataPath<-paste0("../../data/wlen/data_wlen.csv")
metadataPath<-paste0("../../data/metadata.txt")
dataNL<-read.delim(dataPath, row.names=1, stringsAsFactors=TRUE)
metadata<-read.delim(metadataPath, header=T, row.names=1, stringsAsFactors=TRUE)

colnames(metadata)<-c("specie", "quality", "tissue_abv", "rep", "location")
length_vec<-dataNL$Length

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

jointData<-dataMatcher(dataNL, metadata)

dataNL<-jointData[[1]]
metadata<-jointData[[2]]

################
################
library(edgeR)

dfA<-metadata
dfB<-dataNL

#Total Counts per tissue
rsumer<-function(data, metadata, tissue_name){ #calculates the sum of all columns that belong to a location (ex: mean of all COB columns)
  
  loc_mdata<-metadata[metadata$location == tissue_name, ] #filter metadata tissue (get metadata of location only)
  
  data<-data[,colnames(data) %in% rownames(loc_mdata)] #get data of location only, based on metadata
  
  data<-rowSums(data) #calculate mean for each gene out of the locations(replicates)
  
  return(as.data.frame(data))
}

tissue_data<-levels(dfA$location) #get list of tissue names

d_joint<-sapply(tissue_data, function(tissue_name) rsumer(dfB, dfA, tissue_name)) #returns an array where each entry is a column with the mean data of the replicates (rows are genes)

repl_data<-as.data.frame(d_joint) #data joint by replicate

colnames(repl_data) = gsub(pattern = "*.data", replacement = "", x = tolower(colnames(repl_data))) #get column names to be only location

rownames(repl_data)<-rownames(dfB) #rename rows to be genes again

length_vec<-data.frame(Length=length_vec) #convert to dataframe
dge <- DGEList(repl_data,genes=length_vec) #use edgeR for normalization
dge <- calcNormFactors(dge)
Nrepl_data <- rpkm(dge, log=TRUE)
Nrepl_data<-as.data.frame(Nrepl_data)

dataTableMaker<-function(data, gene, tissue){
  cell<-data[gene, tissue]
  geneNT<-sum(data[gene,])
  NgeneT<-sum(data[,tissue])
  exclusiondf<-Nrepl_data[rownames(Nrepl_data)!=gene, colnames(Nrepl_data)!=tissue] #all rows but gene, all column but tissue
  NgeneNT<-sum(exclusiondf)
                                                        #format on the table
  dd<-matrix(c(cell, NgeneT, geneNT, NgeneNT), nrow=2)  # cell    geneNT
  return(dd)                                            # NgeneT  NgeneNT
}

results<-list()

for (tissue in colnames(Nrepl_data)){
  results[[tissue]]<-lapply(rownames(Nrepl_data), function(gene){ #use lapply to iterate ober all rows, get argument as gene and pass to function
    dataTableMaker(Nrepl_data, gene, tissue)})
}
saveRDS(file="results.RDS", results)
s=fisher.test(results[["cob"]])
print(s)
