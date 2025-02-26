## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(purl = TRUE)


## ----loadLib---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(WGCNA)
allowWGCNAThreads()
library(randomcoloR)
library(edgeR)
library(tidyverse)
library(dplyr)
library(gridExtra)
#devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
library(ggpubr)


## ----dtLoad----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

d<-read.delim("../data/metadata.txt", row.names=1, stringsAsFactors=TRUE)
spls<-levels(d$V2)

for(specie in spls){
dataPath<-paste0("../data/wlen/", specie, ".csv")
metadataPath<-paste0("../data/wlen/", specie, "_m.txt")

dataNL<-read.delim(dataPath, row.names=1, stringsAsFactors=TRUE)
metadata<-read.delim(metadataPath, header=T, row.names=1, stringsAsFactors=TRUE)
#dataNL<-read.delim("../data/data_nolen.csv", row.names=1, stringsAsFactors=TRUE)
#metadata<-read.delim("../data/metadata.txt", header=T, row.names=1, stringsAsFactors=TRUE)

colnames(metadata)<-c("specie", "quality", "tissue_abv", "rep", "location")


## ----getLen----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ("Length" %in% colnames(dataNL)){
  length_vec<-dataNL$Length
}


## ----dataSearch------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim(dataNL)
dim(metadata)


## ----dataMatch-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ----metadtaCheck----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Should be all 0 due to preprocessing filtering
levels(as.factor(metadata$quality))

#Mapped abreviations
levels(metadata$tissue_abv)

#reps 1,2,3 and 4, is there an imbalance?
table(metadata$rep)

levels(metadata$location)
#Different rep aoumts indicate different amount of each replicate
table(metadata$location) 
#there are different total numbers of tissue replicates
#Solve it by using 1 replicate per tissue (mean of exisitng replicates)



## ----OutlierGenes----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#@@ may not be a good idea to remove genes
outDetect<-goodSamplesGenes(t(dataNL))

table(outDetect$goodGenes) #False genes are outliers
table(outDetect$goodSamples) #All samples are True = no outliers


## ----OutlierGeRemove-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataNL<-dataNL[outDetect$goodGenes==TRUE,] #remove ouliers

if (exists("length_vec")){ #only if it exists
  length_vec<-length_vec[outDetect$goodGenes==TRUE] #if length_vec exists remove outliers from there as well
}



## ----replJoin--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rsumer<-function(data, metadata, tissue_name){ #calculates the mean of all columns that belong to a location (mean of all COB columns)
  
  loc_mdata<-metadata[metadata$location == tissue_name, ] #filter metadata tissue (get metadata of location only)
  
  data<-data[,colnames(data) %in% rownames(loc_mdata)] #get data of lcoation only, based on metadata
  
  data<-rowMeans(data) #calculate mean for each gene out of the locations(replicates)
  
  return(as.data.frame(data))
}

tissue_data<-levels(metadata$location) #get list of tissue names

d_joint<-sapply(tissue_data, function(tissue_name) rsumer(dataNL, metadata, tissue_name)) #returns an array where each entry is a column with the mean data of the replicates (rows are genes)

repl_data<-as.data.frame(d_joint) #data joint by replicate

colnames(repl_data) = gsub(pattern = "*.data", replacement = "", x = tolower(colnames(repl_data))) #get column names to be only location

rownames(repl_data)<-rownames(dataNL) #rename rows to be genes again



## ----replMeta--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
repl_meta<-as.data.frame(colnames(repl_data))
colnames(repl_meta)<-c("location")


## ----rawPlotPrep-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
multiColHist<-function(data, location, color_var){ #gets the data table, the concrete column(tissue) and the color for the tissue
  loc_var<-data[[location]] #get the tissue data
  p<-ggplot(data, aes(x=loc_var))+
    geom_histogram(bins = 30, fill=color_var, color="black")+
    xlab("")+ylab("")+theme_minimal() #plot a colored histogram for all genes of X tissue
  return(p)
}

colorList<-distinctColorPalette(nrow(repl_meta)) #make a list with a color for each tissue

#Create a legend that realates each color to a tissue

#use the repl_meta dataframe for creating the plot, as it contains all tissue(location) names
legend_plot<-ggplot(repl_meta, aes(x=1, y=location, color=location))+
  geom_point()+
  scale_color_manual(values=colorList)+
  guides(color=guide_legend(ncol=1))+
  theme_void()+labs(color="Tissue")+
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=12)) #plot used only to get the legend that associates colors with localizations(tissues)

legend_var<-get_legend(legend_plot) #place the legend into a variable


## ----rawPlotListing--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
RawHistList<-list() #to store the plots

for (i in 1:nrow(repl_meta)){
  color_var<-colorList[i]
  loc<-colnames(repl_data)[i]
  
  temp_plot<-multiColHist(repl_data, loc, color_var)
  RawHistList[[i]]<-temp_plot
} #iterates over each tissue, creates a plot with a distinct color for it and stores it in a list


## ----rawPlot---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw_Ptab<-arrangeGrob(grobs=RawHistList, ncol=nrow(repl_meta)/6) #creates a table that organizes the plots

png(paste0("./DistrPlots/", specie, "_raw_distPlot.png"), width=1600, height=800) #B73_CPM_raw_distPlot
grid.arrange(raw_Ptab, legend_var, widths = c(10, 2.3), ncol=2, top="Raw data distribution") #plots the plot list and legend together
dev.off()


## ----CPMnorm---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#uses edgeR
dge <- DGEList(repl_data)

#Calculate normalization factors
dge <- calcNormFactors(dge)

#Get normalized counts
Nrepl_data <- cpm(dge, log=TRUE) #rpkm with lengths, testo other normalizations too

#Nrepl to data frame
Nrepl_data<-as.data.frame(Nrepl_data) #evetually transpose

NormType<-"CPM"



## ----normPLotFunc----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#@@ analyze CPM and RPKN and choose 1
multiColLine<-function(data, location, color_var){ #gets the data table, the concrete column(tissue) and the color for the tissue
  loc_var<-data[[location]] #get the tissue data
  p<-ggplot(data, aes(x=loc_var))+geom_density(fill=color_var, color="black")+
    xlab("")+ylab("")+theme_minimal() #plot a colored line graph for all genes of X tissue
  return(p)
}


## ----CPMnormPlotListing----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

NormHistList<-list() #to store the plots

for (i in 1:nrow(repl_meta)){
  color_var<-colorList[i]
  loc<-colnames(repl_data)[i]
  
  temp_plot<-multiColHist(Nrepl_data, loc, color_var)
  NormHistList[[i]]<-temp_plot
} #iterates over each tissue, creates a plot with a distinct color for it and stores it in a list


## ----CPMnormPlot-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plotTitle<-paste0("Normalized data distribution using ", NormType)

norm_Ptab<-arrangeGrob(grobs=NormHistList, ncol=nrow(repl_meta)/6) #creates a table that organizes the plots

png(paste0("./DistrPlots/", specie, "_norm", NormType, "_distPlot.png"), width=1600, height=800) #B73_normCPM_distPlot
grid.arrange(norm_Ptab, legend_var, widths = c(10, 2.3), ncol=2, top=plotTitle) #plots the plot list and legend together
dev.off()


## ----RPKMnorm--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (exists("length_vec")){ #if we have lengths
 length_vec<-data.frame(Length=length_vec) #convert to dataframe

  dge <- DGEList(repl_data,genes=length_vec) #use edgeR for normalization

  dge <- calcNormFactors(dge)

  Nrepl_data <- rpkm(dge, log=TRUE)

  Nrepl_data<-as.data.frame(Nrepl_data)
  
  NormType<-"RPKM"
}


## ----RPKMnormPlotListing---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Plotting the normalized data counts <-allPLot Title
NormHistList<-list() #to store the plots

for (i in 1:nrow(repl_meta)){
  color_var<-colorList[i]
  loc<-colnames(repl_data)[i]
  
  temp_plot<-multiColHist(Nrepl_data, loc, color_var)
  NormHistList[[i]]<-temp_plot
} #iterates over each tissue, creates a plot with a distinct color for it and stores it in a list


## ----RPKMnormPlot----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plotTitle<-paste0("Normalized data distribution using ", NormType)

norm_Ptab<-arrangeGrob(grobs=NormHistList, ncol=nrow(repl_meta)/6) #creates a table that organizes the plots

png(paste0("./DistrPlots/", specie, "_norm", NormType, "_distPlot.png"), width=1600, height=800) #B73_normCPM_distPlot
grid.arrange(norm_Ptab, legend_var, widths = c(10, 2.3), ncol=2, top=plotTitle) #plots the plot list and legend together
dev.off()


## ----nwPowerList-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Nrepl_data<-t(Nrepl_data)

power<-c(c(1:15), seq(17, 50, by=2))

#Network topology analysis
sft <- pickSoftThreshold(Nrepl_data,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)





## ----nwPowerChoose---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Use automatic calculation
softPw <- sft$powerEstimate


## ----nwConstr--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
temp_cor <- cor
cor <- WGCNA::cor

ModNetwork <- blockwiseModules(Nrepl_data,
                 nThreads = 16, #32
                 maxBlockSize = 64000, #Memory dedicated to process (blocksize is 14000 with 16GB ram) (if Ngenes>maxBlocksize then Ngenes will be split into blocks to fit mBs)
                 deepSplit = 2,
                 TOMType = "signed", #unsigned?
                 power = softPw,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 42,
                 verbose = 4)

cor<-temp_cor


## ----EigenColors-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
module_eigengenes <- ModNetwork$MEs

#get number of genes for each module
table(ModNetwork$colors)

#Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(ModNetwork$dendrograms[[1]], cbind(ModNetwork$unmergedColors, ModNetwork$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


## ----binMeta---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
levels(as.factor(repl_meta$location))

bin_metadata <- data.frame(
  loc_cob=as.integer(repl_meta$location=="cob"),
  loc_coleoptile=as.integer(repl_meta$location=="coleoptile"),
  loc_crown_root=as.integer(repl_meta$location=="crown.root"),
  loc_first_elongated_internode=as.integer(repl_meta$location=="first.elongated.internode"),
  loc_flag_leaf=as.integer(repl_meta$location=="flag.leaf"),
  loc_immature_cob=as.integer(repl_meta$location=="immature.cob"),
  loc_immature_tassel=as.integer(repl_meta$location=="immature.tassel"),
  loc_leaf_1=as.integer(repl_meta$location=="leaf.1"),
  loc_leaf_3_blade=as.integer(repl_meta$location=="leaf.3.blade"),
  loc_leaf_3_sheath=as.integer(repl_meta$location=="leaf.3.sheath"),
  loc_leaf_5=as.integer(repl_meta$location=="leaf.5"),
  loc_leaf_5_elongation_zone=as.integer(repl_meta$location=="leaf.5.elongation.zone"),
  loc_leaf_5_mature=as.integer(repl_meta$location=="leaf.5.mature"),
  loc_leaf_5_meristem=as.integer(repl_meta$location=="leaf.5.meristem"),
  loc_leaf_8=as.integer(repl_meta$location=="leaf.8"),
  loc_mature_seed_40_dap=as.integer(repl_meta$location=="mature.seed.40.dap"),
  loc_meotic_tassel=as.integer(repl_meta$location=="meotic.tassel"),
  loc_mesophyll=as.integer(repl_meta$location=="mesophyll"),
  loc_prepollinated_cob=as.integer(repl_meta$location=="prepollinated.cob"),
  loc_primary_root=as.integer(repl_meta$location=="primary.root"),
  loc_primary_root_elongation_zone=as.integer(repl_meta$location=="primary.root.elongation.zone"),
  loc_primary_root_meristematic_zone=as.integer(repl_meta$location=="primary.root.meristematic.zone"),
  loc_root_hair_zone=as.integer(repl_meta$location=="root.hair.zone"),
  loc_seed_10_dap=as.integer(repl_meta$location=="seed.10.dap"),
  loc_seed_15_dap=as.integer(repl_meta$location=="seed.15.dap"),
  loc_seed_20_dap=as.integer(repl_meta$location=="seed.20.dap"),
  loc_seed_25_dap=as.integer(repl_meta$location=="seed.25.dap"),
  loc_seed_30_dap=as.integer(repl_meta$location=="seed.30.dap"),
  loc_seminal_root=as.integer(repl_meta$location=="seminal.root"),
  loc_silk=as.integer(repl_meta$location=="silk")
)

rownames(bin_metadata)<-rownames(Nrepl_data)



## ----corCalc---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nTissues <- nrow(Nrepl_data)
nGenes <- ncol(Nrepl_data)

MT_cor<-cor(module_eigengenes, bin_metadata, use="p")
module.trait.corr.pvals <- corPvalueStudent(MT_cor, nTissues)



## ----visCoexp--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
heatmap.data <- merge(module_eigengenes, bin_metadata, by = 'row.names')

head(heatmap.data)

rownames(heatmap.data)<-heatmap.data$Row.names
heatmap.data$Row.names<-NULL

head(heatmap.data)

plotTitle<-paste0("Correlation plot of ", specie, " using ", NormType)
png(paste0("./CorPlots/", specie, "_", NormType, "_corplot.png"), width=1600, height=800) #B73_CPM_corplot
CorLevelPlot(heatmap.data,
             x = names(bin_metadata),
             y = names(module_eigengenes),
             titleX="Modules", titleY="Traits", main=plotTitle,
             rotLabX = 45, rotTitleY = 90,
             cexCorval = 0.7, cexLabY = 0.7, cexLabX = 0.7
            )
dev.off()


}