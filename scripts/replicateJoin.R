dataNL<-read.delim("../data/B73.csv", row.names=1, stringsAsFactors=TRUE)
metadata<-read.delim("../data/B73_m.txt", header=T, row.names=1, stringsAsFactors=TRUE)

colnames(metadata)<-c("specie", "quality", "tissue_abv", "rep", "location")

rsumer<-function(data, metadata, tissue_name){ #calculates the mean of all columns thet belong to a location (mean of all COB columns)
  
  loc_mdata<-metadata[metadata$location == tissue_name, ] #filter metadata tissue (get metadata of location only)
  
  data<-data[,colnames(data) %in% rownames(loc_mdata)] #get data of lcoation only, based on metadata
  
  data<-rowMeans(data) #calculate mean for each gene out of the locations(replicates)
  
  return(as.data.frame(data))
}

tissue_data<-levels(metadata$location) #get list of tissue names

d_joint<-sapply(tissue_data, function(tissue_name) rsumer(dataNL, metadata, tissue_name)) #returns an array where each entry is a column with the mean data of the replicates (rows are genes)

repl_data<-as.data.frame(d_joint)

colnames(repl_data) = gsub(pattern = "*.data", replacement = "", x = colnames(repl_data))

rownames(repl_data)<-rownames(dataNL)

all(rownames(repl_data) == rownames(dataNL))

View(repl_data)
