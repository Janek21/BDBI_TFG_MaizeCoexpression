library(edgeR)
library(WGCNA)

dataPath<-paste0("../data/wlen/data_wlen.csv")
moduledata<-"./geneModule.txt"

dataNL<-read.delim(dataPath, row.names=1, stringsAsFactors=TRUE)
genemodule<-read.table(moduledata)

#Data selection
ME<-"blue"
genemodule$genes<-rownames(genemodule)
MEgeneList<-rownames(genemodule[genemodule$modules==ME,])
geneInfo<-dataNL[MEgeneList,]

#len vector
if ("Length" %in% colnames(geneInfo)){
  length_vec<-geneInfo$Length
}

##Join samples by replicate?

#Normalization
length_vec<-data.frame(Length=length_vec) #convert to dataframe

dge <- DGEList(geneInfo,genes=length_vec) #use edgeR for normalization

dge <- calcNormFactors(dge)
Nrepl_data <- rpkm(dge, log=TRUE)

Nrepl_data<-as.data.frame(Nrepl_data)

#filter data
keep<-apply(Nrepl_data, 1, max)>=0 #keep genes where the counts for at least one replicate are of at least 1 (0 because of log)
Nrepl_data<-Nrepl_data[keep,]

#Nrepl to data frame
Nrepl_data<-as.data.frame(Nrepl_data) #evetually transpose

#nw parameters
Nrepl_data<-t(Nrepl_data)

power <- c(c(1:50)) #more detailed in lower values

#Network topology analysis
sft <- pickSoftThreshold(Nrepl_data,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
softPw <- min(sft$powerEstimate, 30)
softPw<-30

#block build
temp_cor<-cor
cor<-WGCNA::cor

ModNetwork<-blockwiseModules(Nrepl_data,
                             nThreads = 32, #16
                             maxBlockSize = 64000, #directly related to memory, if maxBlockSize<total genes, multiple blocks will have to be used -> worse clustering
                             deepSplit = 4,
                             TOMType = "unsigned", #unsigned?
                             power = softPw,
                             mergeCutHeight = 0.3,#0.3->4799#0.8->4799(low module granularity) #0.1->4799(high granularity(lots of colors))
                             minModuleSize = 20,
                             numericLabels = FALSE,
                             pamRespectsDendro = FALSE,
                             minKMEtoStay = 0.3,#0.3->4799 #0.8->7781 #<0.3 stays the same
                             randomSeed = 42,
                             verbose = 4)


cor<-temp_cor


###################################################################

# Load necessary libraries
library(tidyverse)
library(igraph)

# Step 1: Load the count table
counts <- geneInfo
counts$Length<-NULL

# Step 2: Log-normalize the data
log_counts <- log1p(counts)

# Step 3: Compute correlation matrix (genes x genes)
cor_matrix <- cor(t(log_counts), method = "pearson")

# Step 4: Threshold the correlation matrix to form edges
threshold <- 0.4
adjacency_matrix <- abs(cor_matrix) > threshold
diag(adjacency_matrix) <- 0  # remove self-loops

# Step 5: Build igraph network
gene_network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)

# Step 6: Add absolute correlation as edge weights
edge_list <- as.data.frame(as_edgelist(gene_network))
weights <- mapply(function(g1, g2) abs(cor_matrix[g1, g2]), edge_list$V1, edge_list$V2)
E(gene_network)$weight <- weights

# Step 7: Plot the network (no error now)
plot(
  gene_network,
  vertex.label = NA,
  vertex.size = 3,
  edge.color = "gray",
  edge.width = E(gene_network)$weight * 2,  # optional: scale by weight
  main = "Gene Co-Expression Network"
)

