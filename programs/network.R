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
