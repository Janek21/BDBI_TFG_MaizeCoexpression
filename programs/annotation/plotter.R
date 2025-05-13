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
ME<-"skyblue3"
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

##Find smallest modules
