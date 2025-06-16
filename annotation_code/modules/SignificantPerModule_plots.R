
#Plots significance per function in each module, highlights most significant in total and writes them in a table
#Step 3

library(ggplot2)
library(plotly)
library(ggnewscale)

#read data
dChoice<-"anT_mercator" #anT_mercator, anT_prot-scriber, anT_swissprot
dChoice<-"d_1"
pvdata<-read.table(paste0("../modules/Pres/", dChoice, "_sheet.txt"), sep='\t', header=TRUE)
#Choose relevant
pvdata$Relevancy<-as.logical(pvdata$Relevancy)
pvdata<-pvdata[pvdata$Relevancy,]

#Get top x significant values
highSig<-sort(pvdata$Pvalue)[1:10]
pvdata$TopSig<-pvdata$Pvalue %in% highSig

#plot all values
bareplot<-ggplot(pvdata, aes(x=Function, y=Module))+
  geom_point(aes(size=Pvalue, color=Pvalue))+
  scale_size(range = c(0, 5), trans="reverse")+
  scale_color_continuous(trans="reverse")+
  labs(size="P-value", color="P-value")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        legend.key.width=unit(1, "cm"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=16),
        plot.background=element_rect(fill="transparent",color="transparent"),
        panel.border=element_rect(linewidth=1),
        panel.grid.major=element_line(color="grey"),
        axis.ticks=element_line(colour="black"),
        axis.ticks.length = unit(8, "pt"))+
  guides(color=guide_legend(order=1),
         size=guide_legend(order=1))


#highlight significant values(all equal)
hplot<-bareplot+
  geom_point(data=subset(pvdata, TopSig), aes(fill=TopSig), color="darkorange", size=7)+
  labs(fill="")+scale_fill_hue(labels=c("Highly\nsignificant"))+
  guides(fill=guide_legend(order=2))+
  theme(legend.text=element_text(size=14))

png(paste0("../modules/SigPlots/sig_", dChoice, ".png"), width=1600, height=800)
hplot#+ggtitle(dChoice)
dev.off()
svg(paste0("../modules/SigPlots/sig_", dChoice, ".svg"), width=1600/60, height=800/60)
hplot#+ggtitle(dChoice)
dev.off()

#ggplotly(hplot)

#Set up table for significant functions in each module, remove unneeded data
TopSig<-as.data.frame(subset(pvdata, TopSig))
TopSig$TotalCounts<-TopSig$ModuleCounts<-TopSig$Relevancy<-TopSig$TopSig<-NULL
TopSig$Pvalue<-as.numeric(formatC(TopSig$Pvalue, digits=2))

TopPlot<-ggplot(TopSig, aes(x=Function, y=Module, color=Pvalue, size=Pvalue))+
  geom_point()+
  scale_color_continuous(trans="reverse")+#, breaks=c(min(TopSig$Pvalue)))+
  scale_size(range = c(0, 7), trans="reverse")+
  labs(size="P-value", color="P-value")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        legend.key.width=unit(1, "cm"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=16),
        plot.background=element_rect(fill="transparent",color="transparent"),
        panel.border=element_rect(linewidth=1),
        panel.grid.major=element_line(color="grey"),
        axis.ticks=element_line(colour="black"),
        axis.ticks.length = unit(8, "pt"))+
  guides(color=guide_legend(order=1),
         size=guide_legend(order=1))

png(paste0("../modules/SigPlots/top/top_", dChoice, ".png"), width=1600, height=800)
TopPlot#+ggtitle(dChoice)
dev.off()
svg(paste0("../modules/SigPlots/top/top_", dChoice, ".svg"), width=1600/60, height=800/60)
TopPlot#+ggtitle(dChoice)
dev.off()

#write top significant to table
write.table(TopSig, paste0("../modules/SigPlots/top/top_", dChoice, ".txt"), sep='\t', row.names=FALSE)

