library(ggplot2)
library(stringr)
library(plotly)

#read data
pvdata<-read.table("./Pres/d_1_sheet.txt", sep='\t', header=TRUE)
#Choose relevant
pvdata$Relevancy<-as.logical(pvdata$Relevancy)
pvdata<-pvdata[pvdata$Relevancy,]

pvdata$significance<- -log10(pvdata$Pvalue)

lgn<-c(min(pvdata$Pvalue)+0.0000000001, max(pvdata$Pvalue)/2, max(pvdata$Pvalue))
                      
rr<-ggplot(pvdata, aes(x=Function, y=Module, size=significance, color=significance))+
  geom_point()+
  scale_size_continuous(name="P-value",
                        breaks=-log10(lgn),
                        labels=lgn)+
  scale_color_gradient(name="P-value",
                         breaks=-log10(lgn),
                         labels=lgn,
                       low = "blue", high = "red")+
  labs(x="Genes", y="Functions")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust =1, size=10, color = "black"),
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
        axis.ticks=element_line(colour="black"),
        axis.ticks.length = unit(8, "pt"))

ggplotly(rr)

##??
ggplot(pvdata, aes(x=Function, y=Module, size=Pvalue, color=Pvalue))+
  geom_point()+
  scale_size(trans="reverse")+
  scale_color_continuous(trans="reverse")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust =1, size=10, color = "black"),
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
        axis.ticks=element_line(colour="black"),
        axis.ticks.length = unit(8, "pt"))


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