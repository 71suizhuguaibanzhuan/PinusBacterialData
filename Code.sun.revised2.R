##Figure.2a.Barplot.phylum.otu.genus.treatment####
library(vegan)
library(betapart)
library(colorRamps)
library(agricolae)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
#library(learnasreml)
##Bacteria-AIR,ROOT,LEAF
setwd("~/Documents/sunR/Bacteria882")

rm(list = ls())

##1.data
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##2.phylum
bact.phy <- aggregate(otu_all1,by=list(ID$phylum) , sum) 
rownames(bact.phy) <- bact.phy[,1]; #fungi.lev rownames
bact.phy <- bact.phy[,-1] 
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

bact.phy0 <- bact.phy0[,order(-colSums(bact.phy0))] ## 
lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.phy.lev <- aggregate(bact.phy0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

library(splitstackshape)
bact.phy1 <- bact.phy.lev[,c(1,2:11)] # the domiant OTUs 前十的门
bact.phy1 <- melt(bact.phy1,id.vars = "Group.1")
bact.phy1 <- cSplit(bact.phy1, "Group.1", ":") 

names(bact.phy1)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")

bact.phy2 <- bact.phy.lev[,c(1,12:ncol(bact.phy.lev))] # the rare OTUs, combined as 'others'
bact.phy2 <- melt(bact.phy2,id.vars = "Group.1")
bact.phy2 <- cSplit(bact.phy2, "Group.1", ":")
names(bact.phy2) <- c("Bact","Relative_Abundance", "Time", "Type","Treatment")
bact.phy2$Bact <- "Other"

bact.bind1 <- rbind(bact.phy1,bact.phy2) # combine the domiant and rare OTUs

col17<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378",
          "lightcoral", "rosybrown", "#43A743","#97D1A0",  "grey",
          "#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", 
          "darkgreen","maroon3", "black", "bisque", "grey")

p1 <- ggplot(bact.bind1, aes(x = Time, y = Relative_Abundance, fill=Bact)) +
  geom_bar(stat='identity', position = "fill")+
  #geom_flow(aes(alluvium = Bacteria),alpha=0.5)+
  facet_grid(Treatment~Type)+
  scale_fill_manual(values= col17)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "Phylum"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 15,face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.y=element_text(colour="black",size=10,face="bold"),
        axis.text.x=element_text(colour="black",size=12,face="bold"),
        axis.title=element_text(colour="black",size=16,face="bold"))


p1
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-treatment.barplot.phylum.pdf", p1, height = 7, width = 12)


##7.otu
ID$otu <- paste(ID$names, ID$genus, sep = ".")

bact.otu <- aggregate(otu_all1,by=list(ID$otu) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )

bact.otu0 <- bact.relabu.otu

bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))] ## 
lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.otu.lev <- aggregate(bact.otu0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

bact.otu1 <- bact.otu.lev[,c(1,2:11)] # 
bact.otu1 <- melt(bact.otu1,id.vars = "Group.1")
bact.otu1 <- cSplit(bact.otu1, "Group.1", ":") 

names(bact.otu1)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")

bact.otu2 <- bact.otu.lev[,c(1,12:ncol(bact.otu.lev))] # the rare OTUs, combined as 'others'
bact.otu2 <- melt(bact.otu2,id.vars = "Group.1")
bact.otu2<-cSplit(bact.otu2, "Group.1", ":")
names(bact.otu2)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")
bact.otu2$Bact <- "Other"

bact.bind6 <- rbind(bact.otu1,bact.otu2) # combine the domiant and rare OTUs


p2 <- ggplot(bact.bind6, aes(x = Time, y = Relative_Abundance, fill=Bact)) +
  geom_bar(stat='identity', position = "fill")+
  #geom_flow(aes(alluvium = Bacteria),alpha=0.5)+
  facet_grid(Treatment~Type)+
  scale_fill_manual(values= col17)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "OTU"))+
  scale_y_continuous(labels = scales::percent)+
  scale_y_cut(breaks = 0.65, scales = 20, which = 1, expand = T) +
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 15,face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.y=element_text(colour="black",size=10,face="bold"),
        axis.text.x=element_text(colour="black",size=12,face="bold"),
        axis.title=element_text(colour="black",size=16,face="bold"))


p2
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-treatment.barplot.otu.pdf", p2, height = 7, width = 12)

##2.genus
bact.phy <- aggregate(otu_all1,by=list(ID$genus) , sum) 
rownames(bact.phy) <- bact.phy[,1]; 
bact.phy <- bact.phy[,-1] 
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

bact.phy0 <- bact.phy0[,order(-colSums(bact.phy0))] ## 
lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.phy.lev <- aggregate(bact.phy0, by = list(lev), FUN = "mean") # generate the mean for each factor level
library(dplyr)
bact.phy.lev1 <- bact.phy.lev %>% select(-uncultured, everything())

library(splitstackshape)
bact.phy1 <- bact.phy.lev1[,c(1,2:15)] 
bact.phy1 <- melt(bact.phy1,id.vars = "Group.1")
bact.phy1 <- cSplit(bact.phy1, "Group.1", ":") 

names(bact.phy1)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")

bact.phy2 <- bact.phy.lev1[,c(1,16:ncol(bact.phy.lev1))] # the rare OTUs, combined as 'others'
bact.phy2 <- melt(bact.phy2,id.vars = "Group.1")
bact.phy2 <- cSplit(bact.phy2, "Group.1", ":")
names(bact.phy2) <- c("Bact","Relative_Abundance", "Time", "Type","Treatment")
bact.phy2$Bact <- "Other"

bact.bind1 <- rbind(bact.phy1,bact.phy2) # combine the domiant and rare OTUs

col17<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378",
          "lightcoral", "rosybrown", "#43A743","#97D1A0",
          "maroon3", "navy", "deepskyblue","darkgreen", "grey",
          "#ff00ff","#00ff00","gold",  "red", 
          "black", "bisque", "grey")

p7 <- ggplot(bact.bind1, aes(x = Time, y = Relative_Abundance, fill=Bact)) +
  geom_bar(stat='identity', position = "fill")+
  #geom_flow(aes(alluvium = Bacteria),alpha=0.5)+
  facet_grid(Treatment~Type)+
  scale_fill_manual(values= col17)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "Genus"))+
  scale_y_continuous(labels = scales::percent)+
  scale_y_cut(breaks = 0.45, scales = 20, which = 1, expand = T) +
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 15,face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.y=element_text(colour="black",size=10,face="bold"),
        axis.text.x=element_text(colour="black",size=12,face="bold"),
        axis.title=element_text(colour="black",size=16,face="bold"))


p7
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-treatment.barplot.genus.pdf", p7, height = 7, width = 12)

##2.family
bact.phy <- aggregate(otu_all1,by=list(ID$family) , sum) 
rownames(bact.phy) <- bact.phy[,1]; 
bact.phy <- bact.phy[,-1]
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

bact.phy0 <- bact.phy0[,order(-colSums(bact.phy0))] ## 
lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.phy.lev <- aggregate(bact.phy0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

library(splitstackshape)
bact.phy1 <- bact.phy.lev[,c(1,2:15)] 
bact.phy1 <- melt(bact.phy1,id.vars = "Group.1")
bact.phy1 <- cSplit(bact.phy1, "Group.1", ":") 

names(bact.phy1)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")

bact.phy2 <- bact.phy.lev[,c(1,16:ncol(bact.phy.lev))] # the rare OTUs, combined as 'others'
bact.phy2 <- melt(bact.phy2,id.vars = "Group.1")
bact.phy2 <- cSplit(bact.phy2, "Group.1", ":")
names(bact.phy2) <- c("Bact","Relative_Abundance", "Time", "Type","Treatment")
bact.phy2$Bact <- "Other"

bact.bind1 <- rbind(bact.phy1,bact.phy2) # combine the domiant and rare OTUs

col17<-c( "#807EBA","#A7B7DF", "#ABDAEC","#E9CEE5","#FDD5C0","#FDD378",
          "lightcoral", "rosybrown", "#43A743","#97D1A0",  
          "#ff00ff","#00ff00", "deepskyblue", "gold", "grey","red", "navy", 
          "darkgreen","maroon3", "black", "bisque", "grey")

p9 <- ggplot(bact.bind1, aes(x = Time, y = Relative_Abundance, fill=Bact)) +
  geom_bar(stat='identity', position = "fill")+
  #geom_flow(aes(alluvium = Bacteria),alpha=0.5)+
  facet_grid(Treatment~Type)+
  scale_fill_manual(values= col17)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "Family"))+
  scale_y_continuous(labels = scales::percent)+
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 15,face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.y=element_text(colour="black",size=10,face="bold"),
        axis.text.x=element_text(colour="black",size=12,face="bold"),
        axis.title=element_text(colour="black",size=16,face="bold"))


p9
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-treatment.barplot.family.pdf", p1, height = 7, width = 12)


##7.order
bact.otu <- aggregate(otu_all1,by=list(ID$order) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )

bact.otu0 <- bact.relabu.otu

bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))] ## 
lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.otu.lev <- aggregate(bact.otu0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

bact.otu1 <- bact.otu.lev[,c(1,2:15)] 
bact.otu1 <- melt(bact.otu1,id.vars = "Group.1")
bact.otu1 <- cSplit(bact.otu1, "Group.1", ":") 

names(bact.otu1)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")

bact.otu2 <- bact.otu.lev[,c(1,16:ncol(bact.otu.lev))] # the rare OTUs, combined as 'others'
bact.otu2 <- melt(bact.otu2,id.vars = "Group.1")
bact.otu2<-cSplit(bact.otu2, "Group.1", ":")
names(bact.otu2)<-c("Bact","Relative_Abundance", "Time", "Type","Treatment")
bact.otu2$Bact <- "Other"

bact.bind6 <- rbind(bact.otu1,bact.otu2) # combine the domiant and rare OTUs


p2 <- ggplot(bact.bind6, aes(x = Time, y = Relative_Abundance, fill=Bact)) +
  geom_bar(stat='identity', position = "fill")+
  #geom_flow(aes(alluvium = Bacteria),alpha=0.5)+
  facet_grid(Treatment~Type)+
  scale_fill_manual(values= col17)+
  theme_bw()+
  labs(x="",y = "Relative abundance")+
  guides(fill=guide_legend(title= "Order"))+
  scale_y_continuous(labels = scales::percent)+
  #scale_y_cut(breaks = 0.2, scales = 20, which = 1, expand = T) +
  theme(panel.spacing = unit(0, "lines"),
        strip.text = element_text(size = 15,face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text.y=element_text(colour="black",size=10,face="bold"),
        axis.text.x=element_text(colour="black",size=12,face="bold"),
        axis.title=element_text(colour="black",size=16,face="bold"))


p2
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-treatment.barplot.order.pdf", p2, height = 7, width = 12)




##Figure.2b.PCoA####
library(vegan)
library(rmarkdown)
library(colorRamps)
library(splitstackshape)
library(ggnewscale)
library(ggplot2)
library(ape)
library(umap)
library(plyr)
library(plotly)
library(tidyverse)
library(readxl)
library(writexl)

##Bacteria-PCoA-7
setwd("~/Documents/sunR/Bacteria882/")
rm(list = ls())

##sample_metadata
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),] 
group1 <- group[, c("sample.id","sampleType","TP0","TP00","Treatment")]
group2 <- na.omit(group1)
group2$t1 <- paste(group2$sampleType,group2$Treatment,sep="_")

##otu table
otu_all0 <- read.csv("otu_bactFlattening825.csv",check.names = F,row.names = 1,header = T)

otu_all0 <- otu_all0[, group2$sample.id]
otu_all0$sum <- rowSums(otu_all0)
zero <- which(otu_all0$sum == 0)
otu_all <- otu_all0[-zero, ]
otu_all <- otu_all[, -826]

otu_all <- as.data.frame(t(otu_all))
otu_all$sample_id <- row.names(otu_all)

otu1 <- otu_all <- merge(otu_all, group2, by.x = "sample_id", by.y = "sample.id")
#####################################################################PCoA
pcoa <- cmdscale(vegdist(decostand(otu1[,2:5212], "hellinger"), method = 'bray'), 
                 k = (nrow(otu1) - 1), eig = TRUE)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste0(round(100*pcoa_exp[1], 2) ,'%')
pcoa2 <- paste0(round(100*pcoa_exp[2], 2) ,'%')
pcoa3 <- paste0(round(100*pcoa_exp[3], 2) ,'%')

permanova <- adonis2(decostand(otu1[,2:5212],"hellinger") ~ otu1$TP00+otu1$sampleType+otu1$Treatment,
                     by="margin",distance = 'bray', permutations = 999,parallel=60)###perm anova### 
print(permanova)

permanova1 <- adonis2(decostand(otu1[,2:5212],"hellinger") ~ 
                        otu1$TP00*otu1$sampleType + otu1$TP00*otu1$Treatment +
                        otu1$sampleType*otu1$Treatment,
                      by="margin",distance = 'bray', permutations = 999,parallel=60)###perm anova### 
print(permanova1)

permanova2 <- adonis2(decostand(otu1[,2:5212],"hellinger") ~ 
                        otu1$TP00*otu1$sampleType*otu1$Treatment,
                      by="margin",distance = 'bray', permutations = 999,parallel=60)###perm anova### 
print(permanova2)


time1 <- data.frame(pcoa$point)[1:2]
time1$sample <- rownames(time1)
time1$sample <- otu1$sample_id
time1 <- merge(time1, group2, by.x = 'sample',by.y = "sample.id")
names(time1)[2:3] <- c('pcoa1', 'pcoa2')

h1 <- max(time1$pcoa1)-0.77*(range(time1$pcoa1)[2]-range(time1$pcoa1)[1])

#save.image("fungi-pcoa.RData")
##########################################################################PCoA
polygon <- time1[,c("pcoa1", "pcoa2", "t1")]
find_hull <- function(df1) df1[chull(df1$pcoa1, df1$pcoa2), ]
hulls <- ddply(polygon, "t1", find_hull)

p7 <- ggplot(data = time1, aes(pcoa1, pcoa2)) +
  geom_polygon(data = hulls, alpha = 0.15, aes(fill = t1)) +
  geom_point(aes(color = TP00, shape = t1)) +
  scale_shape_manual(name="Treatment",values = c(0,1,2,1,2),guide=guide_legend(order=2)) +
  scale_color_manual(values = blue2red(7)) +
  scale_fill_manual(values = c("blue","darkgreen","#5DCA3B","darkred","chocolate")) +
  theme_bw()+
  guides(color = guide_legend(title = "Time"), 
         shape = guide_legend(title = "Treatment"), 
         fill = guide_legend(title = "Treatment"))+
  theme(panel.grid = element_blank(),  #
        #panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(colour="black", size=13, face="bold", hjust = 0.5), ##legend.position = "none" 
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text = element_text(color = 'black',size=10, face="bold"), 
        axis.title = element_text(color = 'black',size=16, face="bold"),
        axis.ticks = element_line(color = 'black')) +
  geom_vline(xintercept = 0, color = 'black', size = 0.5, linetype = 2) +
  geom_hline(yintercept = 0, color = 'black', size = 0.5, linetype = 2) +
  labs(x=paste("PCo 1 (", pcoa1, ")", sep=""),
       y=paste("PCo 2 (", pcoa2, ")", sep="")) +
  annotate('text', label = paste("Time: R^2 = ", round(permanova$R2[1],3), "***", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.29*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Compartment: R^2 = ",   round(permanova$R2[2],3), "***", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.25*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Treatment: R^2 = ",   round(permanova$R2[3],3), "**", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.21*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Compartment: R^2 = ", round(permanova1$R2[1],3), "***", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.17*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Treatment: R^2 = ", round(permanova1$R2[2],3), "", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.13*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Compartment*Treatment: R^2 = ", round(permanova1$R2[3],3), "*", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.09*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Compartment*Treatment: R^2 = ", round(permanova2$R2[1],3), "", sep=""), 
           x = h1, y = min(time1$pcoa2)+0.05*(range(time1$pcoa2)[2]-range(time1$pcoa2)[1]), size = 4,color="black",fontface="bold") 



p7
#setwd("~/Documents/sunR/Bacteria882/n.treatment/")
#ggsave("B-pcoa.n.pdf",p7,height = 7, width = 9)
##########################################################################PCo3
time2 <- data.frame(pcoa$point)[1:3]
time2$sample <- rownames(time2)
time2$sample <- otu1$sample_id
time2 <- merge(time2, group2, by.x = 'sample',by.y = "sample.id")
names(time2)[2:4] <- c('pcoa1', "pcoa2",'pcoa3')

h2 <- max(time2$pcoa1)-0.75*(range(time2$pcoa1)[2]-range(time2$pcoa1)[1])
polygon <- time2[,c("pcoa1", "pcoa3", "t1")]
find_hull <- function(df1) df1[chull(df1$pcoa1, df1$pcoa3), ]
hulls <- ddply(polygon, "t1", find_hull)

p3 <- ggplot(data = time2, aes(pcoa1, pcoa3)) +
  geom_polygon(data = hulls, alpha = 0.15, aes(fill = t1)) +
  geom_point(aes(color = TP00, shape = t1)) +
  #stat_ellipse(aes(fill = sampleTime2), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +      
  scale_shape_manual(name="Treatment",values = c(0,1,2,1,2),guide=guide_legend(order=2)) +
  scale_color_manual(values = blue2red(7)) +
  scale_fill_manual(values = c("blue","darkgreen","#5DCA3B","darkred","chocolate")) +
  theme_bw()+
  guides(color = guide_legend(title = "Time"), 
         shape = guide_legend(title = "Treatment"), 
         fill = guide_legend(title = "Treatment"))+
  theme(panel.grid = element_blank(),  #
        #panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(colour="black", size=13, face="bold", hjust = 0.5), ##legend.position = "none" 
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        axis.text = element_text(color = 'black',size=10, face="bold"), 
        axis.title = element_text(color = 'black',size=16, face="bold"),
        axis.ticks = element_line(color = 'black')) +
  geom_vline(xintercept = 0, color = 'black', size = 0.5, linetype = 2) +
  geom_hline(yintercept = 0, color = 'black', size = 0.5, linetype = 2) +
  labs(x=paste("PCo 1 (", pcoa1, ")", sep=""),
       y=paste("PCo 3 (", pcoa3, ")", sep="")) +
  annotate('text', label = paste("Time: R^2 = ", round(permanova$R2[1],3), "***", sep=""), 
           x = h2, y = min(time2$pcoa2)+1.13*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Compartment: R^2 = ",   round(permanova$R2[2],3), "***", sep=""), 
           x = h2, y = min(time2$pcoa2)+1.09*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Treatment: R^2 = ",   round(permanova$R2[3],3), "**", sep=""), 
           x = h2, y = min(time2$pcoa2)+1.05*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Compartment: R^2 = ", round(permanova1$R2[1],3), "***", sep=""), 
           x = h2, y = min(time2$pcoa2)+1.01*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Treatment: R^2 = ", round(permanova1$R2[2],3), "", sep=""), 
           x = h2, y = min(time2$pcoa2)+0.97*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Compartment*Treatment: R^2 = ", round(permanova1$R2[3],3), "*", sep=""), 
           x = h2, y = min(time2$pcoa2)+0.93*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") +
  annotate('text', label = paste("Time*Compartment*Treatment: R^2 = ", round(permanova2$R2[1],3), "", sep=""), 
           x = h2, y = min(time2$pcoa2)+0.89*(range(time2$pcoa2)[2]-range(time2$pcoa2)[1]), size = 4,color="black",fontface="bold") 

p3

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-pcoa.n.轴3.pdf",p3,height = 7, width = 9)


##Figure.2c.Mantel test####
##2025.3.13
library(vegan)
library(ggplot2)
library(ecodist)
library(patchwork)
library(readxl)
library(writexl)
###Bacteria-Bray-curtis dissimilarity-Mantel test
rm(list=ls())

setwd("~/Documents/sunR/Bacteria882")
##
otu_all0 <- read.csv("otu_bactFlattening825.csv",check.names = F,row.names = 1,header = T)
otu_all0 <- as.data.frame(t(otu_all0))  
##
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
otu_all0$sample.id <- rownames(otu_all0)
otu_all <- merge(otu_all0,group,by= "sample.id")

group1 <- group[, c("sample.id","sampleType","TP0","Treatment")]
group2 <- na.omit(group1)
otu7 <- merge(otu_all0, group2, by = "sample.id")


####air.7
beta7 <- vegdist(decostand(otu7[which(otu7$sampleType == "Air"),2:5295], "hellinger"), method = 'bray')
td7 <- distance(otu7[which(otu7$sampleType == "Air"), 5297])

df7 <- as.data.frame(cbind(beta7,td7))
df7$beta7 <- as.numeric(df7$beta7)
df7$td7 <- as.factor(df7$td7)

ma7 <- mantel(beta7~td7, nperm = 999)

df.mantel7 <- data.frame(
  label = sprintf(" 'Mantel' ~ italic(R) ~ '=' ~ %.3g ~ '; ' ~ italic(p) ~ '=' ~ %.3g", ma7[1],ma7[2]), 
  x = 1,
  y=0.85 )


p7 <- ggplot(df7,aes(x=td7,y=beta7)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.3, alpha = 0.2) +
  labs(x = "", y = "Bray-Curtis dissimilarity") +
  ylim(0,1)+
  #ggtitle("Air") +
  #geom_text(data = df.mantel7,aes(label = label, x = 1, y = 0.1),parse = T,size = 5,color = "red", hjust = 0)+
  theme_bw() +
  theme(axis.title = element_text(size = 16,color="black",face = "bold"),
        axis.text = element_text(size = 11,color="black",face = "bold")) +
  annotate("text", x = 1, y = 0.1, label = "Air: R = 0.488, P = 0.001", size=5.5, hjust = 0, color="red")

p7


####root
root1 <- otu7 %>% filter(sampleType=="Root") %>% filter(Treatment=="Control")

beta2 <- vegdist(decostand(root1[,2:5295], "hellinger"), method = 'bray')
td2 <- distance(root1[, 5297])

df2 <- as.data.frame(cbind(beta2,td2))
df2$beta2 <- as.numeric(df2$beta2)
df2$td2 <- as.factor(df2$td2)

ma2 <- mantel(beta2~td2, nperm = 999)

df.mantel2 <- data.frame(
  label = sprintf(" 'Mantel' ~ italic(R) ~ '=' ~ %.3g ~ '; ' ~ italic(p) ~ '=' ~ %.3g", ma2[1],ma2[2]), 
  x = 1,
  y=0.85 )


p2 <- ggplot(df2,aes(x=td2,y=beta2)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.3, alpha = 0.2) +
  labs(x = "", y = "") +
  ylim(0,1)+
  #ggtitle("Root") +
  #geom_text(data = df.mantel2,aes(label = label, x = 1, y = 0.1),parse = T,size = 5,color = "red", hjust = 0)+
  theme_bw() +
  theme(axis.title = element_text(size = 16,color="black",face = "bold"),
        axis.text = element_text(size = 11,color="black",face = "bold")) +
  annotate("text", x = 1, y = 0.1, label = "Root_Control: \n          R = 0.724, P = 0.001",size=5.5, hjust = 0, color="red")

p2


root2 <- otu7 %>% filter(sampleType=="Root") %>% filter(Treatment=="N-addition")

beta2 <- vegdist(decostand(root2[,2:5295], "hellinger"), method = 'bray')
td2 <- distance(root2[, 5297])

df2 <- as.data.frame(cbind(beta2,td2))
df2$beta2 <- as.numeric(df2$beta2)
df2$td2 <- as.factor(df2$td2)

ma2 <- mantel(beta2~td2, nperm = 999)

df.mantel2 <- data.frame(
  label = sprintf(" 'Mantel' ~ italic(R) ~ '=' ~ %.3g ~ '; ' ~ italic(p) ~ '=' ~ %.3g", ma2[1],ma2[2]), 
  x = 1,
  y=0.85 )


p2n <- ggplot(df2,aes(x=td2,y=beta2)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.3, alpha = 0.2) +
  labs(x = "", y = "") +
  ylim(0,1)+
  #ggtitle("Root") +
  #geom_text(data = df.mantel2,aes(label = label, x = 1, y = 0.1),parse = T,size = 5,color = "red", hjust = 0)+
  theme_bw() +
  theme(axis.title = element_text(size = 16,color="black",face = "bold"),
        axis.text = element_text(size = 11,color="black",face = "bold")) +
  annotate("text", x = 1, y = 0.1, label = "Root_N-addition: \n          R = 0.554, P = 0.001",size=5.5, hjust = 0, color="red")

p2n



####leaf
leaf1 <- otu7 %>% filter(sampleType=="Leaf") %>% filter(Treatment=="Control")

beta3 <- vegdist(decostand(leaf1[,2:5295], "hellinger"), method = 'bray')
td3 <- distance(leaf1[, 5297])

df3 <- as.data.frame(cbind(beta3,td3))
df3$beta3 <- as.numeric(df3$beta3)
df3$td3 <- as.factor(df3$td3)

ma3 <- mantel(beta3~td3, nperm = 999)

df.mantel3 <- data.frame(
  label = sprintf(" 'Mantel' ~ italic(R) ~ '=' ~ %.3g ~ '; ' ~ italic(p) ~ '=' ~ %.3g", ma3[1],ma3[2]), 
  x = 1,
  y=0.85 )


p3 <- ggplot(df3,aes(x=td3,y=beta3)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.3, alpha = 0.2) +
  labs(x = "", y = "") +
  ylim(0,1)+
  #ggtitle("Leaf") +
  #geom_text(data = df.mantel3,aes(label = label, x = 1, y = 0.1),parse = T,size = 5,color = "red", hjust = 0)+
  theme_bw() +
  theme(axis.title = element_text(size = 16,color="black",face = "bold"),
        axis.text = element_text(size = 11,color="black",face = "bold")) +
  annotate("text", x = 1, y = 0.1, label = "Leaf_Control: \n         R = 0.455, P = 0.001",size=5.5, hjust = 0, color="red")

p3



leaf2 <- otu7 %>% filter(sampleType=="Leaf") %>% filter(Treatment=="N-addition")

beta3 <- vegdist(decostand(leaf2[,2:5295], "hellinger"), method = 'bray')
td3 <- distance(leaf2[, 5297])

df3 <- as.data.frame(cbind(beta3,td3))
df3$beta3 <- as.numeric(df3$beta3)
df3$td3 <- as.factor(df3$td3)

ma3 <- mantel(beta3~td3, nperm = 999)

df.mantel3 <- data.frame(
  label = sprintf(" 'Mantel' ~ italic(R) ~ '=' ~ %.3g ~ '; ' ~ italic(p) ~ '=' ~ %.3g", ma3[1],ma3[2]), 
  x = 1,
  y=0.85 )


p3n <- ggplot(df3,aes(x=td3,y=beta3)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 0.3, alpha = 0.2) +
  labs(x = "Temporal distance (month)", y = "") +
  ylim(0,1)+
  #ggtitle("Leaf") +
  #geom_text(data = df.mantel3,aes(label = label, x = 1, y = 0.1),parse = T,size = 5,color = "red", hjust = 0)+
  theme_bw() +
  theme(axis.title = element_text(size = 16,color="black",face = "bold"),
        axis.text = element_text(size = 11,color="black",face = "bold")) +
  annotate("text", x = 1, y = 0.1, label = "Leaf_N-addition: \n         R = 0.431, P = 0.001",size=5.5, hjust = 0, color="red")

p3n


library(gridExtra)

grid.arrange(
  p7, p3, p2, 
  p3n, p2n,    
  ncol = 3,    
  widths = c(1, 1, 1) 
)
empty_grob <- grid::nullGrob()
plots <- list(
  p7,p3,p2,  
  empty_grob,p3n,p2n  
)

p <- grid.arrange(grobs = plots, ncol = 3, widths = c(1, 1, 1))


#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-mantel.n.pdf", p, height = 8, width = 12)




##Figure.2d.Treatment.Root.otu18.otu36####
library(vegan)
library(betapart)
library(colorRamps)
library(agricolae)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(dplyr)

setwd("~/Documents/sunR/Bacteria882")

rm(list = ls())

##1.data
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]
group1 <- group[,c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
bact.otu <- aggregate(otu_all1,by=list(ID$names) , sum) 
rownames(bact.otu) <- bact.otu[,1]
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )

bact.otu0 <- bact.relabu.otu
bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))]

colnames(bact.otu0)
bact.otu10 <- bact.otu0[, c(1:10)]
bact.otu10$sample.id <- rownames(bact.otu10)

b.otu10 <- merge(bact.otu10, group2, by="sample.id")

b.otu <- gather(b.otu10,key = "OTU",value = "Abundance",-"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")


b.otu$value <- b.otu$Abundance*100

###OTU36###
otu36 <- b.otu[b.otu$OTU=="OTU36",]
otu36.r <- otu36[otu36$sampleTyp=="Root",]

root1 <- otu36.r %>% filter(Treatment=="Control")

lm1 <- lm(value~TP0, root1)
shapiro.test(residuals(lm1))
hist(residuals(lm1))

root2 <- otu36.r %>% filter(Treatment=="N-addition")

lm2 <- lm(value~TP0, root2)
shapiro.test(residuals(lm2))
hist(residuals(lm2))

root_otu1 <- otu36.r[otu36.r$TP00=="TP6",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu36.r$value, otu36.r$Treatment, p.adj = "bonferroni")
kw

root_otu2 <- otu36.r[otu36.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 23.87689, df = 1, p-value = 1.026972e-06

p.36r1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU36") + 
  ylim(-1, 26)+
  geom_text(x=1.5,y=22,label=expression(paste(chi^2, " = ", 23.878, ", ", "Df = ", 1, ", ", "P < ", 0.001, sep="")),
            size=6.9,color="black",fontface="bold")+
  annotate("text",x=1,y=7.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=7.2,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=19,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p.36r1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p.36r2 <- ggplot(otu36.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(-1, 26)+
  annotate("text",x=2,y=22,label="*",size=12,color="black")+
  annotate("text",x=4,y=6,label="*",size=12,color="black")+
  annotate("text",x=5,y=5,label="*",size=12,color="black")+
  annotate("text",x=6,y=5,label="*",size=12,color="black")+
  annotate("text", x=2.3, y=21, label="OTU36", 
           size=7, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.3, y=18, label="Genome Size: 5.778 Mb", 
           size=7, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.3, y=15, label="RRN: 3", 
           size=7, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p.36r2

##OTU36.lm
otu36 <- b.otu[b.otu$OTU=="OTU36",]
otu36.r <- otu36[otu36$sampleTyp=="Root",]
otu36.r1 <- otu36.r[otu36.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, otu36.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu36.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# slope
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu36.r1$TP0
y <- otu36.r1$value
cor.test(x, y, method = "spearman")


otu36.r2 <- otu36.r[otu36.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu36.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu36.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu36.r2$TP0
y <- otu36.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.36r <- ggplot(data=otu36.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(-1,26)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color = "black"),
        axis.title=element_text(size=20,face="bold",color = "black")) +
  annotate("text", x=0.5, y=24, label="Root_Control: R = -0.682, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=21, label="Root_N-addition: R = -0.417, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.36r



root <- p.36r1+p.36r2+p.36r
root
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-otu36.root.pdf",root,height = 4.5, width = 13.5)


###OTU18###
otu18 <- b.otu[b.otu$OTU=="OTU18",]
otu18.r <- otu18[otu18$sampleTyp=="Root",]

root_otu1 <- otu18.r[otu18.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu18.r$value, otu18.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- otu18.r[otu18.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 0.04761749, df = 1, p-value = 0.8272621

p_gs3 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "",y = "Root - OTU18") + 
  ylim(0, 16)+
  geom_text(x=1.5,y=15,label=expression(paste(chi^2, " = ", 0.048, ", ", "Df = ", 1, ", ", "P = ", 0.827, sep="")),
            size=7,color="black",fontface="bold")+
  annotate("text",x=1,y=12,label="a",size=7,color="black")+
  annotate("text",x=2,y=13,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=19,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p_gs3


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_gs1 <- ggplot(otu18.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "",y = "") + 
  ylim(0, 16)+
  annotate("text", x=1, y=15, label="OTU18", 
           size=7, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=13, label="Genome Size: 8.817 Mb", 
           size=7, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=11, label="RRN: 2", 
           size=7, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color="black"),
        axis.title=element_text(size=18,face="bold",color="black"))

p_gs1

##OTU18.lm
otu18 <- b.otu[b.otu$OTU=="OTU18",]
otu18.r <- otu18[otu18$sampleTyp=="Root",]
otu18.r1 <- otu18.r[otu18.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, otu18.r1)))

lm_model1 <- lm(value~TP0, otu18.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu18.r1$TP0
y <- otu18.r1$value

cor.test(x, y, method = "spearman")



otu18.r2 <- otu18.r[otu18.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu18.r2)))
lm_model1 <- lm(value~TP0, otu18.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu18.r2$TP0
y <- otu18.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.18r <- ggplot(data=otu18.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,16)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color = "black"),
        axis.title=element_text(size=18,face="bold",color="black")) +
  annotate("text", x=0.5, y=15, label="Root_Control: R = 0.914, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=13, label="Root_N-addition: R = 0.702, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.18r



root.18 <- p_gs3+p_gs1+p.18r
root.18

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-otu18.root.pdf",root.18,height = 4.5, width = 13.5)

root <- p_gs3+p_gs1+p.18r+p.36r1+p.36r2+p.36r
root
setwd("~/Documents/sunR/Bacteria882/n.treatment/Fig2修改")
#ggsave("B-otu18.otu36.root.revised.pdf",root,height = 9, width = 13.5)


##Figure.3.Genomic.traits identity>95% & covs>95%####
####gs.rrn.compartment.boxplot.no-N####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)

#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/GTDB")
#
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

#otu table
otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")

#blast
blast <- read_xlsx("blast.uotus.bacteria.gtdb.xlsx", sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)


#GTDB metadata
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

#
blast1 <- blast0[blast0$qseqid %in% otu_table$OTU_id,]

blast1 <- blast1 %>% separate(sseqid, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$qseqid),] 


#metadata
blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")

#
otu_table2 <- merge(otu_table,blast2,by.x = "OTU_id",by.y = "qseqid")

otu_table3 <- otu_table2
otu_table3[,2:826] <- apply(otu_table3[,2:826],2,function(x) x/sum(x))  
colSums(otu_table3[,2:826]) 

##
name = "genome_size" 
for(i in 2:826){  
  otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
}

#
genomesize_ave <- as.data.frame(colSums(otu_table3[,2:826]))
genomesize_ave$sample.id <- rownames(genomesize_ave)

names(genomesize_ave) <- c("value","sample.id")
genomesize <- merge(group0, genomesize_ave, by = "sample.id")
#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#write.csv(genomesize,"Bact_genomesize.csv")

genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","Treatment","value")]
genomesize0 <- na.omit(genomesize0)
genomesize0 <- genomesize0 %>% filter(Treatment=="Control")

gs_root <- genomesize0 %>% filter(sampleType=="Root")
gs_leaf <- genomesize0 %>% filter(sampleType=="Leaf")
gs_air <- genomesize0 %>% filter(sampleType=="Air")

shapiro.test(residuals(lm(value~TP0, gs_root)))
hist(residuals(lm(value~TP0, gs_root)))

shapiro.test(residuals(lm(value~TP0, gs_leaf)))
hist(residuals(lm(value~TP0, gs_leaf)))

shapiro.test(residuals(lm(value~TP0, gs_air)))
hist(residuals(lm(value~TP0, gs_air)))


kw <- kruskal(genomesize0$value, genomesize0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 357.7459, df = 2, p-value = 0
plot(kw)
library(dunn.test)
dunn.test(genomesize0$value, genomesize0$sampleType, method = "bonferroni")


pa <- ggplot(genomesize0, aes(x = sampleType, y = value/1000000)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Average Genome Size (Mb)") + 
  ylim(2, 8.3)+
  geom_text(x=2,y=2.6,label=expression(paste(chi^2, " = ", 357.746, ", ", "Df = ", 2, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=6.6,label="b",size=7,color="black")+
  annotate("text",x=2,y=6.9,label="c",size=7,color="black")+
  annotate("text",x=3,y=7.5,label="a",size=7,color="black")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))

pa

###rrnDB
setwd("~/Documents/sunR/Bacteria882/rrndb")
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
# 
otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")
# 
blast <- read_xlsx("blast.uotus.bacteria.rrndb.xlsx",sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)

#rrndb metadata
metadata <- read.delim("rrnDB-5.8.tsv",sep = "\t",header = T)

blast1 <- blast0[blast0$qseqid %in% otu_table$OTU_id,] ##5182 OTU
blast2 <- blast1[!duplicated(blast1$qseqid),]  
blast0 <- concat.split(data = blast2,split.col = "sseqid",sep = "\\|",drop = TRUE) 
#
blast3 <- merge(blast0, metadata, by.x = "sseqid_2", by.y="Data.source.record.id")

# 
otu_table2 <- merge(otu_table, blast3, by.x = "OTU_id", by.y = "qseqid")

otu_table3 <- otu_table2
otu_table3[,2:826] <- apply(otu_table3[,2:826],2,function(x) x/sum(x))  
colSums(otu_table3[,2:826])  

##RRN
name = "X16S.gene.count" #gc_percentage
for(i in 2:826){  
  otu_table3[,i] <- otu_table3[,i]*otu_table3[,name] ##function name?
}

# 
rrna_ave <- as.data.frame(colSums(otu_table3[,2:826]))
rrna_ave$sample.id <- rownames(rrna_ave)

names(rrna_ave) <- c("value","sample.id")
rrna <- merge(group0, rrna_ave, by = "sample.id")
#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#write.csv(rrna,"Bact_rrn.csv")

rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","Treatment","value")]
rrna0 <- na.omit(rrna0)

rrna0 <- rrna0 %>% filter(Treatment=="Control") 

rrn_root <- rrna0 %>% filter(sampleType=="Root")
rrn_leaf <- rrna0 %>% filter(sampleType=="Leaf")
rrn_air <- rrna0 %>% filter(sampleType=="Air")

shapiro.test(residuals(lm(value~TP0, rrn_root)))
hist(residuals(lm(value~TP0, rrn_root)))

shapiro.test(residuals(lm(value~TP0, rrn_leaf)))
hist(residuals(lm(value~TP0, rrn_leaf)))

shapiro.test(residuals(lm(value~TP0, rrn_air)))
hist(residuals(lm(value~TP0, rrn_air)))


kw <- kruskal(rrna0$value, rrna0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 60.61102, df = 2, p-value = 6.894485e-14
#Air,b; leaf,a; Root,b

plot(kw)

pe <- ggplot(rrna0, aes(x = sampleType, y = value)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Average rRNA Operon Copy Number")+
  ylim(2,7)+
  geom_text(x=2,y=2.3,label=expression(paste(chi^2, " = ", 60.611, ", ", "Df = ", 2, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=6,label="b",size=7,color="black")+
  annotate("text",x=2,y=6.6,label="a",size=7,color="black")+
  annotate("text",x=3,y=5,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))

pe


#setwd("~/Documents/sunR/Bacteria882/results2")
#ggsave("B-genome.compartment.boxplot.gs.rrn.pdf",p_comp,height = 4.5, width = 9)


####gs.rrn.temporal.no-N####
genomesize0 <- genomesize[, c("sample.id","sampleType","Treatment","sampleSite","TP0","TP00","value")]
genomesize0 <- na.omit(genomesize0)
genomesize0 <- genomesize0 %>% filter(Treatment=="Control")

model <- lmer(value/1000000~TP0*sampleType + (1|sampleSite), data = genomesize0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

gs_root <- genomesize0 %>% filter(sampleType=="Root")
shapiro.test(residuals(lm(value~TP0, gs_root)))
cor.test(gs_root$value/1000000, gs_root$TP0, method = "spearman")

gs_leaf <- genomesize0 %>% filter(sampleType=="Leaf")
shapiro.test(residuals(lm(value~TP0, gs_leaf)))
cor.test(gs_leaf$value/1000000, gs_leaf$TP0, method = "spearman")

gs_air <- genomesize0[genomesize0$sampleType == "Air",]
shapiro.test(residuals(lm(value~TP0, gs_air)))
cor.test(gs_air$value/1000000, gs_air$TP0, method = "spearman")


pb <- ggplot(data=genomesize0) + 
  geom_jitter(aes(x= TP0, y = value/1000000, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= value/1000000, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(2, 8.3)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=1, y=3.4, label="Air: R = 0.005, P = 0.924", 
           size=5, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1, y=2.9, label="Leaf: R = 0.174, P = 0.032", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1, y=2.4, label="Root: R = 0.607, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") 


pb


###rrnDB
rrna0 <- rrna[, c("sample.id","sampleType","sampleSite","TP0","TP00","Treatment","value")]
rrna0 <- na.omit(rrna0)
rrna0 <- rrna0 %>% filter(Treatment=="Control")

model <- lmer(value~TP0*sampleType + (1|sampleSite), data = rrna0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)  
shapiro.test(res)


rrn_root <- rrna0[rrna0$sampleType == "Root",]
shapiro.test(residuals(lm(value~TP0, rrn_root)))
cor.test(rrn_root$value, rrn_root$TP0, method = "spearman")

rrn_leaf <- rrna0[rrna0$sampleType == "Leaf",]
shapiro.test(residuals(lm(value~TP0, rrn_leaf)))
cor.test(rrn_leaf$value, rrn_leaf$TP0, method = "spearman")

rrn_air <- rrna0[rrna0$sampleType == "Air",]
shapiro.test(residuals(lm(value~TP0, rrn_air)))
cor.test(rrn_air$value, rrn_air$TP0, method = "spearman")


#"#1F71E0", "#0F7C41","#7F4E49"
pf <- ggplot() + 
  geom_jitter(data=rrna0,  height = 0, 
              aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna0,
              aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),
                     guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(2, 7)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=3)),
         shape = guide_legend(title = "Compartment"))+
  #guides(color = guide_legend(override.aes=list(size = 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none",
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))+
  annotate("text", x=1, y=6.7, label="Air: R = -0.002, P = 0.964", 
           size=5, hjust = 0, color="blue",fontface="bold")+
  annotate("text", x=1, y=6.3, label="Leaf: R = 0.200, P = 0.013", 
           size=5, hjust = 0, color="forestgreen",fontface="bold")+
  annotate("text", x=1, y=5.9, label="Root: R = -0.613, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") 

pf


####gs.rrn.leaf.root.Treatment####
genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
genomesize0 <- na.omit(genomesize0)

genomesize1 <- genomesize0 %>% filter(sampleType=="Root")

model <- lmer(value/1000000~TP0*Treatment + (1|sampleSite), data = genomesize1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

root <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="Control")
shapiro.test(residuals(lm(value/1000000~TP0, root)))
cor.test(root$value/1000000, root$TP0, method = "spearman")

root <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
shapiro.test(residuals(lm(value/1000000~TP0, root)))
cor.test(root$value/1000000, root$TP0, method = "spearman")


pd <- ggplot() + 
  geom_jitter(data=genomesize1,  height = 0, 
              aes(x= TP0, y = value/1000000, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=genomesize1,
              aes(x= TP0, y= value/1000000, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,8.3)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=3.4, label="Root_Control: R = 0.607, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=2.9, label="Root_N-addition: R = 0.457, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


pd



genomesize2 <- genomesize0 %>% filter(sampleType=="Leaf")
model <- lmer(value/1000000~TP0*Treatment + (1|sampleSite), data = genomesize2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)   
shapiro.test(res)

leaf <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
shapiro.test(residuals(lm(value/1000000~TP0, leaf)))
cor.test(leaf$value/1000000, leaf$TP0, method = "spearman")

leaf <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
shapiro.test(residuals(lm(value/1000000~TP0, leaf)))
cor.test(leaf$value/1000000, leaf$TP0, method = "spearman")

pc <- ggplot() + 
  geom_jitter(data=genomesize2,  height = 0, 
              aes(x= TP0, y = value/1000000, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=genomesize2,
              aes(x= TP0, y= value/1000000, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("forestgreen","black"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,8.3)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=3.4, label="Leaf_Control: R = 0.174, P = 0.032", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=0.5, y=2.9, label="Leaf_N-addition: R = 0.008, P = 0.953", 
           size=5, hjust = 0, color="black",fontface="bold") 


pc



###rrnDB
rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
rrna0 <- na.omit(rrna0)

rrna1 <- rrna0 %>% filter(sampleType=="Root")

model <- lmer(value~TP0*Treatment + (1|sampleSite), data = rrna1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

root <- rrna0 %>% filter(sampleType=="Root"&Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, root)))
cor.test(root$value, root$TP0, method = "spearman")

root <- rrna0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, root)))
cor.test(root$value, root$TP0, method = "spearman")


ph <- ggplot() + 
  geom_jitter(data=rrna1,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna1,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,7)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=6.7, label="Root_Control: R = -0.613, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=6.3, label="Root_N-addition: R = -0.118, P = 0.332", 
           size=5, hjust = 0, color="black",fontface="bold") 


ph



rrna2 <- rrna0 %>% filter(sampleType=="Leaf")

model <- lmer(value~TP0*Treatment + (1|sampleSite), data = rrna2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

leaf <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, leaf)))
cor.test(leaf$value, leaf$TP0, method = "spearman")

leaf <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, leaf)))
cor.test(leaf$value, leaf$TP0, method = "spearman")

pg <- ggplot() + 
  geom_jitter(data=rrna2,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna2,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("forestgreen","black"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,7)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=6.7, label="Leaf_Control: R = 0.200, P = 0.013", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=0.5, y=6.3, label="Leaf_N-addition: R = 0.341, P = 0.012", 
           size=5, hjust = 0, color="black",fontface="bold") 


pg

#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")

library(gridExtra)

p <- grid.arrange(
  pa,pb,pc,pd,
  pe,pf,pg,ph,
  ncol = 4,    
  widths = c(1, 1, 1,1) 
)


#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")

#ggsave("B.otu-genomic.ident0.95&covs0.95.pdf",p,height = 9, width = 18)

##Figure4.nifH####
###nifH.gene.qpcr####
setwd("~/Documents/sunR/Bacteria882/n.treatment/nifH")
data_nif <- read_xlsx("nifH.data2.xlsx",sheet=5)
data_qubit <- read_xlsx("DNA-qubit.xlsx")
data1 <- merge(data_nif,data_qubit,by="sample.id")
data1$value <- data1$copy/data1$qubit

group <- read_xlsx("../../sample_metadata.xlsx", sheet = 3)
group1 <- group %>% filter(sampleType=="Root")
group1 <- group1 %>% dplyr::select("sample.id","sampleType","TP0","TP00","sampleSite","Treatment")

data2 <- merge(data1,group1,by="sample.id")
data3 <- data2 %>% dplyr::select("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")

model <- lmer(value ~ TP0*Treatment + (1|sampleSite), data = data3)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="nifH - res")
qqnorm(res); qqline(res)    
par(mfrow=c(1,2))

shapiro.test(res)

root1 <- data3 %>% filter(Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, root1)))
cor.test(root1$value, root1$TP0, method = "spearman")

model <- lmer(value ~ TP0 + (1|sampleSite), data = root1)
summary(model)
anova(model)

root2 <- data3 %>% filter(Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, root2)))
cor.test(root2$value, root2$TP0, method = "spearman")

model <- lmer(value ~ TP0 + (1|sampleSite), data = root2)
summary(model)
anova(model)

pd <- ggplot() + 
  geom_jitter(data=data3,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=data3,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "nifH copy number", title = "Based on qPCR")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  #scale_y_continuous(limits = c(2,8.3)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0, size = 15, face = "bold",color = "black"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=17,face="bold",color = "black")) +
  annotate("text", x=0.5, y=720000, label="Root_Control: R = 0.116, P = 0.109", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=640000, label="Root_N-addition: R = -0.386, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


pd

#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")

#ggsave("B-nifH.pdf",pd,height = 5, width = 5.5)

###nifH.Picrust2####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

group2 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group2)
group3 <- group2 %>% filter(Treatment == "Control")


ko <- read_xlsx("ko_gene.xlsx", sheet = 1)
ko_gene <- read_xlsx("ko_gene.xlsx", sheet=2)
ko_ncycle <- read_xlsx("ko_gene.xlsx",sheet=3)

ko1 <- ko[,colnames(ko) %in% c("KO",group2$sample.id)]
ko1[, 2:826]  <- apply(ko1[, 2:826], 2, function(x) x / sum(x))  

ko2 <- ko1[ko1$KO %in% ko_ncycle$KO,]
data1 <- merge(ko2,ko_ncycle,by="KO")
rownames(data1) <- data1$gene
data1 <- data1[,!(names(data1) %in% c("KO","gene","KO_name","EC"))]

data10 <- data.frame(t(data1))
data10$sample.id <- rownames(data10)
data11 <- merge(data10,group2,by="sample.id")

data12 <- gather(data11,key = "Ngene",value = "value", -"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")

nifh <- data12 %>% filter(Ngene=="nifH"&sampleType=="Root")

root1 <- nifh %>% filter(Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, root1)))
cor.test(root1$value, root1$TP0, method = "pearson")


root2 <- nifh %>% filter(Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, root2)))
cor.test(root2$value, root2$TP0, method = "spearman")

pd <- ggplot() + 
  geom_jitter(data=nifh,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=nifh,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "Relative abundance of nifH",title = "Based on PICRUSt2")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  #scale_y_continuous(limits = c(2,8.3)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0, size = 15, face = "bold",color = "black"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=17,face="bold",color = "black")) +
  annotate("text", x=0.5, y=.0001, label="Root_Control: R = 0.397, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=.00009, label="Root_N-addition: R = 0.607, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


pd
#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#ggsave("B-nifH(Picrust2).pdf",pd,height = 5, width = 5.5)




##Supplementary.Figure#################
##Supplementary.Figure1.Alpha.diversity-samples825####
library(vegan)
library(picante)
library(ggplot2)
library(car)
library(ggpubr)
library(splitstackshape)
library(dplyr)
library(stringr)
library(broom)
library(agricolae)  
library(readxl)
library(patchwork)
##Bacteria-air7, root, leaf
rm(list = ls())
setwd("~/Documents/sunR/Bacteria882")
group <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group[order(group$sample.id), ] 
group1 <- group[, c("sample.id","TP00","TP0","sampleType","sampleSite","Treatment")]
group2 <- group1 %>% filter(Treatment=="Control")

otu_all0 <- read.csv("otu_bactFlattening825.csv", check.names = F, row.names = 1, header = T)
otu_all <- otu_all0 %>% dplyr::select(all_of(group1$sample.id))

otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_all1$OTU.ID <- rownames(otu_all1)
otu_all1 <- otu_all1[, 1:825]  
otu_all1 <- t(otu_all1)

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = "shannon", base = base)
  Simpson <- diversity(x, index = "simpson") # Gini-Simpson 
  Pielou <- Shannon / log(ACE, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- "PD_whole_tree"
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}

alpha_all <- alpha(otu_all1, base = 2)
alpha_all$sample.id <- rownames(alpha_all)

sum(alpha_all$sample.id %in% group1$sample.id == FALSE) 
which(alpha_all$sample.id %in% group1$sample.id == FALSE)

alpha_all0 <- merge(alpha_all, group1, by = "sample.id", all.x = TRUE) 
#write.csv(alpha_all0, "Bacteria_alpha.csv")
alpha_all0$t <- paste(alpha_all0$sampleType,alpha_all0$Treatment,sep="_")

####alpha.diversity.Shannon.temporal####
alpha_root <- alpha_all0[alpha_all0$t == "Root_Control", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_root)))
cor.test(alpha_root$TP0, alpha_root$Shannon, method = "spearman")

alpha_root <- alpha_all0[alpha_all0$t == "Root_N-addition", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_root)))
cor.test(alpha_root$TP0, alpha_root$Shannon, method = "spearman")

alpha_leaf <- alpha_all0[alpha_all0$t == "Leaf_Control", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_leaf)))
cor.test(alpha_leaf$TP0, alpha_leaf$Shannon, method = "spearman")

alpha_leaf <- alpha_all0[alpha_all0$t == "Leaf_N-addition", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_leaf)))
cor.test(alpha_leaf$TP0, alpha_leaf$Shannon, method = "spearman")

alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Shannon
cor.test(x, y, method = "spearman")

library(lmerTest)
model <- lmer(Shannon ~ TP0*sampleType*Treatment + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)


res <- resid(model, type = "pearson")
fit <- fitted(model)
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals")
qqnorm(res); qqline(res)                        
par(mfrow=c(1,1))
shapiro.test(res)


p.sh <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Shannon, color = t, shape = t), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Shannon, group = t, color = t), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16,16,17, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Time")+
  labs(x = "Time",y = "Shannon index")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0.5, 10)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=3.2, label="Air: R = -0.334, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=2.6, label="Leaf_C: R = 0.194, P = 0.016", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=2, label="Leaf_N: R = -0.137, P = 0.324", 
           size=6, hjust = 0, color="darkolivegreen",fontface="bold") +
  annotate("text", x=1.5, y=1.4, label="Root_C: R = 0.562, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=1.5, y=0.8, label="Root_N: R = 0.005, P = 0.970", 
           size=6, hjust = 0, color="sienna",fontface="bold") 


p.sh


#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.shannon.pdf",p1,height = 5, width = 5)

kw <- kruskal(alpha_all0$Shannon, alpha_all0$t, p.adj = "bonferroni")
kw
plot(kw)

p_shannon <- ggplot(alpha_all0, aes(x = t, y = Shannon)) + 
  geom_boxplot(aes(fill = t),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = t),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Compartment") +
  labs(x = "Group",y = "Shannon index") + 
  geom_text(x=3,y=2,label=expression(paste(chi^2, " = ", 495.130, ", ", "Df = ", 4, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=9.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=8,label="c",size=7,color="black")+
  annotate("text",x=3,y=8,label="bc",size=7,color="black")+
  annotate("text",x=4,y=8,label="b",size=7,color="black")+
  annotate("text",x=5,y=8,label="b",size=7,color="black")+
  ylim(0.5, 10)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_shannon

####alpha.diversity.Richness.temporal####
alpha_root <- alpha_all0[alpha_all0$t == "Root_Control", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_root)))
cor.test(alpha_root$TP0, alpha_root$Richness, method = "spearman")

alpha_root <- alpha_all0[alpha_all0$t == "Root_N-addition", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_root)))
cor.test(alpha_root$TP0, alpha_root$Richness, method = "spearman")


alpha_leaf <- alpha_all0[alpha_all0$t == "Leaf_Control", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_leaf)))
cor.test(alpha_leaf$TP0, alpha_leaf$Richness, method = "spearman")

alpha_leaf <- alpha_all0[alpha_all0$t == "Leaf_N-addition", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_leaf)))
cor.test(alpha_leaf$TP0, alpha_leaf$Richness, method = "spearman")


alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Richness
cor.test(x, y, method = "spearman")

model <- lmer(Richness ~ TP0*sampleType*Treatment + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)


res <- resid(model, type = "pearson")
fit <- fitted(model)
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals")
qqnorm(res); qqline(res)                        
par(mfrow=c(1,1))

shapiro.test(res)

p.rich <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Richness, color = t, shape = t), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Richness, group = t, color = t), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 16, 17, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Time")+
  labs(x = "Time",y = "OTU richness")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0,770)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=180, label="Air: R = -0.396, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=135, label="Leaf_C: R = 0.285, P < 0.001", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=90, label="Leaf_N: R = -0.097, P = 0.484", 
           size=6, hjust = 0, color="darkolivegreen",fontface="bold") +
  annotate("text", x=1.5, y=45, label="Root_C: R = 0.697, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=1.5, y=1, label="Root_N: R = 0.277, P = 0.020", 
           size=6, hjust = 0, color="sienna",fontface="bold") 


p.rich

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.richness.pdf",p2,height = 5, width = 5)


kw <- kruskal(alpha_all0$Richness, alpha_all0$t, p.adj = "bonferroni")
kw
plot(kw)

p_rich <- ggplot(alpha_all0, aes(x = t, y = Richness)) + 
  geom_boxplot(aes(fill = t),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = t),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkolivegreen","darkred","sienna"), name="Compartment") +
  labs(x = "Group",y = "OTU richness") + 
  geom_text(x=3,y=20,label=expression(paste(chi^2, " = ", 611.376, ", ", "Df = ", 4, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=760,label="a",size=7,color="black")+
  annotate("text",x=2,y=470,label="c",size=7,color="black")+
  annotate("text",x=3,y=470,label="b",size=7,color="black")+
  annotate("text",x=4,y=470,label="b",size=7,color="black")+
  annotate("text",x=5,y=470,label="b",size=7,color="black")+
  ylim(0, 770)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_rich


combined_plot <- p.sh + p_shannon + p.rich + p_rich +
  plot_layout(ncol = 4) 
combined_plot
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-alpha.shannon.richness.825.revised.pdf",combined_plot,height = 5, width = 20)




##Supplementary Figure.3.Treatment.Root.Bradyrhizobium & Herbaspirillum####
library(vegan)
library(betapart)
library(colorRamps)
library(agricolae)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(dplyr)
#rm(list = ls())
setwd("~/Documents/sunR/Bacteria882")

##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

###genus####
bact.genus <- aggregate(otu_all1,by=list(ID$genus) , sum) 
rownames(bact.genus) <- bact.genus[,1]; 
bact.genus <- bact.genus[,-1] 
data6 <- data.frame(t(bact.genus))

total6 <- apply(data6, 1, sum); 
bact.relabu.genus <- data.frame(lapply(data6, function(x) {  x / total6  }) )
bact.genus0 <- bact.relabu.genus
bact.genus0 <- bact.genus0[,order(-colSums(bact.genus0))]

colnames(bact.genus0)
bact.genus15 <- bact.genus0[, c(1:15)]
bact.genus15$sample.id <- rownames(bact.genus15)

b.genus15 <- merge(bact.genus15, group2, by="sample.id")

b.genus <- gather(b.genus15,key = "genus",value = "Abundance",-"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")


b.genus$value <- b.genus$Abundance*100

###Bradyrhizobium####
gb <- b.genus[b.genus$genus=="Bradyrhizobium",]
gb.r <- gb[gb$sampleTyp=="Root",]

root1 <- gb.r %>% filter(Treatment=="Control")

lm1 <- lm(value~TP0, root1)
shapiro.test(residuals(lm1))
hist(residuals(lm1))

root2 <- gb.r %>% filter(Treatment=="N-addition")

lm2 <- lm(value~TP0, root2)
shapiro.test(residuals(lm2))
hist(residuals(lm2))

kw <- kruskal(gb.r$value, gb.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- gb.r[gb.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw


p.gbr1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "",y = "Root - Bradyrhizobium") + 
  ylim(0, 16)+
  geom_text(x=1.5,y=15,label=expression(paste(chi^2, " = ", 0.137, ", ", "Df = ", 1, ", ", "P = ", 0.711, sep="")),
            size=6.9,color="black",fontface="bold")+
  annotate("text",x=1,y=13,label="a",size=7,color="black")+
  annotate("text",x=2,y=13.5,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=19,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p.gbr1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)
root_otu1 <- gb.r[gb.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

p.gbr2 <- ggplot(gb.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "",y = "") + 
  ylim(0, 15)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p.gbr2

##
gb <- b.genus[b.genus$genus=="Bradyrhizobium",]
gb.r <- gb[gb$sampleTyp=="Root",]
gb.r1 <- gb.r[gb.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, gb.r1)))

x <- gb.r1$TP0
y <- gb.r1$value
cor.test(x, y, method = "spearman")


gb.r2 <- gb.r[gb.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, gb.r2)))

x <- gb.r2$TP0
y <- gb.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.gbr <- ggplot(data=gb.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,15)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color = "black"),
        axis.title=element_text(size=20,face="bold",color = "black")) +
  annotate("text", x=0.5, y=14, label="Root_Control: R = 0.916, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=12.5, label="Root_N-addition: R = 0.705, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.gbr



root <- p.gbr1+p.gbr2+p.gbr
root
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-root.Bradyrhizobium.pdf",root,height = 4.5, width = 13.5)


###Herbaspirillum####
gh <- b.genus[b.genus$genus=="Herbaspirillum",]
gh.r <- gh[gh$sampleTyp=="Root",]

kw <- kruskal(gh.r$value, gh.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- gh.r[gh.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw


p_ghr1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - Herbaspirillum") + 
  ylim(-1, 27)+
  geom_text(x=1.5,y=24,label=expression(paste(chi^2, " = ", 22.717, ", ", "Df = ", 1, ", ", "P < ", 0.001, sep="")),
            size=7,color="black",fontface="bold")+
  annotate("text",x=1,y=12,label="a",size=7,color="black")+
  annotate("text",x=2,y=8,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=19,color="black"),
        axis.title=element_text(size=20,face="bold",color="black"))

p_ghr1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)

root_otu1 <- gh.r[gh.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw


p_ghr2 <- ggplot(gh.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(-1, 27)+
  annotate("text",x=2,y=24,label="*",size=12,color="black")+
  annotate("text",x=5,y=5,label="*",size=12,color="black")+
  annotate("text",x=6,y=5,label="*",size=12,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color="black"),
        axis.title=element_text(size=18,face="bold",color="black"))

p_ghr2

##
gh <- b.genus[b.genus$genus=="Herbaspirillum",]
gh.r <- gh[gh$sampleTyp=="Root",]
gh.r1 <- gh.r[gh.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, gh.r1)))

x <- gh.r1$TP0
y <- gh.r1$value

cor.test(x, y, method = "spearman")



gh.r2 <- gh.r[gh.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, gh.r2)))

x <- gh.r2$TP0
y <- gh.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.ghr <- ggplot(data=gh.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(-1,27)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=18,color = "black"),
        axis.title=element_text(size=18,face="bold",color="black")) +
  annotate("text", x=0.5, y=24, label="Root_Control: R = -0.663, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=21, label="Root_N-addition: R = -0.409, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.ghr



root.gh <- p_ghr1+p_ghr2+p.ghr
root.gh

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-otu18.root.pdf",root.gh,height = 4.5, width = 13.5)

root <- p.gbr1+p.gbr2+p.gbr+p_ghr1+p_ghr2+p.ghr

root

setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-genus.B&H.root.pdf",root,height = 9, width = 14)




##Supplementary.Figure4.ASV.Genomic.traits identity>95% & covs>95%####
####gs.rrn.compartment.boxplot.no-N####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)

#rm(list=ls())
setwd("~/Documents/sunR/bacteria.ASV")
# 
group0 <- read_xlsx("group.asv.xlsx", sheet = 1)
ID <- read_xlsx("ID_bacteria.xlsx",sheet=2)
#
otu_table <- otu_table <- read.csv("asv_Flattening811.csv")
colnames(otu_table)[1] <- "asv_id"

asv <- merge(otu_table,ID, by="asv_id")
# 
setwd("~/Documents/sunR/bacteria.ASV/genomic")
blast <- read_xlsx("blast.asv.bacteria.gtdb.xlsx", sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)


# 
setwd("~/Documents/sunR/Bacteria882/GTDB")
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

#
blast1 <- blast0[blast0$qseqid %in% asv$qseqid,]
blast1 <- merge(blast1, ID, by="qseqid")

blast1 <- blast1 %>% separate(sseqid, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$asv_id),] 

blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")

otu_table2 <- merge(otu_table,blast2,by = "asv_id")

otu_table3 <- otu_table2
otu_table3[,2:812] <- apply(otu_table3[,2:812],2,function(x) x/sum(x))  
colSums(otu_table3[,2:812])  


##average genome size 
name = "genome_size" 
for(i in 2:812){  
  otu_table3[,i] <- otu_table3[,i]*otu_table3[,name]
}

genomesize_ave <- as.data.frame(colSums(otu_table3[,2:812]))
genomesize_ave$sample.id <- rownames(genomesize_ave)

names(genomesize_ave) <- c("value","sample.id")
genomesize <- merge(group0, genomesize_ave, by = "sample.id")
#setwd("~/Documents/sunR/bacteria.ASV/genomic")
#write.csv(genomesize,"Bact_genomesize(identity>95%).csv")

genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","Treatment","value")]
genomesize0 <- na.omit(genomesize0)
genomesize0 <- genomesize0 %>% filter(Treatment=="Control")

gs_root <- genomesize0 %>% filter(sampleType=="Root")
gs_leaf <- genomesize0 %>% filter(sampleType=="Leaf")
gs_air <- genomesize0 %>% filter(sampleType=="Air")

shapiro.test(residuals(lm(value~TP0, gs_root)))
hist(residuals(lm(value~TP0, gs_root)))

shapiro.test(residuals(lm(value~TP0, gs_leaf)))
hist(residuals(lm(value~TP0, gs_leaf)))

shapiro.test(residuals(lm(value~TP0, gs_air)))
hist(residuals(lm(value~TP0, gs_air)))


kw <- kruskal(genomesize0$value, genomesize0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 357.7459, df = 2, p-value = 0
plot(kw)
library(dunn.test)
dunn.test(genomesize0$value, genomesize0$sampleType, method = "bonferroni")


pa <- ggplot(genomesize0, aes(x = sampleType, y = value/1000000)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Average Genome Size (Mb)") + 
  ylim(2, 8.5)+
  #geom_text(x=2,y=2.5,label=expression(~chi^2==414.677~~Df==2~~P==0),size=6,color="black")+
  geom_text(x=2,y=2.6,label=expression(paste(chi^2, " = ", 299.719, ", ", "Df = ", 2, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=6.6,label="b",size=7,color="black")+
  annotate("text",x=2,y=6.9,label="c",size=7,color="black")+
  annotate("text",x=3,y=7.5,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))

pa

###rrnDB
setwd("~/Documents/sunR/bacteria.ASV/genomic")

group0 <- read_xlsx("../group.asv.xlsx", sheet = 1)
ID <- read_xlsx("../ID_bacteria.xlsx",sheet=2)
otu_table <- otu_table <- read.csv("../asv_Flattening811.csv")
colnames(otu_table)[1] <- "asv_id"

asv <- merge(otu_table,ID, by="asv_id")
setwd("~/Documents/sunR/bacteria.ASV/genomic")
blast <- read_xlsx("blast.asv.bacteria.rrndb.xlsx", sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)

setwd("~/Documents/sunR/Bacteria882/rrndb")
metadata <- read.delim("rrnDB-5.8.tsv",sep = "\t",header = T)
blast1 <- blast0[blast0$qseqid %in% asv$qseqid,]
blast1 <- merge(blast1, ID, by="qseqid")
blast1 <- blast1[!duplicated(blast1$asv_id),] 

blast2 <- concat.split(data = blast1,split.col = "sseqid",sep = "\\|",drop = TRUE)
blast3 <- merge(blast2, metadata, by.x = "sseqid_2", by.y="Data.source.record.id")
otu_table2 <- merge(otu_table,blast3,by = "asv_id")

otu_table3 <- otu_table2
otu_table3[,2:812] <- apply(otu_table3[,2:812],2,function(x) x/sum(x))  
colSums(otu_table3[,2:812]) 


name = "X16S.gene.count" #gc_percentage
for(i in 2:812){  
  otu_table3[,i] <- otu_table3[,i]*otu_table3[,name] ##function name?
}


rrna_ave <- as.data.frame(colSums(otu_table3[,2:812]))
rrna_ave$sample.id <- rownames(rrna_ave)

names(rrna_ave) <- c("value","sample.id")
rrna <- merge(group0, rrna_ave, by = "sample.id")
#setwd("~/Documents/sunR/bacteria.ASV/genomic")
#write.csv(rrna,"Bact_rrn(identity>95%).csv")

rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","Treatment","value")]
rrna0 <- na.omit(rrna0)

rrna0 <- rrna0 %>% filter(Treatment=="Control") 

rrn_root <- rrna0 %>% filter(sampleType=="Root")
rrn_leaf <- rrna0 %>% filter(sampleType=="Leaf")
rrn_air <- rrna0 %>% filter(sampleType=="Air")

shapiro.test(residuals(lm(value~TP0, rrn_root)))
hist(residuals(lm(value~TP0, rrn_root)))

shapiro.test(residuals(lm(value~TP0, rrn_leaf)))
hist(residuals(lm(value~TP0, rrn_leaf)))

shapiro.test(residuals(lm(value~TP0, rrn_air)))
hist(residuals(lm(value~TP0, rrn_air)))


kw <- kruskal(rrna0$value, rrna0$sampleType, p.adj = "bonferroni")
kw

plot(kw)

pe <- ggplot(rrna0, aes(x = sampleType, y = value)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "forestgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Average rRNA Operon Copy Number")+
  ylim(2,7)+
  geom_text(x=2,y=2.1,label=expression(paste(chi^2, " = ", 60.611, ", ", "Df = ", 2, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=6,label="b",size=7,color="black")+
  annotate("text",x=2,y=6.6,label="a",size=7,color="black")+
  annotate("text",x=3,y=5,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))

pe


#setwd("~/Documents/sunR/Bacteria882/results2")
#ggsave("B-genome.compartment.boxplot.gs.rrn.pdf",p_comp,height = 4.5, width = 9)


####gs.rrn.temporal.no-N####
genomesize0 <- genomesize[, c("sample.id","sampleType","Treatment","sampleSite","TP0","TP00","value")]
genomesize0 <- na.omit(genomesize0)
genomesize0 <- genomesize0 %>% filter(Treatment=="Control")

model <- lmer(value/1000000~TP0*sampleType + (1|sampleSite), data = genomesize0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

gs_root <- genomesize0 %>% filter(sampleType=="Root")
shapiro.test(residuals(lm(value~TP0, gs_root)))
cor.test(gs_root$value/1000000, gs_root$TP0, method = "spearman")

gs_leaf <- genomesize0 %>% filter(sampleType=="Leaf")
shapiro.test(residuals(lm(value~TP0, gs_leaf)))
cor.test(gs_leaf$value/1000000, gs_leaf$TP0, method = "spearman")

gs_air <- genomesize0[genomesize0$sampleType == "Air",]
shapiro.test(residuals(lm(value~TP0, gs_air)))
cor.test(gs_air$value/1000000, gs_air$TP0, method = "spearman")


pb <- ggplot(data=genomesize0) + 
  geom_jitter(aes(x= TP0, y = value/1000000, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= value/1000000, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(2, 8.5)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=1, y=3.4, label="Air: R = -0.060, P = 0.257", 
           size=5, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1, y=2.9, label="Leaf: R = 0.073, P = 0.389", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1, y=2.4, label="Root: R = 0.681, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") 


pb


###rrnDB
rrna0 <- rrna[, c("sample.id","sampleType","sampleSite","TP0","TP00","Treatment","value")]
rrna0 <- na.omit(rrna0)
rrna0 <- rrna0 %>% filter(Treatment=="Control")

model <- lmer(value~TP0*sampleType + (1|sampleSite), data = rrna0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)   
shapiro.test(res)


rrn_root <- rrna0[rrna0$sampleType == "Root",]
shapiro.test(residuals(lm(value~TP0, rrn_root)))
cor.test(rrn_root$value, rrn_root$TP0, method = "spearman")

rrn_leaf <- rrna0[rrna0$sampleType == "Leaf",]
shapiro.test(residuals(lm(value~TP0, rrn_leaf)))
cor.test(rrn_leaf$value, rrn_leaf$TP0, method = "spearman")

rrn_air <- rrna0[rrna0$sampleType == "Air",]
shapiro.test(residuals(lm(value~TP0, rrn_air)))
cor.test(rrn_air$value, rrn_air$TP0, method = "spearman")


#"#1F71E0", "#0F7C41","#7F4E49"
pf <- ggplot() + 
  geom_jitter(data=rrna0,  height = 0, 
              aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna0,
              aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),
                     guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(2, 7)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=3)),
         shape = guide_legend(title = "Compartment"))+
  #guides(color = guide_legend(override.aes=list(size = 3)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none",
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold"))+
  annotate("text", x=1, y=6.7, label="Air: R = 0.089, P = 0.093", 
           size=5, hjust = 0, color="blue",fontface="bold")+
  annotate("text", x=1, y=6.3, label="Leaf: R = 0.113, P = 0.181", 
           size=5, hjust = 0, color="forestgreen",fontface="bold")+
  annotate("text", x=1, y=5.9, label="Root: R = -0.647, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") 

pf


####gs.rrn.leaf.root.Treatment####
genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
genomesize0 <- na.omit(genomesize0)

genomesize1 <- genomesize0 %>% filter(sampleType=="Root")

model <- lmer(value/1000000~TP0*Treatment + (1|sampleSite), data = genomesize1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

root <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="Control")
shapiro.test(residuals(lm(value/1000000~TP0, root)))
cor.test(root$value/1000000, root$TP0, method = "spearman")

root <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
shapiro.test(residuals(lm(value/1000000~TP0, root)))
cor.test(root$value/1000000, root$TP0, method = "spearman")


pd <- ggplot() + 
  geom_jitter(data=genomesize1,  height = 0, 
              aes(x= TP0, y = value/1000000, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=genomesize1,
              aes(x= TP0, y= value/1000000, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,8.5)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=3.4, label="Root_Control: R = 0.681, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=2.9, label="Root_N-addition: R = 0.564, P = 3.75e-07", 
           size=5, hjust = 0, color="black",fontface="bold") 


pd



genomesize2 <- genomesize0 %>% filter(sampleType=="Leaf")
model <- lmer(value/1000000~TP0*Treatment + (1|sampleSite), data = genomesize2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="GS - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

leaf <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
shapiro.test(residuals(lm(value/1000000~TP0, leaf)))
cor.test(leaf$value/1000000, leaf$TP0, method = "spearman")

leaf <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
shapiro.test(residuals(lm(value/1000000~TP0, leaf)))
cor.test(leaf$value/1000000, leaf$TP0, method = "spearman")

pc <- ggplot() + 
  geom_jitter(data=genomesize2,  height = 0, 
              aes(x= TP0, y = value/1000000, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=genomesize2,
              aes(x= TP0, y= value/1000000, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("forestgreen","black"), name="Time")+
  labs(x = "Time",y = "Average Genome Size (Mb)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,8.5)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=3.4, label="Leaf_Control: R = 0.073, P = 0.389", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=0.5, y=2.9, label="Leaf_N-addition: R = -0.157, P = 0.276", 
           size=5, hjust = 0, color="black",fontface="bold") 


pc



###rrnDB
rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
rrna0 <- na.omit(rrna0)

rrna1 <- rrna0 %>% filter(sampleType=="Root")

model <- lmer(value~TP0*Treatment + (1|sampleSite), data = rrna1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

root <- rrna0 %>% filter(sampleType=="Root"&Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, root)))
cor.test(root$value, root$TP0, method = "spearman")

root <- rrna0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, root)))
cor.test(root$value, root$TP0, method = "spearman")


ph <- ggplot() + 
  geom_jitter(data=rrna1,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna1,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,7)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=6.7, label="Root_Control: R = -0.647, P < 2.2e-16", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=6.3, label="Root_N-addition: R = -0.181, P = 0.135", 
           size=5, hjust = 0, color="black",fontface="bold") 


ph



rrna2 <- rrna0 %>% filter(sampleType=="Leaf")

model <- lmer(value~TP0*Treatment + (1|sampleSite), data = rrna2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN - res")
qqnorm(res); qqline(res)    
shapiro.test(res)

leaf <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
shapiro.test(residuals(lm(value~TP0, leaf)))
cor.test(leaf$value, leaf$TP0, method = "spearman")

leaf <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
shapiro.test(residuals(lm(value~TP0, leaf)))
cor.test(leaf$value, leaf$TP0, method = "spearman")

pg <- ggplot() + 
  geom_jitter(data=rrna2,  height = 0, 
              aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1) + 
  geom_smooth(data=rrna2,
              aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 1.5, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("forestgreen","black"), name="Time")+
  labs(x = "Time",y = "Average rRNA Operon Copy Number")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(2,7)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=0.5, y=6.7, label="Leaf_Control: R = 0.113, P = 0.181", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=0.5, y=6.3, label="Leaf_N-addition: R = 0.539, P = 5.40e-05", 
           size=5, hjust = 0, color="black",fontface="bold") 


pg

#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")



library(gridExtra)

p <- grid.arrange(
  pa,pb,pc,pd,
  pe,pf,pg,ph,
  ncol = 4,    
  widths = c(1, 1, 1,1) 
)


#setwd("~/Documents/sunR/bacteria.ASV/results")
#ggsave("Bact.asv_genomic.ident0.95&covs0.95.pdf",p,height = 9, width = 18)






##Supplementary.Figure5.Treatment.Root.top10otu####
library(vegan)
library(betapart)
library(colorRamps)
library(agricolae)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(dplyr)
##Bacteria-AIR,ROOT,LEAF-7
setwd("~/Documents/sunR/Bacteria882")
rm(list = ls())
##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
#ID$otu <- paste(ID$names, ID$phylum, ID$genus, sep = ".")

bact.otu <- aggregate(otu_all1,by=list(ID$names) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )
#rowSums(bact.relabu.otu[,1:5211])
bact.otu0 <- bact.relabu.otu
bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))]

colnames(bact.otu0)
bact.otu10 <- bact.otu0[, c(1:10)]
bact.otu10$sample.id <- rownames(bact.otu10)

b.otu10 <- merge(bact.otu10, group2, by="sample.id")

b.otu <- gather(b.otu10,key = "OTU",value = "Abundance",-"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")


b.otu$value <- b.otu$Abundance*100

###OTU36###
otu36 <- b.otu[b.otu$OTU=="OTU36",]
otu36.r <- otu36[otu36$sampleTyp=="Root",]
otu36.r <- otu36.r[otu36.r$value!="0",]
root_otu1 <- otu36.r[otu36.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu36.r$value, otu36.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- otu36.r[otu36.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 23.87689, df = 1, p-value = 1.026972e-06

p.36r1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU36") + 
  ylim(-1, 26)+
  geom_text(x=1.5,y=22,label=expression(paste(chi^2, " = ", 23.878, ", ", "Df = ", 1, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=7.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=7.2,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p.36r1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p.36r2 <- ggplot(otu36.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(-1, 26)+
  annotate("text",x=2,y=22,label="*",size=12,color="black")+
  annotate("text",x=4,y=6,label="*",size=12,color="black")+
  annotate("text",x=5,y=5,label="*",size=12,color="black")+
  annotate("text",x=6,y=5,label="*",size=12,color="black")+
  annotate("text", x=3, y=24, label="OTU36", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=3, y=21, label="Genome Size: 5.778 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=3, y=18, label="rrn copy number: 3", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p.36r2

##OTU36.lm
otu36 <- b.otu[b.otu$OTU=="OTU36",]
otu36.r <- otu36[otu36$sampleTyp=="Root",]
otu36.r1 <- otu36.r[otu36.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, otu36.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu36.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#slope
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu36.r1$TP0
y <- otu36.r1$value
cor.test(x, y, method = "spearman")


otu36.r2 <- otu36.r[otu36.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu36.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model2 <- lm(value~TP0, otu36.r2)
lm_summary2 <- summary(lm_model2)
print(lm_summary2)
# 
slope2 <- lm_summary2$coefficients["TP0", "Estimate"]
print(slope2)

x <- otu36.r2$TP0
y <- otu36.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

library(simba)
set.seed(123)
diffslope(otu36.r1$TP0,otu36.r1$value,otu36.r2$TP0,otu36.r2$value)


p.36r <- ggplot(data=otu36.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(-1,26)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=24, label="Root_Control: R = -0.682, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=21, label="Root_N-addition: R = -0.417, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.36r


root <- p.36r1+p.36r2+p.36r
root
#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu36.root.pdf",root,height = 4.5, width = 13.5)



###OTU18###
otu18 <- b.otu[b.otu$OTU=="OTU18",]
otu18.r <- otu18[otu18$sampleTyp=="Root",]
otu18.r <- otu18.r[otu18.r$value!="0",]
root_otu1 <- otu18.r[otu18.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu18.r$value, otu18.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- otu18.r[otu18.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 0.04761749, df = 1, p-value = 0.8272621

p.18r1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU18") + 
  ylim(0, 16)+
  geom_text(x=1.5,y=15,label=expression(paste(chi^2, " = ", 0.048, ", ", "Df = ", 1, ", ", "P = ", 0.827, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=12,label="a",size=7,color="black")+
  annotate("text",x=2,y=13,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p.18r1

#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p.18r2 <- ggplot(otu18.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 16)+
  annotate("text", x=1, y=15, label="OTU18", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=13, label="Genome Size: 8.817 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=11, label="rrn copy number: 2", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p.18r2

##OTU18.lm
otu18 <- b.otu[b.otu$OTU=="OTU18",]
otu18.r <- otu18[otu18$sampleTyp=="Root",]
otu18.r1 <- otu18.r[otu18.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, otu18.r1)))

lm_model1 <- lm(value~TP0, otu18.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu18.r1$TP0
y <- otu18.r1$value
cor.test(x, y, method = "spearman")


otu18.r2 <- otu18.r[otu18.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu18.r2)))

lm_model1 <- lm(value~TP0, otu18.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu18.r2$TP0
y <- otu18.r2$value
cor.test(x, y, method = "spearman")

library(simba)
set.seed(123)
diffslope(otu18.r1$TP0,otu18.r1$value,otu18.r2$TP0,otu18.r2$value)


p.18r <- ggplot(data=otu18.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,16)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=15, label="Root_Control: R = 0.914, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=13, label="Root_N-addition: R = 0.702, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.18r


root.18 <- p.18r1+p.18r2+p.18r
root.18

#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu18.root.pdf",root.18,height = 4.5, width = 13.5)


###OTU7###
OTU7 <- b.otu[b.otu$OTU=="OTU7",]
otu7.r <- OTU7[OTU7$sampleTyp=="Root",]
otu7.r <- otu7.r[otu7.r$value!="0",]
root_otu1 <- otu7.r[otu7.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu7.r$value, otu7.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- otu7.r[otu7.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 1.028683, df = 1, p-value = 0.3104682

p_7.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU7") + 
  ylim(0, 37)+
  geom_text(x=1.5,y=33,label=expression(paste(chi^2, " = ", 1.029, ", ", "Df = ", 1, ", ", "P = ", 0.310, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=14,label="a",size=7,color="black")+
  annotate("text",x=2,y=15,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_7.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_7.2 <- ggplot(otu7.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 37)+
  annotate("text", x=1, y=34, label="OTU7", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=30, label="Genome Size: 6.467 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=1, y=26, label="rrn copy number: 6", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_7.2

##OTU7.lm
OTU7 <- b.otu[b.otu$OTU=="OTU7",]
otu7.r <- OTU7[OTU7$sampleTyp=="Root",]
otu7.r1 <- otu7.r[otu7.r$Treatment=="Control",]
shapiro.test(residuals(lm(value~TP0, otu7.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu7.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu7.r1$TP0
y <- otu7.r1$value
cor.test(x, y, method = "spearman")


otu7.r2 <- otu7.r[otu7.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu7.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu7.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu7.r2$TP0
y <- otu7.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.7r <- ggplot(data=otu7.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,37)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=34, label="Root_Control: R = 0.122, P = 0.091", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=30, label="Root_N-addition: R = -0.107, P = 0.380", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.7r


root.7 <- p_7.1+p_7.2+p.7r
root.7

#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu7.root.pdf",root.7,height = 4.5, width = 13.5)


###OTU5###
OTU5 <- b.otu[b.otu$OTU=="OTU5",]
otu5.r <- OTU5[OTU5$sampleTyp=="Root",]
otu5.r <- otu5.r[otu5.r$value!="0",]
root_otu1 <- otu5.r[otu5.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(otu5.r$value, otu5.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- otu5.r[otu5.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_5.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU5") + 
  ylim(0, 16)+
  geom_text(x=1.5,y=14,label=expression(paste(chi^2, " = ", 1.948, ", ", "Df = ", 1, ", ", "P = ", 0.163, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=10,label="a",size=7,color="black")+
  annotate("text",x=2,y=9,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_5.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_5.2 <- ggplot(otu5.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 16)+
  annotate("text", x=2.5, y=15, label="OTU5", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=13, label="Genome Size: 7.523 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=11, label="rrn copy number: 7", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text",x=2,y=9,label="*",size=12,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_5.2

##OTU7.lm
OTU5 <- b.otu[b.otu$OTU=="OTU5",]
otu5.r <- OTU5[OTU5$sampleTyp=="Root",]
otu5.r1 <- otu5.r[otu5.r$Treatment=="Control",]
otu5.r3 <- otu5.r1[otu5.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, otu5.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu5.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu5.r1$TP0
y <- otu5.r1$value
cor.test(x, y, method = "spearman")


otu5.r2 <- otu5.r[otu5.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, otu5.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, otu5.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- otu5.r2$TP0
y <- otu5.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.5r <- ggplot(data=otu5.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,16)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=15, label="Root_Control: R = 0.032, P = 0.6641", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=13, label="Root_N-addition: R = 0.429, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.5r


root.5 <- p_5.1+p_5.2+p.5r
root.5

#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu5.root.pdf",root.5,height = 4.5, width = 13.5)


###OTU2###
OTU2 <- b.otu[b.otu$OTU=="OTU2",]
OTU2.r <- OTU2[OTU2$sampleTyp=="Root",]
OTU2.r <- OTU2.r[OTU2.r$value!="0",]
root_otu1 <- OTU2.r[OTU2.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU2.r$value, OTU2.r$Treatment, p.adj = "bonferroni") ##
kw

root_otu2 <- OTU2.r[OTU2.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") ##
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_2.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU2") + 
  ylim(0, 27)+
  geom_text(x=1.5,y=24,label=expression(paste(chi^2, " = ", 14.726, ", ", "Df = ", 1, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=18,label="b",size=7,color="black")+
  annotate("text",x=2,y=19,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_2.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_2.2 <- ggplot(OTU2.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 27)+
  annotate("text",x=2,y=16,label="*",size=12,color="black")+
  annotate("text",x=7,y=16,label="*",size=12,color="black")+
  annotate("text", x=2.5, y=26, label="OTU2", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=23, label="Genome Size: 5.084 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=20, label="rrn copy number: 4", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_2.2

##OTU2.lm
OTU2 <- b.otu[b.otu$OTU=="OTU2",]
OTU2.r <- OTU2[OTU2$sampleTyp=="Root",]
OTU2.r1 <- OTU2.r[OTU2.r$Treatment=="Control",]
OTU2.r3 <- OTU2.r1[OTU2.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU2.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU2.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU2.r1$TP0
y <- OTU2.r1$value
cor.test(x, y, method = "spearman")


OTU2.r2 <- OTU2.r[OTU2.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU2.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU2.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU2.r2$TP0
y <- OTU2.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.2r <- ggplot(data=OTU2.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,27)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.2, y=26, label="Root_Control: R = -0.068, P = 0.352", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.2, y=23, label="Root_N-addition: R = -0.204, P = 0.090", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.2r


root.2 <- p_2.1+p_2.2+p.2r
root.2

#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu2.root.pdf",root.2,height = 4.5, width = 13.5)


###OTU21###
OTU21 <- b.otu[b.otu$OTU=="OTU21",]
OTU21.r <- OTU21[OTU21$sampleTyp=="Root",]
OTU21.r <- OTU21.r[OTU21.r$value!="0",]
root_otu1 <- OTU21.r[OTU21.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU21.r$value, OTU21.r$Treatment, p.adj = "bonferroni") ##
kw

root_otu2 <- OTU21.r[OTU21.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") ##
kw


p_21.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU21") + 
  ylim(0, 6.1)+
  geom_text(x=1.5,y=5.5,label=expression(paste(chi^2, " = ", 21.859, ", ", "Df = ", 1, ", ", "P < ", 0.001, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=2.1,label="a",size=7,color="black")+
  annotate("text",x=2,y=1.5,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_21.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_21.2 <- ggplot(OTU21.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 6.1)+
  annotate("text",x=2,y=3,label="*",size=12,color="black")+
  annotate("text",x=6,y=1.8,label="*",size=12,color="black")+
  annotate("text", x=2.5, y=5.5, label="OTU21", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=4.9, label="Genome Size: 5.233 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2.5, y=4.3, label="rrn copy number: 7", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_21.2

##OTU21.lm
OTU21 <- b.otu[b.otu$OTU=="OTU21",]
OTU21.r <- OTU21[OTU21$sampleTyp=="Root",]
OTU21.r <- OTU21.r[OTU21.r$value!="0",]
OTU21.r1 <- OTU21.r[OTU21.r$Treatment=="Control",]
OTU21.r3 <- OTU21.r1[OTU21.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU21.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU21.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU21.r1$TP0
y <- OTU21.r1$value
cor.test(x, y, method = "spearman")


OTU21.r2 <- OTU21.r[OTU21.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU21.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU21.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU21.r2$TP0
y <- OTU21.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.21r <- ggplot(data=OTU21.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,6.1)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=5.5, label="Root_Control: R = -0.252, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=4.9, label="Root_N-addition: R = 0.166, P = 0.204", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.21r


root.21 <- p_21.1+p_21.2+p.21r
root.21

#setwd("~/Documents/sunR/Bacteria882/n.treatment/otu")
#ggsave("B-otu21.root.pdf",root.21,height = 4.5, width = 13.5)


###OTU42###
OTU42 <- b.otu[b.otu$OTU=="OTU42",]
OTU42.r <- OTU42[OTU42$sampleTyp=="Root",]
OTU42.r <- OTU42.r[OTU42.r$value!="0",]
root_otu1 <- OTU42.r[OTU42.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU42.r$value, OTU42.r$Treatment, p.adj = "bonferroni") ##
kw

root_otu2 <- OTU42.r[OTU42.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") ##
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_42.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU42") + 
  ylim(0, 7.2)+
  geom_text(x=1.5,y=6.5,label=expression(paste(chi^2, " = ", 0.356, ", ", "Df = ", 1, ", ", "P = ", 0.551, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=1,label="a",size=7,color="black")+
  annotate("text",x=2,y=1,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_42.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_42.2 <- ggplot(OTU42.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 7.2)+
  annotate("text", x=2, y=6.5, label="OTU42", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=5.9, label="Genome Size: 2.059 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=5.3, label="rrn copy number: 8", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_42.2

##OTU42.lm
OTU42 <- b.otu[b.otu$OTU=="OTU42",]
OTU42.r <- OTU42[OTU42$sampleTyp=="Root",]
OTU42.r <- OTU42.r[OTU42.r$value!="0",]
OTU42.r1 <- OTU42.r[OTU42.r$Treatment=="Control",]
OTU42.r3 <- OTU42.r1[OTU42.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU42.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU42.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU42.r1$TP0
y <- OTU42.r1$value
cor.test(x, y, method = "spearman")


OTU42.r2 <- OTU42.r[OTU42.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU42.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU42.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
#
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU42.r2$TP0
y <- OTU42.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.42r <- ggplot(data=OTU42.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,7.2)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=6.5, label="Root_Control: R = -0.492, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=5.9, label="Root_N-addition: R = -0.222, P = 0.187", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.42r


root.42 <- p_42.1+p_42.2+p.42r
root.42

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-otu42.root.pdf",root.42,height = 4.5, width = 13.5)


###OTU47###
OTU47 <- b.otu[b.otu$OTU=="OTU47",]
OTU47.r <- OTU47[OTU47$sampleTyp=="Root",]
OTU47.r <- OTU47.r[OTU47.r$value!="0",]
root_otu1 <- OTU47.r[OTU47.r$TP00=="TP2",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU47.r$value, OTU47.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- OTU47.r[OTU47.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_47.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU47") + 
  ylim(0, 0.51)+
  geom_text(x=1.5,y=0.45,label=expression(paste(chi^2, " = ", 1.087, ", ", "Df = ", 1, ", ", "P = ", 0.297, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=0.19,label="a",size=7,color="black")+
  annotate("text",x=2,y=0.18,label="a",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_47.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_47.2 <- ggplot(OTU47.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(0, 0.51)+
  annotate("text", x=2, y=0.45, label="OTU47", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=0.39, label="Genome Size: 3.252 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=0.33, label="rrn copy number: 4", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_47.2

##OTU47.lm
OTU47 <- b.otu[b.otu$OTU=="OTU47",]
OTU47.r <- OTU47[OTU47$sampleTyp=="Root",]
OTU47.r <- OTU47.r[OTU47.r$value!="0",]
OTU47.r1 <- OTU47.r[OTU47.r$Treatment=="Control",]
OTU47.r3 <- OTU47.r1[OTU47.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU47.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU47.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU47.r1$TP0
y <- OTU47.r1$value
cor.test(x, y, method = "spearman")


OTU47.r2 <- OTU47.r[OTU47.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU47.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU47.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
# 
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU47.r2$TP0
y <- OTU47.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.47r <- ggplot(data=OTU47.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(0,0.51)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=0.45, label="Root_Control: R = -0.164, P = 0.160", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=0.39, label="Root_N-addition: R = -0.343, P = 0.15", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.47r


root.47 <- p_47.1+p_47.2+p.47r
root.47

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-OTU47.root.pdf",root.47,height = 4.5, width = 13.5)


###OTU15###
OTU15 <- b.otu[b.otu$OTU=="OTU15",]
OTU15.r <- OTU15[OTU15$sampleTyp=="Root",]
OTU15.r <- OTU15.r[OTU15.r$value!="0",]
root_otu1 <- OTU15.r[OTU15.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU15.r$value, OTU15.r$Treatment, p.adj = "bonferroni") ##
kw

root_otu2 <- OTU15.r[OTU15.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") ##
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_15.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU15") + 
  ylim(-1, 7.6)+
  geom_text(x=1.5,y=6,label=expression(paste(chi^2, " = ", 1.267, ", ", "Df = ", 1, ", ", "P = ", 0.260, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=1,label="a",size=7,color="black")+
  annotate("text",x=2,y=1,label="a",size=7,color="black")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_15.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_15.2 <- ggplot(OTU15.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(-1, 7.6)+
  annotate("text", x=2, y=6.5, label="OTU15", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=5.9, label="Genome Size: 6.121 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=5.3, label="rrn copy number: 7", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_15.2

##OTU15.lm
OTU15 <- b.otu[b.otu$OTU=="OTU15",]
OTU15.r <- OTU15[OTU15$sampleTyp=="Root",]
OTU15.r <- OTU15.r[OTU15.r$value!="0",]
OTU15.r1 <- OTU15.r[OTU15.r$Treatment=="Control",]
OTU15.r3 <- OTU15.r1[OTU15.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU15.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU15.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU15.r1$TP0
y <- OTU15.r1$value
cor.test(x, y, method = "spearman")


OTU15.r2 <- OTU15.r[OTU15.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU15.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU15.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU15.r2$TP0
y <- OTU15.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.15r <- ggplot(data=OTU15.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(-1,7.6)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=16,face="bold")) +
  annotate("text", x=0.5, y=6.5, label="Root_Control: R = -0.623, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=5.9, label="Root_N-addition: R = -0.223, P = 0.173", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.15r


root.15 <- p_15.1+p_15.2+p.15r
root.15

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-OTU15.root.pdf",root.15,height = 4.5, width = 13.5)


###OTU12###
OTU12 <- b.otu[b.otu$OTU=="OTU12",]
OTU12.r <- OTU12[OTU12$sampleTyp=="Root",]
OTU12.r <- OTU12.r[OTU12.r$value!="0",]
root_otu1 <- OTU12.r[OTU12.r$TP00=="TP7",]
kw <- kruskal(root_otu1$value, root_otu1$Treatment, p.adj = "bonferroni")
kw

kw <- kruskal(OTU12.r$value, OTU12.r$Treatment, p.adj = "bonferroni") 
kw

root_otu2 <- OTU12.r[OTU12.r$TP00 != "TP1",]

kw <- kruskal(root_otu2$value, root_otu2$Treatment, p.adj = "bonferroni") 
kw
#chi-squared = 1.948454, df = 1, p-value = 0.1627535

p_12.1 <- ggplot(root_otu2, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.7,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(fill = Treatment),alpha = 0.2, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Treatment",y = "Root - OTU12") + 
  ylim(-1, 13)+
  geom_text(x=1.5,y=10,label=expression(paste(chi^2, " = ", 9.337, ", ", "Df = ", 1, ", ", "P = ", 0.002, sep="")),
            size=6,color="black",fontface="bold")+
  annotate("text",x=1,y=2.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=2.3,label="b",size=7,color="black")+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))

p_12.1


#setwd("~/Documents/sunR/Bacteria882/results/genome")
#ggsave("B-genomesize-Ntreatment-root-boxplot.pdf",p_gs,height = 5, width = 5)


p_12.2 <- ggplot(OTU12.r, aes(x = TP00, y = value)) + 
  geom_boxplot(aes(fill = Treatment),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = Treatment),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("darkred", "black"), name="Treatment") +
  scale_fill_manual(values= c("darkred", "black"), name="Treatment") +
  labs(x = "Time",y = "") + 
  ylim(-1, 13)+
  annotate("text",x=2,y=4,label="*",size=12,color="black")+
  annotate("text",x=7,y=1.3,label="*",size=12,color="black")+
  annotate("text", x=2, y=11, label="OTU12", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=10, label="Genome Size: 4.397 Mb", 
           size=5, hjust = 0, color="black",fontface="bold") +
  annotate("text", x=2, y=9, label="rrn copy number: 1", 
           size=5, hjust = 0, color="black",fontface="bold") +
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold"))


p_12.2

##OTU12.lm
OTU12 <- b.otu[b.otu$OTU=="OTU12",]
OTU12.r <- OTU12[OTU12$sampleTyp=="Root",]
OTU12.r <- OTU12.r[OTU12.r$value!="0",]
OTU12.r1 <- OTU12.r[OTU12.r$Treatment=="Control",]
OTU12.r3 <- OTU12.r1[OTU12.r1$TP00 != "TP1",]
shapiro.test(residuals(lm(value~TP0, OTU12.r1)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU12.r1)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU12.r1$TP0
y <- OTU12.r1$value
cor.test(x, y, method = "spearman")


OTU12.r2 <- OTU12.r[OTU12.r$Treatment=="N-addition",]
shapiro.test(residuals(lm(value~TP0, OTU12.r2)))
# lm <- summary(lm(value~TP00, gs_root))
# lm
lm_model1 <- lm(value~TP0, OTU12.r2)
lm_summary1 <- summary(lm_model1)
print(lm_summary1)
slope1 <- lm_summary1$coefficients["TP0", "Estimate"]
print(slope1)

x <- OTU12.r2$TP0
y <- OTU12.r2$value
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

p.12r <- ggplot(data=OTU12.r) + 
  geom_jitter(aes(x= TP0, y = value, color = Treatment, shape = Treatment), 
              alpha = 0.4, size =2.1, height = 0) + 
  geom_smooth(aes(x= TP0, y= value, group = Treatment, color = Treatment), 
              method ="lm", linewidth = 2, show.legend = F)+
  scale_shape_manual(name="Treatment",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("darkred","black"), name="Time")+
  labs(x = "Time",y = "")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  scale_y_continuous(limits = c(-1,13)) + 
  theme_bw()+
  guides(color = guide_legend(title = "Treatment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Treatment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=18,face="bold")) +
  annotate("text", x=0.5, y=11, label="Root_Control: R = -0.649, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") +
  annotate("text", x=0.5, y=10, label="Root_N-addition: R = -0.569, P < 0.001", 
           size=5, hjust = 0, color="black",fontface="bold") 


p.12r


root.12 <- p_12.1+p_12.2+p.12r
root.12

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-OTU12.root.pdf",root.12,height = 4.5, width = 13.5)

##Supplementary.Figure6-8.titan.gs.rrn.abundance+genus####
library(TITAN2)
library(ggplot2)
library(splitstackshape)
library(ggh4x) #facet_grid axis independent
library(colorRamps)
library(RColorBrewer)
library(ggsci)
library(patchwork)
library(readxl)
library(writexl)
rm(list = ls())
##otu.relative.abundance####
setwd("~/Documents/sunR/Bacteria882")
#rm(list = ls())
##1.data
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
ID$otu <- paste(ID$names, sep = ".")

bact.otu <- aggregate(otu_all1,by=list(ID$otu) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )

bact.otu0 <- bact.relabu.otu

bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))] ## 
#lev <- interaction(group2$TP00, group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
lev <- interaction(group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.otu.lev <- aggregate(bact.otu0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

rel <- bact.otu.lev
rel1 <- setNames(as.data.frame(t(rel[,-1])), rel$Group.1)

###root####
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-root-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Root"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
#otu_genome <- otu_genome[, -1]
colnames(otu_rrna)[1] <- "sample.id"
#otu_rrna <- otu_rrna[, -1]
colnames(otu_rrna)[2] <- "rrna"

library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")


genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]
rrna_up <- mean(genome_up$rrna)
genome_up0 <- na.omit(genome_up)
gs_up <- mean(genome_up0$genome_size/1000000)
gc_up <- mean(genome_up0$gc_percentage)
genome_down <- genome[genome$Titan == "z-", ]
rrna_down <- mean(genome_down$rrna)
genome_down0 <- na.omit(genome_down)
gs_down <- mean(genome_down0$genome_size/1000000)
gc_down <- mean(genome_down0$gc_percentage)



color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)


rel1_filtered <- rel1[rownames(rel1) %in% dat3$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)

dat4 <- merge(dat3, rel1_filtered, by="sample.id")



p1 <- ggplot()+
  geom_linerange(data=dat4,
                 aes(x=order.x, y=zenv.cp, ymin=dat4$`5%`,ymax = dat4$`95%`,
                     color = factor(phylum.x), shape=Titan,linetype=Titan),
                 alpha = 0.8)+
  geom_point(data = dat4,
             aes(x = order.x, y = zenv.cp, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  #scale_size_continuous(range = c(.1, 2)) + 
  scale_color_manual(values=color) +
  scale_shape_manual(values=c(16, 1)) + 
  coord_flip() + 
  ylab("Time") + xlab("Root - OTUs") + 
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(size = guide_legend(title= "Abundance"), 
         color = guide_legend(title= "Phylum",order=2,override.aes = list(size = 3)),
         shape = guide_legend(title= "Titan",order=1, nrow = 1,override.aes = list(size = 3)),
         linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=18, face="bold"),
    legend.text = element_text(colour="black", size=15, face="bold"),
    #legend.position = "none",
    axis.text = element_text(size = 12,face="bold",color="black"),
    axis.title=element_text(size=18,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"))

p1


p2 <-   ggplot(dat4, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$genome_size)/1000000, 
               yend = mean(genome_down$genome_size)/1000000,
               x = 0, xend = 55, colour = "red") +
  geom_segment(y = mean(genome_up$genome_size)/1000000, 
               yend = mean(genome_up$genome_size)/1000000,
               x = 56, xend = 130, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = genome_size/1000000, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=18,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p4 <- ggplot(dat4, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$rrna), yend = mean(genome_down$rrna),
               x = 0, xend = 55, colour = "red") +
  geom_segment(y = mean(genome_up$rrna), yend = mean(genome_up$rrna),
               x = 56, xend = 130, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = rrna, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  geom_text(data = dat4,
            aes(x = order.x + 0.5, y = 0, label = genus.y, color = phylum.x), 
            size = 1, hjust = 0) +
  theme_bw() + 
  coord_flip() + 
  scale_shape_manual(values = c(16, 1)) + 
  ylab("RRN") + 
  xlab("") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.line.x = element_line(linetype = "solid"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(color = "white"),
    legend.position = "none"
  ) +
  scale_color_manual(values = color)

p4

p <- p1+p2+p4+ plot_layout(guides = 'collect')
p
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-titan.gs.rrn.abundance.root.pdf",p,height = 6, width = 11)


###leaf####
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-leaf-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Leaf"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
#otu_genome <- otu_genome[, -1]
colnames(otu_rrna)[1] <- "sample.id"
#otu_rrna <- otu_rrna[, -1]
colnames(otu_rrna)[2] <- "rrna"


library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")


genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]
rrna_up <- mean(genome_up$rrna)
genome_up0 <- na.omit(genome_up)
gs_up <- mean(genome_up0$genome_size/1000000)
gc_up <- mean(genome_up0$gc_percentage)
genome_down <- genome[genome$Titan == "z-", ]
rrna_down <- mean(genome_down$rrna)
genome_down0 <- na.omit(genome_down)
gs_down <- mean(genome_down0$genome_size/1000000)
gc_down <- mean(genome_down0$gc_percentage)


color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)

# 
rel1_filtered <- rel1[rownames(rel1) %in% dat3$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)

dat4 <- merge(dat3, rel1_filtered, by="sample.id")


p1 <- ggplot()+
  geom_linerange(data=dat4,
                 aes(x=order.x, y=zenv.cp, ymin=dat4$`5%`,ymax = dat4$`95%`,
                     color = factor(phylum.x), shape=Titan,linetype=Titan),
                 alpha = 0.8)+
  geom_point(data = dat4,
             aes(x = order.x, y = zenv.cp, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  #scale_size_continuous(range = c(.1, 2)) + 
  scale_color_manual(values=color) +
  scale_shape_manual(values=c(16, 1)) + 
  coord_flip() + 
  ylab("Time") + xlab("Leaf - OTUs") + 
  #scale_size_continuous(range=c(1,2))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(size = guide_legend(title= "Abundance"), 
         color = guide_legend(title= "Phylum",order=2,override.aes = list(size = 3)),
         shape = guide_legend(title= "Titan",order=1, nrow = 1,override.aes = list(size = 3)),
         linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  #theme_set(theme_minimal() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))) +
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=18, face="bold"),
    legend.text = element_text(colour="black", size=15, face="bold"),
    #legend.position = "none",
    axis.text = element_text(size = 12,face="bold",color="black"),
    axis.title=element_text(size=18,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"))

p1


p2 <-   ggplot(dat4, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$genome_size)/1000000, 
               yend = mean(genome_down$genome_size)/1000000,
               x = 0, xend = 29, colour = "red") +
  geom_segment(y = mean(genome_up$genome_size)/1000000, 
               yend = mean(genome_up$genome_size)/1000000,
               x = 30, xend = 47, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = genome_size/1000000, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=18,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        #axis.line.x.top = element_line(color = "white"),
        #axis.line.y.right = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p4 <- ggplot(dat4, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$rrna), yend = mean(genome_down$rrna),
               x = 0, xend = 29, colour = "red") +
  geom_segment(y = mean(genome_up$rrna), yend = mean(genome_up$rrna),
               x = 30, xend = 47, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = rrna, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  geom_text(data = dat4,
            aes(x = order.x + 0.5, y = 0, label = genus.y, color = phylum.x), 
            size = 2, hjust = 0) +
  theme_bw() + 
  coord_flip() + 
  scale_shape_manual(values = c(16, 1)) + 
  ylab("RRN") + 
  xlab("") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.line.x = element_line(linetype = "solid"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(color = "white"),
    legend.position = "none"
  ) +
  scale_color_manual(values = color)

p4

p <- p1+p2+p4+ plot_layout(guides = 'collect')
p
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-titan.gs.rrn.abundance.leaf.pdf",p,height = 6, width = 11)




###air####
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-air-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Air"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
colnames(otu_rrna)[1] <- "sample.id"
colnames(otu_rrna)[2] <- "rrna"


library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")

genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]

gs_up <- mean(genome_up$genome_size/1000000, na.rm=T)
gc_up <- mean(genome_up0$gc_percentage, na.rm=T)
rrna_up <- mean(genome_up$rrna, na.rm=T)

genome_down <- genome[genome$Titan == "z-", ]

gs_down <- mean(genome_down0$genome_size/1000000, na.rm=T)
gc_down <- mean(genome_down0$gc_percentage, na.rm=T)
rrna_down <- mean(genome_down$rrna, na.rm=T)



color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)


rel1_filtered <- rel1[rownames(rel1) %in% dat3$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)

dat4 <- merge(dat3, rel1_filtered, by="sample.id")


p1 <- ggplot()+
  geom_linerange(data=dat4,
                 aes(x=order.x, y=zenv.cp, ymin=dat4$`5%`,ymax = dat4$`95%`,
                     color = factor(phylum.x), shape=Titan,linetype=Titan),
                 alpha = 0.8)+
  geom_point(data = dat4,
             aes(x = order.x, y = zenv.cp, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  #scale_size_continuous(range = c(.1, 2)) + 
  scale_color_manual(values=color) +
  scale_shape_manual(values=c(16, 1)) + 
  coord_flip() + 
  ylab("Time") + xlab("Air - OTUs") + 
  #scale_size_continuous(range=c(1,2))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(size = guide_legend(title= "Abundance", nrow=1), 
         color = guide_legend(title= "Phylum",order=2,override.aes = list(size = 3)),
         shape = guide_legend(title= "Titan",order=1, nrow = 1,override.aes = list(size = 3)),
         linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  #theme_set(theme_minimal() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))) +
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=14, face="bold"),
    legend.text = element_text(colour="black", size=12, face="bold"),
    #legend.position = "none",
    axis.text = element_text(size = 12,face="bold",color="black"),
    axis.title=element_text(size=18,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"))

p1


p2 <-   ggplot(dat4, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = gs_down, 
               yend = gs_down,
               x = 0, xend = 198, colour = "red") +
  geom_segment(y = gs_up, 
               yend = gs_up,
               x = 199, xend = 325, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = genome_size/1000000, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) + 
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=18,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p4 <- ggplot(dat4, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = rrna_down, yend = rrna_down,
               x = 0, xend = 198, colour = "red") +
  geom_segment(y = rrna_up, yend = rrna_up,
               x = 199, xend = 325, colour = "blue") +
  geom_point(data = dat4,
             aes(x = order.x, y = rrna, shape = Titan,
                 color = factor(phylum.x), 
                 size = dat4$`Root:Control`),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(1, 5)) +  
  geom_text(data = dat4,
            aes(x = ifelse(order.x %% 2 == 1,order.x - .5,order.x + 0.5),
                y = 10,
                label = genus.y,
                color = phylum.x,
                hjust = ifelse(order.x %% 2 == 1, 1, 0)),size = 1) +
  theme_bw() + 
  coord_flip() + 
  scale_shape_manual(values = c(16, 1)) + 
  ylab("RRN") + 
  xlab("") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 18, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.line.x = element_line(linetype = "solid"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(color = "white"),
    legend.position = "none"
  ) +
  scale_color_manual(values = color)

p4

p <- p1+p2+p4+ plot_layout(guides = 'collect')
p
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-titan.gs.rrn.abundance.air.pdf",p,height = 6, width = 12)




##Supplementary.Figure9.Function####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)

rm(list=ls())

setwd("~/Documents/sunR/Bacteria882/new-picrust2")
data <- read_xlsx("function_proportion_n.xlsx", sheet=1)
data_subdata <- read_xlsx("data_subdata_n.xlsx",sheet=3)

p_d <- ggplot(data=data) + 
  geom_jitter(aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =1.5) + 
  geom_text(data = data_subdata,
            aes(x=2.0,y=anno_y1,label = rlt_sig,color=sampleType),
            parse = F, show.legend = F, size = 5,hjust = 0,
            fontface = "bold")+
  facet_wrap(~level2, scales = "free") +
  geom_smooth(aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Mean proportion (%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(colour = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=5),nrow=1),
         shape = guide_legend(title = "Compartment",nrow=1))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold",color="black"),
        legend.title = element_text(colour="black", size=17, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        legend.position = "inside",
        legend.position.inside = c(0.8,0.1),
        strip.text = element_text(size=11,face="bold"),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=20,face="bold",color="black")) 

p_d

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("S.Figure.function.pdf", p_d, height = 18, width = 24)




##Supplementary.Figure10a.Function.heatmap####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
library(readxl)
library(dplyr)
##air
rm(list=ls())
##
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group1 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- group1 %>% filter(Treatment=="Control")

group3 <- group2 %>% filter(sampleType == "Air")

#kegg.data
kegg <- read.table("kegg htext.txt", sep = "\t",fill = TRUE,header = T,quote = "")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
# KO.data
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)

# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group3$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function")  #
ko_abundance <-ko_abundance[rowSums(abundance) != 0,]  # 

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  #

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group3$sample.id)]  # 
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  # 

# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  # 
table(duplicated(ko_abundance5$pathway_id))  # 

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]  # 

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:359],FUN="sum")  # 
ko_abundance6[, 2:358]  <- apply(ko_abundance6[, 2:358], 2, function(x) x / sum(x))  # 


# 2. heatmap
ko_air <- column_to_rownames(ko_abundance6,"level2")
ko_air <- as.data.frame(t(ko_air))
group3 <- group3[,c(1,4)]
group3 <- column_to_rownames(group3, "sample.id")
group3 <- group3[rownames(ko_air),]

corr_matrix <- corr.test(ko_air, group3, method = 'spearman', adjust = "fdr") 
corr_matrix$r  # r matrix
corr_matrix$p <- corr_matrix$p.adj  # p.adj matrix
corr_matrix$p
corr_matrix$p.adj

p2_1 <- qcorrplot(corr_matrix) +
  geom_square() +
  geom_mark(sep = '\n',sig_thres = 0.05, size = 2, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limit = c(-0.55, 0.55)) +
  guides(fill = guide_colorbar(title = "Spearman's Rho")) +
  theme(plot.margin = margin(0, 0, 0, 0.5, "cm"),
        plot.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"), 
        legend.background = element_rect(fill = "transparent", color = NA))

print(p2_1)

#ggsave("B-picrust-heatmap-air.pdf", p2_1, height = 12, width = 12)

##root
group4 <- group2 %>% filter(sampleType == "Root")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
# 
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)
# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group4$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function")  # 
ko_abundance <-ko_abundance[rowSums(abundance) != 0,]  # 

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  # 

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group4$sample.id)]  # 
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  # 

# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  # 
table(duplicated(ko_abundance5$pathway_id))  # 

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]  # 

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:194],FUN="sum")  # 

ko_abundance6[, 2:193]  <- apply(ko_abundance6[, 2:193], 2, function(x) x / sum(x))  #

# 2. heatmap
ko_root <- column_to_rownames(ko_abundance6,"level2")
ko_root <- as.data.frame(t(ko_root))

group4 <- group4[,c(1,4)]
group4 <- column_to_rownames(group4, "sample.id")

group4 <- group4[rownames(ko_root),]

corr_matrix1 <- corr.test(ko_root, group4, method = 'spearman', adjust = "fdr") 
corr_matrix1$r  # r matrix
corr_matrix1$p <- corr_matrix1$p.adj  # p.adj matrix
corr_matrix1$p
corr_matrix1$p.adj


p2_2 <- qcorrplot(corr_matrix1) +
  geom_square() +
  geom_mark(sep = '\n',sig_thres = 0.05, size = 2, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limit = c(-0.55, 0.55)) +
  #scale_x_discrete(labels = c(Col1 = "Root-Time"), limit = c()) +
  guides(fill = guide_colorbar(title = "Spearman's Rho")) +
  theme(plot.margin = margin(0, 0, 0, 0.5, "cm"),
        plot.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),  
        #legend.position = c(0.8,0.15),  # c(x,y)
        legend.background = element_rect(fill = "transparent", color = NA))

print(p2_2)

#ggsave("B-picrust-heatmap-root.pdf", p2_2, height = 12, width = 12)




##leaf
group5 <- group2 %>% filter(sampleType == "Leaf")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
# 
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)
# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group5$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function")  
ko_abundance <-ko_abundance[rowSums(abundance) != 0,] 

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  #

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group5$sample.id)]  # 
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  # 

# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  
table(duplicated(ko_abundance5$pathway_id)) 

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),] 


ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:154],FUN="sum") 

ko_abundance6[, 2:153]  <- apply(ko_abundance6[, 2:153], 2, function(x) x / sum(x))  

# 2. heatmap
ko_leaf <- column_to_rownames(ko_abundance6,"level2")
ko_leaf <- as.data.frame(t(ko_leaf))

group5 <- group5[,c(1,4)]
group5 <- column_to_rownames(group5, "sample.id")

group5 <- group5[rownames(ko_leaf),]

corr_matrix2 <- corr.test(ko_leaf, group5, method = 'spearman', adjust = "fdr") 
corr_matrix2$r  # r matrix
#corr_matrix2$p <- corr_matrix2$p.adj  # p.adj matrix
corr_matrix2$p
corr_matrix2$p.adj

p2_3 <- qcorrplot(corr_matrix2) +
  geom_square() + 
  geom_mark(sep = '\n',sig_thres = 0.05, size = 2, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                       limit = c(-0.55, 0.55)) +
  guides(fill = guide_colorbar(title = "Spearman's Rho")) +
  theme(plot.margin = margin(0, 0, 0, 0.5, "cm"),
        plot.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),  
        #legend.position = c(0.8,0.15),  # c(x,y)
        legend.background = element_rect(fill = "transparent", color = NA))

print(p2_3)

#ggsave("B-picrust-heatmap-leaf.pdf", p2_3, height = 12, width = 12)


air <- data.frame("Air-Time"=corr_matrix$r, p1=corr_matrix$p.adj)
air$functions <- rownames(air)
leaf <- data.frame("Leaf-Time"=corr_matrix2$r, p2=corr_matrix2$p.adj)
leaf$functions <- rownames(leaf)
root <- data.frame("Root-Time"=corr_matrix1$r, p3=corr_matrix1$p.adj)
root$functions <- rownames(root)

df <- merge(air, leaf, by = "functions")
df <- merge(df, root, by = "functions")
df0 <- df
rownames(df) <- df$functions
df <- df[,-1]
df1 <- df[, c("Air.Time", "Leaf.Time", "Root.Time")]
df2 <- df[, c("p1", "p2", "p3")]
data <- as.matrix(df1)
pvalue <- as.matrix(df2)

#
add_stars <- function(pvalue) {
  ifelse(pvalue < 0.001, "***",
         ifelse(pvalue < 0.01, "**",
                ifelse(pvalue < 0.05, "*", "")))
}

# 
stars_matrix <- matrix("", nrow=nrow(pvalue), ncol=ncol(pvalue))
for (i in 1:nrow(pvalue)) {
  for (j in 1:ncol(pvalue)) {
    stars_matrix[i, j] <- add_stars(pvalue[i, j])
  }
}


p2 <- qcorrplot(data, p.mat = pvalue) + 
  geom_square() + 
  geom_mark(sep = '\n',sig_thres = 0.05, size = 2, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  guides(fill = guide_colorbar(title = "Spearman's Rho")) +
  theme(plot.margin = margin(0, 0, 0, 0.5, "cm"),
        plot.title = element_blank(),
        axis.text = element_text(size = 12,face = "bold"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),  
        #legend.position = c(0.8,0.15),  # c(x,y)
        legend.background = element_rect(fill = "transparent", color = NA))

print(p2)
#ggsave("B-heatmap-level2.pdf", p2, height = 12, width = 12)


library(pheatmap)
pheatmap(data,
         scale = "none", 
         cluster_row = T, 
         cluster_col = F, 
         color = RColorBrewer::brewer.pal(11, "RdBu"),
         display_numbers = stars_matrix, 
         fontsize_row = 10,      
         fontsize_col = 10,     
         border_color = NA,  
         fontsize_number = 10, 
         number_color = "white",
         cellwidth = 20, 
         cellheight = 20)

#ggsave("B-pheatmap-level2.pdf", p, height = 9, width = 9)


# 
library(ggplot2)
library(reshape2)

# 
p_value <- df0[, c("functions", "p1", "p2", "p3")]
names(p_value) <- c("functions", "Air.Time", "Leaf.Time", "Root.Time")
r_value <- df0[, c("functions", "Air.Time", "Leaf.Time", "Root.Time")]
data_melt1 <- melt(r_value, id.vars = "functions")
names(data_melt1)<-c("functions","Type", "rvalue")
data_melt2 <- melt(p_value, id.vars = "functions")
names(data_melt2)<-c("functions","Type", "pvalue")

data_melt <- merge(data_melt1, data_melt2, by = c("functions", "Type"))

# 
star_function <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# 
data_melt$P_star <- sapply(data_melt$pvalue, star_function)

data_melt0 <- data.frame(data_melt)

# 
p <- ggplot(data_melt, aes(x = Type, y = functions, fill = rvalue)) +
  geom_tile() +
  coord_equal() + 
  geom_text(aes(label = paste0("", format(round(rvalue, 3), 3), "\n", P_star)), 
            vjust = 0.6, color = "white", size = 2) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  theme_minimal() +
  ylab("") + xlab("") +
  theme(panel.grid = element_blank(),
        legend.position=c(1.5, 0.5),
        axis.title=element_text(size=12,face="bold",color="black"),
        axis.text = element_text(size = 11, face = "bold",color="black"),
        axis.text.x = element_text(size = 12, face = "bold",color="black",
                                   angle = -90, vjust = 0.5, hjust = 0))

p

#ggsave("B-heatmap-ggplot.pdf", p, height = 10, width = 8)

library(ggtree)

df1.clust<-hclust(dist(df1))
p1.1 <- ggtree(df1.clust) 

p1.1
p_tree <- ggtree(df1.clust) + 
  geom_tiplab() +  
  theme(legend.position = "none") 
p_tree

library(aplot)

p_with_tree <- insert_left(p, p_tree, width = 0)
p_with_tree

setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-function.heatmap.pdf", p_with_tree, height = 10, width = 6)


##Supplementary.Figure10bcd.function.vocano.level2.n####
###root, leaf####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)  # https://github.com/kwstat/pals/issues/3
library(RColorBrewer)
library(readxl)
##1. 
rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
# 
kegg <- read.table("kegg htext.txt", sep = "\t",fill = TRUE,header = T,quote = "")
# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group1 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group1 <- group1 %>% filter(Treatment=="Control")
group2 <- group1 %>% filter(sampleType == "Root" | sampleType == "Leaf")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
#
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)

# 
ko_abundance <- ko_abundance[, colnames(ko_abundance) %in% c("function",group2$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function") 
ko_abundance <-ko_abundance[rowSums(abundance) != 0,]  

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  # 

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group2$sample.id)] 
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  

# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  
table(duplicated(ko_abundance5$pathway_id)) 

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]  

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

#
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:346],FUN="sum")  

##2.DEseq2
library(DESeq2)
library(pheatmap)

rawdata <- ko_abundance6  
rownames(rawdata) <- rawdata$level2
diffcount <- rawdata[,2:345]

diffcount <- diffcount[, group2$sample.id] 

head(diffcount)

dds <- DESeqDataSetFromMatrix(countData = round(diffcount), colData = group2, design = ~sampleType)

dds <- DESeq(dds)

#
res <- results(dds,contrast = c("sampleType","Root","Leaf"), 
               pAdjustMethod = 'fdr', alpha = 0.05)  

summary(res)
plotMA(res, alpha = 0.05, ylim = c(-3, 3))

#
diff_res <- as.data.frame(res)
diff_res$pathway_id <- rownames(diff_res)
diff_res<- na.omit(diff_res)

diff_res[which(diff_res$padj %in% NA),'sig'] <- 'no diff'
diff_res[which(diff_res$log2FoldChange >= 0.1 & diff_res$padj < 0.05),'sig'] <- 'rich (p.adj < 0.05, log2FC >= 0.1)'
diff_res[which(diff_res$log2FoldChange <= -0.1 & diff_res$padj < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -0.1)'
diff_res[which(abs(diff_res$log2FoldChange) < 0.1 | diff_res$padj >= 0.05),'sig'] <- 'no diff'


#setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#write.csv(diff_res,"KO.level2(root-leaf).result.up.down.csv",row.names = FALSE)


library(ggplot2)
library(ggrepel)
library(stringr)
#dataset <- read.csv("~/Documents/sunR/Bacteria882/new-picrust2/KO.Pathway3(root-leaf).result.up.down.csv", head = T)

p_rl <- ggplot(diff_res, aes(x = -log2FoldChange, y = -log10(pvalue), colour=sig)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_manual(values=c("red", "grey", "blue"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.02),lty=4,col="grey",lwd=0.5) +
  xlim(-1.6,1.6)+
  ylim (0,200)+
  scale_y_continuous(breaks=seq(0,200,by=40))+
  geom_text_repel(data=subset(diff_res, abs(log2FoldChange) >= 0.1 & padj < 0.05),
                  vjust="inward",hjust="inward", 
                  aes(label=pathway_id, color = sig), angle = 0, size=3)+
  labs(x="log2 (fold change)", y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_tex3t(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold",color="black")) +
  ggtitle("Root enriched                                         Leaf enriched")

p_rl



###root, air####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)  # https://github.com/kwstat/pals/issues/3
library(RColorBrewer)
library(readxl)
##1.
#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
# 
kegg <- read.table("kegg htext.txt", sep = "\t",fill = TRUE,header = T,quote = "")
# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group1 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group1 <- group1 %>% filter(Treatment=="Control")
group2 <- group1 %>% filter(sampleType == "Root" | sampleType == "Air")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
#
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)

# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group2$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function") 
ko_abundance <-ko_abundance[rowSums(abundance) != 0,] 

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO))) 

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group2$sample.id)]  
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum) 
# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  
table(duplicated(ko_abundance5$pathway_id)) 

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]  

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:551],FUN="sum")  

##2. DEseq2
library(DESeq2)
library(pheatmap)

rawdata <- ko_abundance6 
rownames(rawdata) <- rawdata$level2
diffcount <- rawdata[,2:550]

diffcount <- diffcount[, group2$sample.id] 

head(diffcount)

# 
dds <- DESeqDataSetFromMatrix(countData = round(diffcount), colData = group2, design = ~sampleType)

# 
dds <- DESeq(dds)

# 
res <- results(dds,contrast = c("sampleType","Root","Air"), 
               pAdjustMethod = 'fdr', alpha = 0.05)  

summary(res)
plotMA(res, alpha = 0.05, ylim = c(-3, 3))


diff_res <- as.data.frame(res)
diff_res$pathway_id <- rownames(diff_res)
diff_res<- na.omit(diff_res)

#
diff_res[which(diff_res$padj %in% NA),'sig'] <- 'no diff'
diff_res[which(diff_res$log2FoldChange >= 0.1 & diff_res$padj < 0.05),'sig'] <- 'rich (p.adj < 0.05, log2FC >= 0.1)'
diff_res[which(diff_res$log2FoldChange <= -0.1 & diff_res$padj < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -0.1)'
diff_res[which(abs(diff_res$log2FoldChange) < 0.1 | diff_res$padj >= 0.05),'sig'] <- 'no diff'


#setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#write.csv(diff_res,"KO.level2(root-air).result.up.down.csv",row.names = FALSE)

library(ggplot2)
library(ggrepel)
library(stringr)
#dataset <- read.csv("~/Documents/sunR/Bacteria882/new-picrust2/KO.Pathway3(root-leaf).result.up.down.csv", head = T)

p_ra <- ggplot(diff_res, aes(x = -log2FoldChange, y = -log10(pvalue), colour=sig)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_manual(values=c("red", "grey", "blue"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.02),lty=4,col="grey",lwd=0.5) +
  xlim(-1.6,1.6)+
  ylim (0,200)+
  scale_y_continuous(breaks=seq(0,200,by=40))+
  geom_text_repel(data=subset(diff_res, abs(log2FoldChange) >= 0.1 & padj < 0.05),
                  vjust="inward",hjust="inward", 
                  aes(label=pathway_id, color = sig), angle = 0, size=3)+
  labs(x="log2 (fold change)", y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold",color="black")) +
  ggtitle("Root enriched                                            Air enriched")

p_ra





###leaf, air####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)  # https://github.com/kwstat/pals/issues/3
library(RColorBrewer)
library(readxl)
##1. 
#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#
kegg <- read.table("kegg htext.txt", sep = "\t",fill = TRUE,header = T,quote = "")

# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group1 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group1 <- group1 %>% filter(Treatment=="Control")
group2 <- group1 %>% filter(sampleType == "Leaf" | sampleType == "Air")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
# 
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)

# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group2$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function")  # 
ko_abundance <-ko_abundance[rowSums(abundance) != 0,]  

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group2$sample.id)]  
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  

#
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  
table(duplicated(ko_abundance5$pathway_id))  

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]  

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:511],FUN="sum") 

##2. DEseq2
library(DESeq2)
library(pheatmap)

rawdata <- ko_abundance6  
rownames(rawdata) <- rawdata$level2
diffcount <- rawdata[,2:510]

diffcount <- diffcount[, group2$sample.id] 

head(diffcount)

#
dds <- DESeqDataSetFromMatrix(countData = round(diffcount), colData = group2, design = ~sampleType)

#
dds <- DESeq(dds)

#
res <- results(dds,contrast = c("sampleType","Leaf","Air"), 
               pAdjustMethod = 'fdr', alpha = 0.05)  

summary(res)
plotMA(res, alpha = 0.05, ylim = c(-3, 3))


diff_res <- as.data.frame(res)
diff_res$pathway_id <- rownames(diff_res)
diff_res<- na.omit(diff_res)

diff_res[which(diff_res$padj %in% NA),'sig'] <- 'no diff'
diff_res[which(diff_res$log2FoldChange >= 0.1 & diff_res$padj < 0.05),'sig'] <- 'rich (p.adj < 0.05, log2FC >= 0.1)'
diff_res[which(diff_res$log2FoldChange <= -0.1 & diff_res$padj < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -0.1)'
diff_res[which(abs(diff_res$log2FoldChange) < 0.1 | diff_res$padj >= 0.05),'sig'] <- 'no diff'

#setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#write.csv(diff_res,"KO.level2(leaf-air).result.up.down.csv",row.names = FALSE)


library(ggplot2)
library(ggrepel)
library(stringr)
#dataset <- read.csv("~/Documents/sunR/Bacteria882/new-picrust2/KO.Pathway3(root-leaf).result.up.down.csv", head = T)

p_la <- ggplot(diff_res, aes(x = -log2FoldChange, y = -log10(pvalue), colour=sig)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_manual(values=c("red", "grey", "blue"))+
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="grey",lwd=0.5) +
  geom_hline(yintercept = -log10(0.02),lty=4,col="grey",lwd=0.5) +
  xlim(-1.6,1.6)+
  ylim (0,200)+
  scale_y_continuous(breaks=seq(0,200,by=20))+
  geom_text_repel(data=subset(diff_res, abs(log2FoldChange) >= 0.1 & padj < 0.05),
                  vjust="inward",hjust="inward", 
                  aes(label=pathway_id, color = sig), angle = 0, size=3)+
  labs(x="log2 (fold change)", y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", 
        legend.title = element_blank(),
        axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=15,face="bold",color="black")) +
  ggtitle("Leaf enriched                                            Air enriched")

p_la


p <- p_rl+p_ra+p_la

p

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-function.vocano.level2.n.rl.pdf",p_rl,height = 5, width = 5.5)
#ggsave("B-function.vocano.level2.n.ra.pdf",p_ra,height = 5, width = 5.5)
#ggsave("B-function.vocano.level2.n.la.pdf",p_la,height = 5, width = 5.5)

##Supplementary.Figure11.Ncycle.revised####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

group2 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group2)
group3 <- group2 %>% filter(Treatment == "Control")


ko <- read_xlsx("ko_gene.xlsx", sheet = 1)
ko_gene <- read_xlsx("ko_gene.xlsx", sheet=2)
ko_ncycle <- read_xlsx("ko_gene.xlsx",sheet=3)

ko1 <- ko[,colnames(ko) %in% c("KO",group3$sample.id)]
ko1[, 2:702]  <- apply(ko1[, 2:702], 2, function(x) x / sum(x)) 

ko2 <- ko1[ko1$KO %in% ko_ncycle$KO,]
data1 <- merge(ko2,ko_ncycle,by="KO")
rownames(data1) <- data1$gene
data1 <- data1[,!(names(data1) %in% c("KO","gene","KO_name","EC"))]

data10 <- data.frame(t(data1))
data10$sample.id <- rownames(data10)
data11 <- merge(data10,group2,by="sample.id")

data12 <- gather(data11,key = "Ngene",value = "value", -"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")

lm_rlt_Spcs <- function(df, yval, xval,ypos, subgrp = c("all")) {
  ypos_val <- ypos
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Ngene == subgrp)# change
    
  } else {
    
    df_tmp <- df
    
  }
  
  lmR_spearman <- cor.test(df_tmp[[yval]], df_tmp[[xval]], method = "spearman")
  lmR_spe_rho <- lmR_spearman$estimate %>% round(3)
  lmR_spe_pval <- lmR_spearman$p.value
  
  if(is.na(lmR_spe_pval)) {
    
    Rp_sig <- "NA"
    
  } else if(lmR_spe_pval < 0.001) {
    
    Rp_sig <- "***"
    
  } else if(lmR_spe_pval <= 0.01) {
    
    Rp_sig <- "**"
    
  } else if(lmR_spe_pval <= 0.05){
    
    Rp_sig <- "*"
    
  } else if(lmR_spe_pval <= 1){
    
    Rp_sig <- "NS"
    
  }
  
  data_sub <- data12 %>% filter(Ngene == subgrp)
  
  
  df_rlt <- data.frame(
    Ngene = subgrp,
    anno_x1 = (range(data_sub[[xval]])[2] - range(data_sub[[xval]])[1]) * 0.5 + range(data_sub[[xval]])[1],
    anno_y1 = (range(data_sub[[yval]], na.rm = T)[2] - range(data_sub[[yval]], na.rm = T)[1]) * ypos_val +
      range(data_sub[[yval]], na.rm = T)[1],
    rsig = str_c("R = ", lmR_spe_rho),
    psig = Rp_sig,
    anno_x2 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.6 + range(df_tmp[[xval]])[1],
    anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
      range(df_tmp[[yval]], na.rm = T)[1],
    pval = lmR_spe_pval
  )
  
  pval_sig <- str_c("P = ", lmR_spe_pval %>% round(3))
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, ", ", pval_sig))
  
  return(df_rlt_F)
  
}

data_air <- data12[data12$sampleType=="Air",]
data_leaf <- data12[data12$sampleType=="Leaf",]
data_root <- data12[data12$sampleType=="Root",]

level2_list <- data12$Ngene %>% unique()

data_air_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_air,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 1.1,
                                     subgrp = .))

data_air_subdata$sampleType <- "Air"

data_leaf_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_leaf,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos=1.0,
                                     subgrp = .))

data_leaf_subdata$sampleType <- "Leaf"

data_root_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_root,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 0.9,
                                     subgrp = .))

data_root_subdata$sampleType <- "Root"

data_subdata1 <- rbind(data_air_subdata,data_leaf_subdata,data_root_subdata)
#write.csv(data_subdata1, "nitrogen.gene.r.p_n.csv")
data_subdata <- read_xlsx("nitrogen.gene.r.p_n.xlsx",sheet=1)




p_d <- ggplot(data=data12) + 
  geom_jitter(aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =1.5) + 
  geom_text(data = data_subdata,
            aes(x=1.5,y=anno_y1,label = rlt_sig,color=sampleType),
            parse = F, show.legend = F, size = 5,hjust = 0,fontface = "bold")+
  facet_wrap(~Ngene, scales = "free") +
  geom_smooth(aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Mean proportion (%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))+
  theme_bw()+
  guides(colour = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=5),nrow=1),
         shape = guide_legend(title = "Compartment",nrow=1))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold",color="black"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=17, face="bold"),
        #legend.position = "none", 
        legend.position = "inside",
        legend.position.inside = c(0.9,0.08),
        #legend.box = "vertical",
        strip.text = element_text(size=18,face="bold"),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=25,face="bold",color="black")) 

p_d
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("S.Figure8.ncycle.pdf", p_d, height = 20, width = 22)






##Supplementary.Figure12.GS, RRN(indentity>95%&qcovs>95%)-OTU-relative abundance-ratio####
library(ggplot2)
library(dplyr)
library(tidyr)

# 1.
setwd("~/Documents/sunR/Bacteria882")
otu_table <- otu_table <- read.csv("otu_bactFlattening825.csv")
filtered_otu_data <- otu_table[rowSums(otu_table[, -1] > 0) > 0, ]

group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),] 
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)

otu_rowname <- otu_table$OTU_id
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##2.genus
bact.phy <- aggregate(otu_table[, 2:826],by=list(ID$names) , sum) 
rownames(bact.phy) <- bact.phy[,1]; 
bact.phy <- bact.phy[,-1] 
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

bact.phy0 <- bact.phy0[,order(-colSums(bact.phy0))] ## 
lev <- interaction(group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.phy.lev <- aggregate(bact.phy0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

total_abundance <- rowSums(bact.phy.lev[,-1])
###Genomesize-OTUs####
setwd("~/Documents/sunR/Bacteria882/GTDB")

group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")

blast <- read_xlsx("blast.uotus.bacteria.gtdb.xlsx", sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)


# GTDB-metadata
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

blast1 <- blast0[blast0$qseqid %in% otu_table$OTU_id,]

blast1 <- blast1 %>% separate(sseqid, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$qseqid),] 

blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")

fil_otu <- blast2$qseqid
keep_cols <- unique(c("Group.1", fil_otu))
selected_cols <- bact.phy.lev %>% dplyr::select(keep_cols)

selected_abundance <- rowSums(selected_cols[,-1])

df <- data.frame(
  Group.1 = c("Air:Control", "Leaf:Control", "Root:Control", "Leaf:N-addition", "Root:N-addition"),
  Proportion = c(0.7332451, 0.9306579, 0.9414089, 0.8772315, 0.9211929)
)
df_long <- df %>%
  mutate(Other = 1 - Proportion) %>%
  pivot_longer(cols = c(Proportion, Other), names_to = "Type", values_to = "Value")


library(ggplot2)

p1 <- ggplot(data = df_long, aes(x = Group.1, y = Value*100, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(
    aes(label = round(Value*100, 2)), 
    position = position_stack(vjust = 0.5), 
    size = 3.5,
    color = "black") +
  labs(x = "",y = "Relative abundance",fill = "Type") +
  theme_bw() +
  theme(axis.text = element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold",color="black"),
        #axis.text.x = element_text(angle = 45, hjust = 1),  
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p1

#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#ggsave("B-gs.0.95.abundance.pdf",p1,height = 5, width = 6)

###RRN-OTUs####
setwd("~/Documents/sunR/Bacteria882/rrndb")
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")

#
blast <- read_xlsx("blast.uotus.bacteria.rrndb.xlsx",sheet=1)
blast0 <- blast %>% filter(pident>= 95 & qcovs>= 95)

metadata <- read.delim("rrnDB-5.8.tsv",sep = "\t",header = T)

blast1 <- blast0[blast0$qseqid %in% otu_table$OTU_id,]
blast2 <- blast1[!duplicated(blast1$qseqid),]  
blast0 <- concat.split(data = blast2,split.col = "sseqid",sep = "\\|",drop = TRUE)  
blast3 <- merge(blast0, metadata, by.x = "sseqid_2", by.y="Data.source.record.id")

fil_otu <- blast3$qseqid
keep_cols <- unique(c("Group.1", fil_otu))
selected_cols <- bact.phy.lev %>% dplyr::select(keep_cols)

selected_abundance <- rowSums(selected_cols[,-1])

df <- data.frame(
  Group.1 = c("Air:Control", "Leaf:Control", "Root:Control", "Leaf:N-addition", "Root:N-addition"),
  Proportion = c(0.5486751, 0.8501809, 0.8496693, 0.7732593, 0.8063286)
)
df_long1 <- df %>%
  mutate(Other = 1 - Proportion) %>%
  pivot_longer(cols = c(Proportion, Other), names_to = "Type", values_to = "Value")


p2 <- ggplot(data = df_long1, aes(x = Group.1, y = Value*100, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(
    aes(label = round(Value*100, 2)),  
    position = position_stack(vjust = 0.5),  
    size = 3.5,
    color = "black") +
  labs(x = "",y = "Relative abundance",fill = "Type") +
  theme_bw() +
  theme(axis.text = element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=18,face="bold",color="black"),
        #axis.text.x = element_text(angle = 45, hjust = 1),  
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p2


#setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#ggsave("B-rrn.0.95.abundance.pdf",p2,height = 5, width = 6)


##Supplementary.Figure13.Genomesize.GTDB.cv-abundance####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)
setwd("~/Documents/sunR/Bacteria882/GTDB")
metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)

metadata_bac_split <- metadata_bac %>%
  separate(
    col = gtdb_taxonomy, 
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"), 
    sep = ";",  
    remove = FALSE  
  )

result <- metadata_bac_split %>%
  group_by(genus) %>%
  summarise(
    mean = mean(genome_size),                
    stddev = sd(genome_size),                 
    cv = sd(genome_size)/mean(genome_size), 
    .groups = "drop"                   
  )

print(result)


group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)


otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")


blast <- read.table("blast.uotus.bacteria.gtdb.txt",comment.char = "")  # comment.char = "" 

metadata_bac <- read.delim("bac120_metadata_r214.tsv",sep = "\t",header = T)
metadata_arc <- read.delim("ar53_metadata_r214.tsv",sep = "\t",header = T)
metadata <- rbind(metadata_bac,metadata_arc)

blast1 <- blast[blast$V1 %in% otu_table$OTU_id,]

blast1 <- blast1 %>% separate(V5, into = c("First", "Second"), "~")
blast1 <- blast1[!duplicated(blast1$V1),]  


blast2 <- merge(blast1,metadata,by.x = "First",by.y="accession")
blast3 <- blast2 %>% dplyr::select("V1","gtdb_taxonomy","genome_size",
                                   "gtdb_genome_representative")

blast3_split <- blast3 %>%
  separate(
    col = gtdb_taxonomy,  
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),  
    sep = ";",  
    remove = FALSE
  )

data.cv <- merge(blast3_split, result, by="genus")
data.cv[is.na(data.cv)] <- 0


setwd("~/Documents/sunR/Bacteria882")
##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
ID[ID == "NA"] <- NA
id <- merge(ID, data.cv, by.x="names", by.y="V1")


otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% dplyr::select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(id) <- id$names
id <- id[otu_rowname,]

##2.genus
bact.phy <- aggregate(otu_all1,by=list(id$genus.y) , sum) 
rownames(bact.phy) <- bact.phy[,1]; 
bact.phy <- bact.phy[,-1] 
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

gs.genus <- data.frame(t(bact.phy0))
gs.genus$abundance <- rowMeans(gs.genus)*100
gs.genus$genus <- rownames(gs.genus)

gs.g1 <- gs.genus %>% dplyr::select("genus","abundance")

id.cv <- id %>% dplyr::select("genus.y","cv")
id.cv1 <- distinct_all(id.cv)
id.cv1$genus.y
id.cv2 <- id.cv1 %>% mutate(genus.y=str_replace_all(genus.y, pattern = "-", replacement = "."))

gs.cv <- left_join(gs.g1, id.cv2, by = c("genus" = "genus.y"))

library(ggplot2)

sorted_taxa <- gs.cv[order(-gs.cv$cv), "genus"]
gs.cv$genus <- factor(gs.cv$genus, levels = sorted_taxa)
# 2. 
left_axis_range <- c(0, 100)   
right_axis_range <- c(0, 6)   

# 3. 
scale_factor <- left_axis_range[2] / right_axis_range[2]
colors <- c('#5470C6', '#EE6666','#91CC75', '#EE6666')
# 4. 
combined_plot <- ggplot(gs.cv, aes(x = genus)) +
  geom_col(
    aes(y = cv*100),
    fill = colors[1],
    width = 0.7,
    position = position_nudge(x = -0.2)
  ) +
  geom_col(
    aes(y = abundance * scale_factor), 
    fill = colors[2],
    width = 0.7,
    position = position_nudge(x = 0.2)
  ) +
  geom_hline(yintercept = 25, color = "black", linetype = "solid", size = 0.2) +
  scale_y_continuous(
    name = "Coefficient of variation (%)",
    limits = left_axis_range,  
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / scale_factor,  
      name = "Relative abundance (%)",
      breaks = seq(right_axis_range[1], right_axis_range[2], by = 1)  
    )
  ) +
  labs(x = "Genus level") +
  theme(
    axis.title = element_text(size=15,face="bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    axis.line.x = element_line(color = "black"),
    axis.text.y = element_text(color = colors[1]),
    axis.ticks.y = element_line(color = colors[1]),
    axis.title.y = element_text(color = colors[1], face = "bold"),
    axis.line.y.left = element_line(color = colors[1]),
    axis.text.y.right = element_text(color = colors[2]),
    axis.ticks.y.right = element_line(color = colors[2]),
    axis.title.y.right = element_text(color = colors[2], face = "bold"),
    axis.line.y.right = element_line(color = colors[2]),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(combined_plot)
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-gs.gtdb.cv-abund.otus.pdf",combined_plot,height = 5, width = 11)


bact.phy0 <- bact.relabu.phy

bact.phy0 <- bact.phy0[,order(-colSums(bact.phy0))] ## 
lev <- interaction(group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.phy.lev <- aggregate(bact.phy0, by = list(lev) , FUN = "mean") # generate the mean for each factor level
bact.phy.lev <- data.frame(t(bact.phy.lev))
colnames(bact.phy.lev) <- bact.phy.lev[1, ]
bact.phy.lev <- bact.phy.lev[-1,]
bact.phy.lev_num <- bact.phy.lev %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()  
bact.phy.lev_num$genus <- rownames(bact.phy.lev_num)

library(dplyr)
filtered_genus <- gs.cv %>% mutate (cv = as.numeric (cv)) %>% 
  filter (!is.na (cv), cv < 0.25) 

total_abundance <- colSums(bact.phy.lev_num)

genus25 <- bact.phy.lev_num %>% filter(genus %in% filtered_genus$genus)

selected_abundance <- colSums(genus25[,-6])

proportion <- selected_abundance / total_abundance
print(proportion)


df <- data.frame(
  Group.1 = c("Air:Control", "Leaf:Control", "Root:Control", "Leaf:N-addition", "Root:N-addition"),
  Proportion = c(0.9128189, 0.9551786, 0.9397334, 0.9494956, 0.9408245)
)

library(ggplot2)
library(dplyr)
library(tidyr)

df_long <- df %>%
  mutate(Other = 1 - Proportion) %>%
  pivot_longer(cols = c(Proportion, Other), names_to = "Type", values_to = "Value")

colors <- c("Proportion" = "dodgerblue", "Other" = "gray")

p <- ggplot(df_long, aes(x = "", y = Value, fill = Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(label = scales::percent(Value)), 
            position = position_stack(vjust = 0.5),
            color = "white", size = 3) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Group.1) +
  scale_fill_manual(values = colors) +
  labs(title = "Relative abundance (CV_GS < 0.25)") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black")
  )

p

#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-gs.cv0.25.abundance.pdf",p,height = 5, width = 6)



##Supplementary.Figure14.rrn.RDP.Genus.CV-Abundance####
###rrn.CV####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)

#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/rrndb")
#rrndb <- readLines("rrnDB-5.8_16S_rRNA_new.fasta")
#rrndb1 <- readLines("rrnDB-5.8_16S_rRNA.fasta")
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

otu_table <- otu_table <- read.csv("../otu_bactFlattening825.csv")

blast <- read.table("blast.uotus.bacteria.rrndb.txt",comment.char = "") 
metadata <- read.delim("rrnDB-5.8.tsv",sep = "\t",header = T)

blast1 <- blast[blast$V1 %in% otu_table$OTU_id,] 
blast2 <- blast1[!duplicated(blast1$V1),] 
blast0 <- concat.split(data = blast2,split.col = "V5",sep = "\\|",drop = TRUE) 

blast3 <- merge(blast0, metadata, by.x = "V5_2", by.y="Data.source.record.id")

blast4 <- blast3 %>% dplyr::select("V1","RDP.taxa","X16S.gene.count","RDP.taxonomic.lineage")
blast5 <- concat.split(blast4, split.col="RDP.taxa",sep="\\(",drop=T)
blast5$RDP.taxa_2 <- gsub("\\)", "", blast5$RDP.taxa_2)

##RDP
rdp <- read_xlsx("rrnDB-5.8_pantaxa_stats_RDP.xlsx",sheet=2)

data <- merge(blast5,rdp, by.x="RDP.taxa_1",by.y="name")
data0 <- data[!duplicated(data$RDP.taxa_1),]

data_sorted <- data0 %>%
  arrange(desc(CV)) %>%       
  mutate(Index = 1:n())        

str(data_sorted$CV)
data_sorted$RDP.taxa_1 <- sub("[-_ ].+", "", data_sorted$RDP.taxa_1)


p <- ggplot(data_sorted, aes(x = Index, y = CV*100)) +
  geom_col(fill = "dodgerblue") +
  geom_segment(aes(x = 0, xend = 788, y = 0, yend = 0),color = "dodgerblue",linewidth = 0.7,linetype = "solid") +
  labs(title = "",
       x = "Genus level",y = "Coefficient of variation (%)") +
  scale_y_continuous(limits = c(0, 110),expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12,color="black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 15,color="black",face="bold")
  )
p

#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-rrn.rdp.cv.genus.pdf",p,height = 5, width = 10)

###genus.abundance####
library(vegan)
library(betapart)
library(colorRamps)
library(agricolae)
library(ggplot2)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
#library(learnasreml)
##Bacteria-AIR,ROOT,LEAF-7
setwd("~/Documents/sunR/Bacteria882")

#rm(list = ls())
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

blast6 <- blast5 %>% select("V1","RDP.taxa_1")

id <- merge(ID,blast6,by.x="names",by.y="V1",all.x=TRUE)
##2.genus
bact.phy <- aggregate(otu_all1,by=list(id$RDP.taxa_1) , sum) 
rownames(bact.phy) <- bact.phy[,1]; 
bact.phy <- bact.phy[,-1] 
data1 <- data.frame(t(bact.phy))

total1 <- apply(data1, 1, sum); 
bact.relabu.phy <- data.frame(lapply(data1, function(x) {  x / total1  }) )

bact.phy0 <- bact.relabu.phy

rdp.genus <- data.frame(t(bact.phy0))
rdp.genus$abundance <- rowMeans(rdp.genus)*100
rdp.genus$RDP.taxa_1 <- rownames(rdp.genus)
rdp.g1 <- rdp.genus %>% select("RDP.taxa_1","abundance")

rdp.g1$RDP.taxa_1 <- sub("[-_ ].+", "", rdp.g1$RDP.taxa_1)
rdp.data <- merge(data_sorted,rdp.g1,by="RDP.taxa_1",all.x=TRUE)

rdp.data$abundance[is.na(rdp.data$abundance)] <- 0

library(ggplot2)

sorted_taxa <- rdp.data[order(-rdp.data$CV), "V1"]$V1
rdp.data$V1 <- factor(rdp.data$V1, levels = sorted_taxa)
# 2. 
left_axis_range <- c(0, 105)   
right_axis_range <- c(0, 3.5)   

# 3. 
scale_factor <- left_axis_range[2] / right_axis_range[2]
colors <- c('#5470C6', '#EE6666','#91CC75', '#EE6666')
# 4.
combined_plot <- ggplot(rdp.data, aes(x = V1)) +
  geom_col(
    aes(y = CV * 100),
    fill = colors[1],
    width = 0.7,
    position = position_nudge(x = -0.2)
  ) +
  geom_col(
    aes(y = abundance * scale_factor),  
    fill = colors[2],
    width = 0.7,
    position = position_nudge(x = 0.2)
  ) +
  geom_hline(yintercept = 25, color = "black", linetype = "solid", size = 0.2) +
  scale_y_continuous(
    name = "Coefficient of variation (%)",
    limits = left_axis_range, 
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / scale_factor,  
      name = "Relative abundance (%)",
      breaks = seq(right_axis_range[1], right_axis_range[2], by = 1)  
    )
  ) +
  labs(x = "Genus level") +
  theme(
    axis.title = element_text(size=15,face="bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = NA),
    axis.line.x = element_line(color = "black"),
    axis.text.y = element_text(color = colors[1]),
    axis.ticks.y = element_line(color = colors[1]),
    axis.title.y = element_text(color = colors[1], face = "bold"),
    axis.line.y.left = element_line(color = colors[1]),
    axis.text.y.right = element_text(color = colors[2]),
    axis.ticks.y.right = element_line(color = colors[2]),
    axis.title.y.right = element_text(color = colors[2], face = "bold"),
    axis.line.y.right = element_line(color = colors[2]),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(combined_plot)
#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-rrn.rdp.cv-abund.genus.pdf",combined_plot,height = 5, width = 11)







##Supplementary.Table1-6;Figure17.Treatment.Genomic.OTU.Function####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)
library(lme4)
library(lmerTest)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(ggh4x)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
library(gridExtra)
##Genomic.traits####
rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
###genome.size####
genomesize <- read_xlsx("Bact_genomesize.xlsx",sheet=1)

genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","Treatment","sampleSite","value")]
genomesize0 <- na.omit(genomesize0)

model <- lmer(value/1000000 ~ TP0*sampleType*Treatment + (1|sampleSite), data = genomesize0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Genome size - res")
qqnorm(res); qqline(res)   
par(mfrow=c(1,2))

shapiro.test(res)

plot(resid(model) ~ genomesize0$TP0,
     xlab="Genome size$Time", ylab="Model Residuals")
abline(h = 0, col = "red")

###air.genome.lmer.residuals####
air <- genomesize0[genomesize0$sampleType=="Air",]
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)    
shapiro.test(res)

###leaf.genome.lmer.residuals####
leaf <- genomesize0[genomesize0$sampleType=="Leaf",]
model <- lmer(value/1000000 ~ TP0*Treatment + (1|sampleSite), data = leaf)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Genome size (Leaf) - res")
qqnorm(res); qqline(res)    
par(mfrow=c(1,2))
shapiro.test(res)


leaf1 <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = leaf1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)


leaf2 <- genomesize0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = leaf2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)


###root.genome.lmer.residuals####
root <- genomesize0[genomesize0$sampleType=="Root",]
model <- lmer(value/1000000 ~ TP0*Treatment + (1|sampleSite), data = root)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Genome size (Root) - res")
qqnorm(res); qqline(res)   
par(mfrow=c(1,2))
shapiro.test(res)


root1 <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="Control")
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = root1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)  
shapiro.test(res)

root2 <- genomesize0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = root2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)



###rrnDB####
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
rrna <- read_xlsx("Bact_rrn.xlsx",sheet=1)

rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
rrna0 <- na.omit(rrna0)


model3 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = rrna0)
summary(model3)
anova(model3)

plot(resid(model3) ~ rrna0$TP0,
     xlab="RRN$Time", ylab="Model Residuals")
abline(h = 0, col = "red")

res <- resid(model3, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN - res")
qqnorm(res); qqline(res)  
shapiro.test(res)

###air.RRN.lmer.residuals####
air <- rrna0[rrna0$sampleType=="Air",]
model <- lmer(value ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)    
shapiro.test(res)


###leaf.RRN.lmer.residuals####
leaf <- rrna0[rrna0$sampleType=="Leaf",]
model <- lmer(value ~ TP0*Treatment + (1|sampleSite), data = leaf)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN (Leaf) - res")
qqnorm(res); qqline(res)   
par(mfrow=c(1,2))
shapiro.test(res)


leaf1 <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
model <- lmer(value ~ TP0 + (1|sampleSite), data = leaf1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)


leaf2 <- rrna0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
model <- lmer(value ~ TP0 + (1|sampleSite), data = leaf2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)    
shapiro.test(res)


###root.RRN.lmer.residuals####
root <- rrna0[rrna0$sampleType=="Root",]
model <- lmer(value ~ TP0*Treatment + (1|sampleSite), data = root)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="RRN (Root) - res")
qqnorm(res); qqline(res)    
par(mfrow=c(1,2))
shapiro.test(res)


root1 <- rrna0 %>% filter(sampleType=="Root"&Treatment=="Control")
model <- lmer(value ~ TP0 + (1|sampleSite), data = root1)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)

root2 <- rrna0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
model <- lmer(value ~ TP0 + (1|sampleSite), data = root2)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Residual")
qqnorm(res); qqline(res)   
shapiro.test(res)


###OTU####
setwd("~/Documents/sunR/Bacteria882")
#rm(list = ls())
##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","sampleType","TP00","TP0","Treatment","sampleSite")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
bact.otu <- aggregate(otu_all1,by=list(ID$names) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )
bact.otu0 <- bact.relabu.otu
bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))]

colnames(bact.otu0)
bact.otu10 <- bact.otu0[, c(1:10)]
bact.otu10$sample.id <- rownames(bact.otu10)

b.otu10 <- merge(bact.otu10, group2, by="sample.id")

b.otu <- gather(b.otu10,key = "OTU",
                value = "Abundance",-"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment",-"sampleSite")


b.otu$value <- b.otu$Abundance*100

otu_names <- c("OTU5","OTU2","OTU21","OTU42","OTU7",
               "OTU18","OTU36","OTU47","OTU15","OTU12") 

for(otu_name in otu_names) {
  otu_data <- b.otu[b.otu$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * sampleType * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  #print(summary(model))
  
  print(anova(model))
}

###root.otu.lmer####
root <- b.otu[b.otu$sampleType=="Root",]

for(otu_name in otu_names) {
  otu_data <- root[root$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}



root1 <- b.otu %>% filter(sampleType=="Root"&Treatment=="Control")

for(otu_name in otu_names) {
  otu_data <- root1[root1$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}



root2 <- b.otu %>% filter(sampleType=="Root"&Treatment=="N-addition")

for(otu_name in otu_names) {
  otu_data <- root2[root2$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}


###leaf.otu.lmer####
leaf <- b.otu[b.otu$sampleType=="Leaf",]

for(otu_name in otu_names) {
  otu_data <- leaf[leaf$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}


leaf1 <- b.otu %>% filter(sampleType=="Leaf"&Treatment=="Control")

for(otu_name in otu_names) {
  otu_data <- leaf1[leaf1$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}


leaf2 <- b.otu %>% filter(sampleType=="Leaf"&Treatment=="N-addition")

for(otu_name in otu_names) {
  otu_data <- leaf2[leaf2$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}

###air.otu.lmer####
air <- b.otu[b.otu$sampleType=="Air",]

for(otu_name in otu_names) {
  otu_data <- air[air$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")
  
  print(anova(model))
}





###Functions####
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
data <- read_xlsx("function_proportion.xlsx", sheet=2)
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","Treatment","sampleSite")]
group2 <- na.omit(group1)

data0 <- merge(data,group2,by="sample.id")

func_names <- c("Nitrogen fixation",
                "Amino acid metabolism",
                "Transcription",
                "Energy metabolism",
                "Metabolism of terpenoids and polyketides",
                "Carbohydrate metabolism",
                "Nucleotide metabolism",
                "Xenobiotics biodegradation and metabolism",
                
                "Signal transduction",
                "Glycan biosynthesis and metabolism",
                "Transport and catabolism",
                "Cell motility")  

for(func_name in func_names) {
  func_data <- data0[data0$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * sampleType * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}

###root.function.lmer####
root <- data0[data0$sampleType=="Root",]
for(func_name in func_names) {
  func_data <- root[root$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}


root1 <- data0 %>% filter(sampleType=="Root"&Treatment=="Control")
for(func_name in func_names) {
  func_data <- root1[root1$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}


root2 <- data0 %>% filter(sampleType=="Root"&Treatment=="N-addition")
for(func_name in func_names) {
  func_data <- root2[root2$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}



###leaf.function.lmer####
leaf <- data0[data0$sampleType=="Leaf",]
for(func_name in func_names) {
  func_data <- leaf[leaf$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}

leaf1 <- data0 %>% filter(sampleType=="Leaf"&Treatment=="Control")
for(func_name in func_names) {
  func_data <- leaf1[leaf1$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}


leaf2 <- data0 %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
for(func_name in func_names) {
  func_data <- leaf2[leaf2$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}


###air.function.lmer####
air <- data0[data0$sampleType=="Air",]
for(func_name in func_names) {
  func_data <- air[air$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")
  
  print(anova(model))
}





##Supplementary.Table##########
##Supplementary.Table.Treatment.Genomic.OTU.Function####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)
library(lme4)
library(lmerTest)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(ggh4x)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
library(gridExtra)
###Genomic.traits####
#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
###genome size####
genomesize <- read_xlsx("Bact_genomesize.xlsx",sheet=1)

genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","Treatment","sampleSite","value")]
genomesize0 <- na.omit(genomesize0)

model <- lmer(value/1000000 ~ TP0*sampleType*Treatment + (1|sampleSite), data = genomesize0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="residual")
qqnorm(res); qqline(res)    # 
shapiro.test(res)

# residuals vs. TP0
par(mfrow=c(1,1))
plot(resid(model) ~ genomesize0$TP0,
     xlab="Genome size$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
#acf(resid(model))


air <- genomesize0[genomesize0$sampleType=="Air",]
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)
shapiro.test(resid(model, type = "pearson"))

#sampleType
sampleTypes <- c("Root", "Leaf")

#for-sampleType
for(sampleType in sampleTypes) {
  subset_data <- genomesize0[genomesize0$sampleType == sampleType, ]
  
  model <- lmer(value/1000000 ~ TP0 * Treatment + (1 | sampleSite), data = subset_data)
  
  cat("###", sampleType, "###\n")
  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

root <- genomesize0 %>% filter(sampleType=="Root")
Treatments <- c("Control","N-addition")
for(Treatment in Treatments) {
  # 
  subset_data <- root[root$Treatment == Treatment, ]
  # 
  model <- lmer(value/1000000 ~ TP0 + (1 | sampleSite), data = subset_data)
  #
  cat("###", Treatment, "###\n")
  #
  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

leaf <- genomesize0 %>% filter(sampleType=="Leaf")
Treatments <- c("Control","N-addition")
for(Treatment in Treatments) {
  subset_data <- leaf[leaf$Treatment == Treatment, ]
  model <- lmer(value/1000000 ~ TP0 + (1 | sampleSite), data = subset_data)
  cat("###", Treatment, "###\n")
  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}



###rrnDB####
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
rrna <- read_xlsx("Bact_rrn.xlsx",sheet=1)

rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
rrna0 <- na.omit(rrna0)


model3 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = rrna0)
summary(model3)
anova(model3)

res <- resid(model3, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="residual")
qqnorm(res); qqline(res)   #Q-Q plot
shapiro.test(res)

par(mfrow=c(1,1))
plot(resid(model3) ~ genomesize0$TP0,
     xlab="RRN$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
#acf(resid(model))


air <- rrna0[rrna0$sampleType=="Air",]
model <- lmer(value ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)
shapiro.test(resid(model, type = "pearson"))

sampleTypes <- c("Root", "Leaf") 

for(sampleType in sampleTypes) {
  subset_data <- rrna0[rrna0$sampleType == sampleType, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = subset_data)
  cat("###", sampleType, "###\n")
  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

root <- rrna0 %>% filter(sampleType=="Root")
Treatments <- c("Control","N-addition")
for(Treatment in Treatments) {
  subset_data <- root[root$Treatment == Treatment, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = subset_data)
  
  cat("###", Treatment, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

leaf <- rrna0 %>% filter(sampleType=="Leaf")
Treatments <- c("Control","N-addition")
for(Treatment in Treatments) {
  subset_data <- leaf[leaf$Treatment == Treatment, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = subset_data)
  
  cat("###", Treatment, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}







###OTU####
setwd("~/Documents/sunR/Bacteria882")
#rm(list = ls())
##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),] 
group1 <- group[,c("sample.id","sampleType","TP00","TP0","Treatment","sampleSite")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% dplyr::select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
bact.otu <- aggregate(otu_all1,by=list(ID$names) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )
#rowSums(bact.relabu.otu[,1:5211])
bact.otu0 <- bact.relabu.otu
bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))]

colnames(bact.otu0)
bact.otu10 <- bact.otu0[, c(1:10)]
bact.otu10$sample.id <- rownames(bact.otu10)

b.otu10 <- merge(bact.otu10, group2, by="sample.id")

b.otu <- gather(b.otu10,key = "OTU",
                value = "Abundance",-"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment",-"sampleSite")


b.otu$value <- b.otu$Abundance*100

###OTU.residual-time####
otu <- b.otu %>% filter(OTU=="OTU5")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

# residuals vs. TP0
plot(resid(model4) ~ otu$TP0,
     xlab="OTU5$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU2")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

# residuals vs. TP0
plot(resid(model4) ~ otu$TP0,
     xlab="OTU2$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")

otu <- b.otu %>% filter(OTU=="OTU21")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

# residuals vs. TP0
plot(resid(model4) ~ otu$TP0,
     xlab="OTU21$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU42")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU42$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")

otu <- b.otu %>% filter(OTU=="OTU7")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU7$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU18")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU18$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU36")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU36$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU47")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU47$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU15")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU15$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


otu <- b.otu %>% filter(OTU=="OTU12")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = otu)
summary(model4)
anova(model4)

plot(resid(model4) ~ otu$TP0,
     xlab="OTU12$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


###OTU.lmer(3)####
#OTU
otu_names <- c("OTU5","OTU2","OTU21","OTU42","OTU7",
               "OTU18","OTU36","OTU47","OTU15","OTU12") 

#for-OTU
for(otu_name in otu_names) {
  otu_data <- b.otu[b.otu$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * sampleType * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

root <- b.otu[b.otu$sampleType=="Root",]

for(otu_name in otu_names) {
  otu_data <- root[root$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

root1 <- b.otu %>% filter(sampleType=="Root"&Treatment=="Control")
for(otu_name in otu_names) {
  otu_data <- root1[root1$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}
root2 <- b.otu %>% filter(sampleType=="Root"&Treatment=="N-addition")
for(otu_name in otu_names) {
  otu_data <- root2[root2$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}


leaf <- b.otu[b.otu$sampleType=="Leaf",]

for(otu_name in otu_names) {
  otu_data <- leaf[leaf$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}
leaf1 <- leaf %>% filter(sampleType=="Leaf"&Treatment=="Control")
for(otu_name in otu_names) {
  otu_data <- leaf1[leaf1$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

leaf2 <- leaf %>% filter(sampleType=="Leaf"&Treatment=="N-addition")
for(otu_name in otu_names) {
  otu_data <- leaf2[leaf2$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}



air <- b.otu[b.otu$sampleType=="Air",]

for(otu_name in otu_names) {
  otu_data <- air[air$OTU == otu_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = otu_data)
  
  cat("###", otu_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}


###Functions####
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
data <- read_xlsx("function_proportion.xlsx", sheet=2)
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),]  
group1 <- group[,c("sample.id","Treatment","sampleSite")]
group2 <- na.omit(group1)

data0 <- merge(data,group2,by="sample.id")

###function.residual-time####
func <- data0 %>% filter(level2=="Nitrogen fixation")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Nitrogen fixation$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Amino acid metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Amino acid metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Transcription")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Transcription$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Energy metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Energy metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Metabolism of terpenoids and polyketides")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="M.terpenoids.polyketides$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Carbohydrate metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Carbohydrate.metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Nucleotide metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Nucleotide metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Xenobiotics biodegradation and metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Xenobiotics.BD.metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Signal transduction")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Signal transduction$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Glycan biosynthesis and metabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Glycan.biosyn.metabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Transport and catabolism")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Transport and catabolism$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")


func <- data0 %>% filter(level2=="Cell motility")

model4 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = func)
summary(model4)
anova(model4)

plot(resid(model4) ~ func$TP0,
     xlab="Cell motility$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model4)); qqline(resid(model4), col="red")

###function.lmer(3)####
func_names <- c("Nitrogen fixation",
                "Amino acid metabolism",
                "Transcription",
                "Energy metabolism",
                "Metabolism of terpenoids and polyketides",
                "Carbohydrate metabolism",
                "Nucleotide metabolism",
                "Xenobiotics biodegradation and metabolism",
                
                "Signal transduction",
                "Glycan biosynthesis and metabolism",
                "Transport and catabolism",
                "Cell motility")  


for(func_name in func_names) {
  func_data <- data0[data0$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * sampleType * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}


root <- data0[data0$sampleType=="Root",]
for(func_name in func_names) {
  func_data <- root[root$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

root1 <- root %>% filter(Treatment=="Control")
for(func_name in func_names) {
  func_data <- root1[root1$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}
root2 <- root %>% filter(Treatment=="N-addition")
for(func_name in func_names) {
  func_data <- root2[root2$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}


leaf <- data0[data0$sampleType=="Leaf",]
for(func_name in func_names) {
  func_data <- leaf[leaf$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

leaf1 <- leaf %>% filter(Treatment=="Control")
for(func_name in func_names) {
  func_data <- leaf1[leaf1$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}

leaf2 <- leaf %>% filter(Treatment=="N-addition")
for(func_name in func_names) {
  func_data <- leaf2[leaf2$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}


air <- data0[data0$sampleType=="Air",]
for(func_name in func_names) {
  func_data <- air[air$level2 == func_name, ]
  
  model <- lmer(value ~ TP0 + (1 | sampleSite), data = func_data)
  
  cat("###", func_name, "###\n")

  print(anova(model))
  print(shapiro.test(resid(model, type = "pearson")))
}




##Supplementary.Table6.genomic traits.model1 vs model2--temporal autocorrelation####
library(vegan)
library(betapart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(splitstackshape)
library(patchwork)
library(colorRamps)
library(readxl)
library(lme4)
library(lmerTest)
library(writexl)
library(stringr)
library(tidyverse)
library(reshape2)
library(ggbreak)
library(ggh4x)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
library(gridExtra)
###Genomic.traits####
rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
#
genomesize <- read_xlsx("Bact_genomesize.xlsx",sheet=1)

genomesize0 <- genomesize[, c("sample.id","sampleType","TP0","TP00","Treatment","sampleSite","value")]
genomesize0 <- na.omit(genomesize0)

#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC修改图")
#write.csv(genomesize0, "GS.csv")

model <- lmer(value/1000000 ~ 
                TP0*sampleType*Treatment + (1|sampleSite), 
              data = genomesize0)
summary(model)
anova(model)

plot(resid(model) ~ genomesize0$TP0,
     xlab="Genome size$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model)); qqline(resid(model), col="red")


air <- genomesize0[genomesize0$sampleType=="Air",]
model <- lmer(value/1000000 ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)
sampleTypes <- c("Root", "Leaf") 

for(sampleType in sampleTypes) {
  subset_data <- genomesize0[genomesize0$sampleType == sampleType, ]
  model <- lmer(value/1000000 ~ TP0 * Treatment + (1 | sampleSite), data = subset_data)
  
  cat("###", sampleType, "###\n")
  print(anova(model))
}

library(nlme)
##root####
root <- genomesize0 %>% filter(sampleType=="Root")
root1 <- root %>% filter(Treatment=="Control")
root2 <- root %>% filter(Treatment=="N-addition")

root_mean <- root1 %>%
  group_by(sampleSite, TP0) %>%  
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()  # 
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value/1000000 ~ TP0,  
  random = ~1|sampleSite,   
  correlation = corAR1(form = ~TP0|sampleSite),   
  data = root_mean,
  method = "REML"  # 
)
# 
summary(model_with_autocor)
anova(model_with_autocor)
# 
model_no_autocor <- lme(
  mean_value/1e6 ~ TP0,
  random = ~1|sampleSite,
  data = root_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)
# 
model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)   



root_mean <- root2 %>%
  group_by(sampleSite, TP0) %>%   
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()   
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value/1000000 ~ TP0,   
  random = ~1|sampleSite,   
  correlation = corAR1(form = ~TP0|sampleSite),   
  data = root_mean,
  method = "REML"   
)
# 
summary(model_with_autocor)
anova(model_with_autocor)
# 
model_no_autocor <- lme(
  mean_value/1e6 ~ TP0,
  random = ~1|sampleSite,
  data = root_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)
# 
model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  


##leaf####
leaf <- genomesize0 %>% filter(sampleType=="Leaf")
leaf1 <- leaf %>% filter(Treatment=="Control")
leaf2 <- leaf %>% filter(Treatment=="N-addition")

leaf_mean <- leaf1 %>%
  group_by(sampleSite, TP0) %>%   
  summarise(
    mean_value = mean(value, na.rm = TRUE),   
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value/1000000 ~ TP0,   
  random = ~1|sampleSite,   
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = leaf_mean,
  method = "REML"  
)
# 
summary(model_with_autocor)
anova(model_with_autocor)
# 
model_no_autocor <- lme(
  mean_value/1e6 ~ TP0,
  random = ~1|sampleSite,
  data = leaf_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)
# 
model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  



leaf_mean <- leaf2 %>%
  group_by(sampleSite, TP0) %>%  
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value/1000000 ~ TP0, 
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = leaf_mean,
  method = "REML"  
)
# 
summary(model_with_autocor)
anova(model_with_autocor)
# 
model_no_autocor <- lme(
  mean_value/1e6 ~ TP0,
  random = ~1|sampleSite,
  data = leaf_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)
# 
model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)   


##air####
air <- genomesize0 %>% filter(sampleType=="Air")

air_mean <- air %>%
  group_by(sampleSite, TP0) %>% 
  summarise(
    mean_value = mean(value, na.rm = TRUE), 
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value/1000000 ~ TP0,  
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = air_mean,
  method = "REML" 
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value/1e6 ~ TP0,
  random = ~1|sampleSite,
  data = air_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  


###rrnDB####
setwd("~/Documents/sunR/Bacteria882/n.treatment/re-genomic")
rrna <- read_xlsx("Bact_rrn.xlsx",sheet=1)

rrna0 <- rrna[, c("sample.id","sampleType","TP0","TP00","sampleSite","Treatment","value")]
rrna0 <- na.omit(rrna0)

model3 <- lmer(value ~ TP0*sampleType*Treatment + (1|sampleSite), 
               data = rrna0)
summary(model3)
anova(model3)

plot(resid(model3) ~ rrna0$TP0,
     xlab="RRN$Time", ylab="Model Residuals")
abline(h = 0, col = "red")
qqnorm(resid(model3)); qqline(resid(model3), col="red")


air <- rrna0[rrna0$sampleType=="Air",]
model <- lmer(value ~ TP0 + (1|sampleSite), data = air)
summary(model)
anova(model)


sampleTypes <- c("Root", "Leaf") 

for(sampleType in sampleTypes) {
  subset_data <- rrna0[rrna0$sampleType == sampleType, ]
  
  model <- lmer(value ~ TP0 * Treatment + (1 | sampleSite), data = subset_data)
  
  cat("###", sampleType, "###\n")

  print(anova(model))
}

library(nlme)

##root####
root <- rrna0 %>% filter(sampleType=="Root")
root1 <- root %>% filter(Treatment=="Control")
root2 <- root %>% filter(Treatment=="N-addition")

root_mean <- root1 %>%
  group_by(sampleSite, TP0) %>%  
  summarise(
    mean_value = mean(value, na.rm = TRUE), 
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value ~ TP0,  
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite), 
  data = root_mean,
  method = "REML" 
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value ~ TP0,
  random = ~1|sampleSite,
  data = root_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  



root_mean <- root2 %>%
  group_by(sampleSite, TP0) %>% 
  summarise(
    mean_value = mean(value, na.rm = TRUE),  # mean, na
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value ~ TP0,  
  random = ~1|sampleSite, 
  correlation = corAR1(form = ~TP0|sampleSite),
  data = root_mean,
  method = "REML" 
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value ~ TP0,
  random = ~1|sampleSite,
  data = root_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  


##leaf####
leaf <- rrna0 %>% filter(sampleType=="Leaf")
leaf1 <- leaf %>% filter(Treatment=="Control")
leaf2 <- leaf %>% filter(Treatment=="N-addition")

leaf_mean <- leaf1 %>%
  group_by(sampleSite, TP0) %>% 
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value ~ TP0,  
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = leaf_mean,
  method = "REML"  
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value ~ TP0,
  random = ~1|sampleSite,
  data = leaf_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML) 



leaf_mean <- leaf2 %>%
  group_by(sampleSite, TP0) %>% 
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value ~ TP0,  
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = leaf_mean,
  method = "REML"  
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value ~ TP0,
  random = ~1|sampleSite,
  data = leaf_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML)  


##air####
air <- rrna0 %>% filter(sampleType=="Air")

air_mean <- air %>%
  group_by(sampleSite, TP0) %>%  
  summarise(
    mean_value = mean(value, na.rm = TRUE),  
    n = n()  
  ) %>%
  ungroup()

model_with_autocor <- lme(
  mean_value ~ TP0,  
  random = ~1|sampleSite,  
  correlation = corAR1(form = ~TP0|sampleSite),  
  data = air_mean,
  method = "REML"  
)

summary(model_with_autocor)
anova(model_with_autocor)

model_no_autocor <- lme(
  mean_value ~ TP0,
  random = ~1|sampleSite,
  data = air_mean,
  method = "REML"
)
summary(model_no_autocor)
anova(model_no_autocor)

model_with_autocor_ML <- update(model_with_autocor, method = "ML")
model_no_autocor_ML <- update(model_no_autocor, method = "ML")
anova(model_with_autocor_ML, model_no_autocor_ML) 





##Others############
##Supplementary.Figure.Alpha.diversity--QQplot####
library(vegan)
library(picante)
library(ggplot2)
library(car)
#library(ggcor)
library(ggpubr)
library(splitstackshape)
library(dplyr)
library(stringr)
library(broom)
library(agricolae)  
library(readxl)
library(patchwork)
##Bacteria-air7, root, leaf
rm(list = ls())
setwd("~/Documents/sunR/Bacteria882")
group <- read_xlsx("sample_metadata.xlsx", sheet = 3)

group <- group[order(group$sample.id), ]
group1 <- group[, c("sample.id","TP00","TP0","sampleType","sampleSite","Treatment")]
group2 <- group1 %>% filter(Treatment=="Control")
otu_all0 <- read.csv("otu_bactFlattening825.csv", check.names = F, row.names = 1, header = T)
otu_all <- otu_all0 %>% select(all_of(group2$sample.id))

otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -702]

otu_all1$OTU.ID <- rownames(otu_all1)
otu_all1 <- otu_all1[, 1:701] 
otu_all1 <- t(otu_all1) 

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = "shannon", base = base)
  Simpson <- diversity(x, index = "simpson") 
  Pielou <- Shannon / log(ACE, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- "PD_whole_tree"
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}

alpha_all <- alpha(otu_all1, base = 2)
alpha_all$sample.id <- rownames(alpha_all)

sum(alpha_all$sample.id %in% group2$sample.id == FALSE) 
which(alpha_all$sample.id %in% group2$sample.id == FALSE)

alpha_all0 <- merge(alpha_all, group2, by = "sample.id", all.x = TRUE) 

####alpha.diversity.Shannon.temporal.no-N####
library(lmerTest)
model <- lmer(Shannon ~ TP0*sampleType + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Shannon Index - res")
qqnorm(res); qqline(res)                       
par(mfrow=c(1,2))
shapiro.test(res)


alpha_root <- alpha_all0[alpha_all0$sampleType == "Root", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_root)))

x <- alpha_root$TP0
y <- alpha_root$Shannon
cor.test(x, y, method = "spearman")


alpha_leaf <- alpha_all0[alpha_all0$sampleType == "Leaf", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_leaf)))

x <- alpha_leaf$TP0
y <- alpha_leaf$Shannon
cor.test(x, y, method = "spearman")

alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Shannon
cor.test(x, y, method = "spearman")


p.sh <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Shannon, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Shannon, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Shannon index")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0.5, 10)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=2, label="Air: R = -0.334, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=1.4, label="Leaf: R = 0.194, P = 0.016", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=0.8, label="Root: R = 0.562, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") 


p.sh


p1 <- p.sh + labs(subtitle = "Time: F = 7.100, P = 0.008 \nCompartment: F = 187.389, P < 0.001 \nTime × Compartment: F = 22.794, P < 0.001")
p1 <- p1 + theme(plot.subtitle = element_text(size = 14,face = "bold",color="black")) 
p1

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.shannon.pdf",p1,height = 5, width = 5)

kw <- kruskal(alpha_all0$Shannon, alpha_all0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 430.0052, df = 2, p-value = 0
plot(kw)

p_shannon <- ggplot(alpha_all0, aes(x = sampleType, y = Shannon)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Shannon index") + 
  geom_text(x=2,y=2,label=expression(paste(chi^2, " = ", 430.005, ", ", "Df = ", 2, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=9.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=8,label="c",size=7,color="black")+
  annotate("text",x=3,y=8,label="b",size=7,color="black")+
  ylim(0.5, 10)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_shannon

####alpha.diversity.Richness.temporal.no-N####
model <- lmer(Richness ~ TP0*sampleType + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)

res <- resid(model, type = "pearson")
par(mfrow=c(1,2))
hist(res, breaks=30, main="Histogram of residuals", xlab="Richness - res")
qqnorm(res); qqline(res)                      
par(mfrow=c(1,2))
shapiro.test(res)

alpha_root <- alpha_all0[alpha_all0$sampleType == "Root", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_root)))

x <- alpha_root$TP0
y <- alpha_root$Richness
cor.test(x, y, method = "spearman")


alpha_leaf <- alpha_all0[alpha_all0$sampleType == "Leaf", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_leaf)))

x <- alpha_leaf$TP0
y <- alpha_leaf$Richness
cor.test(x, y, method = "spearman")

alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Richness
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")


p.rich <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Richness, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Richness, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "OTU richness")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0,770)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=120, label="Air: R = -0.396, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=70, label="Leaf: R = 0.285, P < 0.001", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=20, label="Root: R = 0.697, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") 


p.rich


p2 <- p.rich + labs(subtitle = "Time: F = 10.658, P = 0.001 \nCompartment: F = 646.375, P < 0.001 \nTime × Compartment: F = 65.233, P < 0.001")
p2 <- p2 + theme(plot.subtitle = element_text(size = 14,face = "bold",color="black")) 
p2

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.richness.pdf",p2,height = 5, width = 5)


kw <- kruskal(alpha_all0$Richness, alpha_all0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 524.9096, df = 2, p-value = 0
plot(kw)

p_rich <- ggplot(alpha_all0, aes(x = sampleType, y = Richness)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "OTU richness") + 
  geom_text(x=2,y=20,label=expression(paste(chi^2, " = ", 524.910, ", ", "Df = ", 2, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=760,label="a",size=7,color="black")+
  annotate("text",x=2,y=470,label="c",size=7,color="black")+
  annotate("text",x=3,y=470,label="b",size=7,color="black")+
  ylim(0, 770)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_rich


p <- p1+p_shannon+p2+p_rich
p

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.shannon.richness.non.01.revised.pdf",p,height = 10, width = 10)




##Figure9.data####
###data.functions####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)  # https://github.com/kwstat/pals/issues/3
library(RColorBrewer)
library(readxl)
library(ggpmisc)

rm(list=ls())
##
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)
group2 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group2)
group3 <- group2 %>% filter(Treatment == "Control")

# 
kegg <- read.table("kegg htext.txt", sep = "\t",fill = TRUE,header = T,quote = "")

setwd("~/Documents/sunR/Bacteria882/new-picrust2/picrust2_out/KO_metagenome_out")
#
ko_abundance <- read.table("pred_metagenome_unstrat.tsv", header = T, check.names = F)

# 
ko_abundance <- ko_abundance[,colnames(ko_abundance) %in% c("function",group3$sample.id)]
abundance = ko_abundance %>% column_to_rownames("function")  
ko_abundance <-ko_abundance[rowSums(abundance) != 0,] 

# 
ko_abundance2 <- merge(kegg,ko_abundance,by.x = "KO",by.y="function")
table(duplicated(paste0(ko_abundance2$pathway_id,ko_abundance2$KO)))  

# 
ko_abundance3 <- ko_abundance2[,c("pathway_id",group3$sample.id)]  
ko_abundance4 <- aggregate(. ~ pathway_id, data = ko_abundance3, FUN = sum)  

# 
ko_abundance5 <- merge(ko_abundance4,kegg[,c("pathway_id","level1","level2","level3")],
                       by.x="pathway_id",by.y="pathway_id")  
table(duplicated(ko_abundance5$pathway_id))  

ko_abundance5 <- ko_abundance5[-which(duplicated(ko_abundance5$pathway_id)),]

# 
ko_abundance5 <- ko_abundance5 %>%
  filter(level1 != "Human Diseases" & level1 != "Organismal Systems"& level2 != "Cellular community - eukaryotes") %>%
  select(-level1, -level3)  

# 
ko_abundance6 <- aggregate(.~level2,ko_abundance5[,2:703],FUN="sum") 
ko_abundance6[, 2:702]  <- apply(ko_abundance6[, 2:702], 2, function(x) x / sum(x))  
#rownames(ko_abundance6) <- ko_abundance6$level2
#write.csv(ko_abundance6, "ko_abundance6_air7.csv")
#write.csv(group2, "group_ko_air7.csv")

# read the data
f_table <- ko_abundance6
f_table

f_grp <- group3


# get the level index
f_lev <- f_table$level2

# create an data.frame for splitted result
f_spl <- data.frame()

# for loop
for(i in f_lev) {
  
  df_tmp <- f_table %>% filter(level2 == i) %>% select(-1) %>% t() %>% as.data.frame()
  colnames(df_tmp) <- "value"
  
  df_tmp <- df_tmp %>% mutate(sample.id = rownames(df_tmp)) %>% 
    mutate(level2 = i) %>% left_join(f_grp, by = "sample.id")
  
  f_spl <- rbind(f_spl, df_tmp)
  
}

f_spl

#setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#write.csv(f_spl, "function_proportion_n.csv")

p_all <- ggplot() + 
  geom_jitter(data=f_spl,  height = 0, 
              aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =2.6) + 
  facet_wrap(~level2, scales = "free") +
  geom_smooth(data=f_spl,
              aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 1, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("red", "blue","black"), name="Time")+
  labs(x = "Time",y = "Mean proportion(%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        #legend.position = "none", 
        axis.text=element_text(size=8),
        axis.title=element_text(size=12,face="bold")) 

p_all

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-function.no.n.lm.pdf", p_all, height = 15, width = 22)

###data.nitreogen.fixation####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
# 
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

group2 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group2)
group3 <- group2 %>% filter(Treatment == "Control")


ko <- read_xlsx("ko_gene.xlsx", sheet = 1)
ko_gene <- read_xlsx("ko_gene.xlsx", sheet=2)
ko_ncycle <- read_xlsx("ko_gene.xlsx",sheet=3)
ko.n <- read_xlsx("ko_gene.xlsx",sheet=4)

ko1 <- ko[,colnames(ko) %in% c("KO",group3$sample.id)]
# 
#ko4 <- aggregate(.~KO,ko,FUN="sum") 
ko1[, 2:702]  <- apply(ko1[, 2:702], 2, function(x) x / sum(x)) 

ko2 <- ko1[ko1$KO %in% ko.n$KO,]

nsums <- colSums(ko2[,2:702])
ko3 <- data.frame(FirstColumn = "Nitrogen fixation",nsums)
ko3$sample.id <- rownames(ko3)
ko3 <- merge(ko3, group3, by="sample.id")

#write.csv(ko3,"nitrogen.fixation.control.csv")

root <- ko3[ko3$sampleType == "Root",]
shapiro.test(residuals(lm(nsums~TP0, root)))

x <- root$TP0
y <- root$nsums
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

leaf <- ko3[ko3$sampleType == "Leaf",]
shapiro.test(residuals(lm(nsums~TP0, leaf)))

x <- leaf$TP0
y <- leaf$nsums
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

air <- ko3[ko3$sampleType == "Air",]
shapiro.test(residuals(lm(nsums~TP0, air)))

x <- air$TP0
y <- air$nsums
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")


p1 <- ggplot(data=ko3) + 
  geom_jitter(aes(x= TP0, y = nsums, color = sampleType, shape = sampleType), 
              height = 0,alpha = 0.4, size =2.1) + 
  geom_smooth(aes(x= TP0, y= nsums, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Nitrogen fixation")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,face="bold",color = "black"),
        axis.title=element_text(size=15,face="bold")) +
  annotate("text", x=1, y=.005, label="Air: R = -0.262, P < 0.001", 
           size=5, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1, y=.0045, label="Leaf: R = 0.153, P = 0.060", 
           size=5, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1, y=.004, label="Root: R = 0.589, P < 0.001", 
           size=5, hjust = 0, color="darkred",fontface="bold") 


p1

#setwd("~/Documents/sunR/Bacteria882/results2")
#ggsave("B-functions.fix.nif.no.naddition.pdf",p1,height = 5, width = 5)



###data.r.p####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
data <- read_xlsx("function_proportion_n.xlsx", sheet=2)


data_lm <- data %>% filter(level2=="Cell motility")
shapiro.test(residuals(lm(value~TP0, data_lm)))
hist(residuals(lm(value~TP0, data_lm)))

lm_rlt_Spcs <- function(df, yval, xval,ypos, subgrp = c("all")) {
  ypos_val <- ypos
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(level2 == subgrp)# change
    
  } else {
    
    df_tmp <- df
    
  }
  
  lmR_spearman <- cor.test(df_tmp[[yval]], df_tmp[[xval]], method = "spearman")
  lmR_spe_rho <- lmR_spearman$estimate %>% round(3)
  lmR_spe_pval <- lmR_spearman$p.value
  
  if(is.na(lmR_spe_pval)) {
    
    Rp_sig <- "NA"
    
  } else if(lmR_spe_pval < 0.001) {
    
    Rp_sig <- "***"
    
  } else if(lmR_spe_pval <= 0.01) {
    
    Rp_sig <- "**"
    
  } else if(lmR_spe_pval <= 0.05){
    
    Rp_sig <- "*"
    
  } else if(lmR_spe_pval <= 1){
    
    Rp_sig <- "NS"
    
  }
  
  data_sub <- data %>% filter(level2 == subgrp)
  
  
  df_rlt <- data.frame(
    level2 = subgrp,
    anno_x1 = (range(data_sub[[xval]])[2] - range(data_sub[[xval]])[1]) * 0.5 + range(data_sub[[xval]])[1],
    anno_y1 = (range(data_sub[[yval]], na.rm = T)[2] - range(data_sub[[yval]], na.rm = T)[1]) * ypos_val +
      range(data_sub[[yval]], na.rm = T)[1],
    rsig = str_c("R = ", lmR_spe_rho),
    psig = Rp_sig,
    anno_x2 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.6 + range(df_tmp[[xval]])[1],
    anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
      range(df_tmp[[yval]], na.rm = T)[1],
    pval = lmR_spe_pval
  )
  
  pval_sig <- str_c("P = ", lmR_spe_pval %>% round(3))
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, ", ", pval_sig))
  
  return(df_rlt_F)
  
}


data_air <- data[data$sampleType=="Air",]
data_leaf <- data[data$sampleType=="Leaf",]
data_root <- data[data$sampleType=="Root",]

level2_list <- data_air$level2 %>% unique()

data_air_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_air,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 1.1,
                                     subgrp = .))

data_air_subdata$sampleType <- "Air"


data_leaf_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_leaf,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos=1.0,
                                     subgrp = .))

data_leaf_subdata$sampleType <- "Leaf"


data_root_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_root,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 0.9,
                                     subgrp = .))

data_root_subdata$sampleType <- "Root"


data1_subdata <- rbind(data_air_subdata,data_leaf_subdata,data_root_subdata)
#write.csv(data1_subdata, "data_subdata_n.csv")



p_d <- ggplot(data=data) + 
  geom_jitter(aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =1.5) + 
  geom_text(data = data1_subdata,aes(x=anno_x1,y=anno_y1,label = rlt_sig,color=sampleType),parse = F, show.legend = F)+
  facet_wrap(~level2, scales = "free") +
  geom_smooth(aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Mean proportion(%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))+
  theme_bw()+
  guides(colour = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=5),nrow=1),
         shape = guide_legend(title = "Compartment",nrow=1))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold",color="black"),
        legend.title = element_text(colour="black", size=15, face="bold"),
        legend.text = element_text(colour="black", size=15, face="bold"),
        #legend.position = "none", 
        legend.position = "inside",
        legend.position.inside = c(0.8,0.1),
        #legend.box = "vertical",
        strip.text = element_text(size=10,face="bold"),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=17,face="bold",color="black")) 

p_d
#setwd("~/Documents/sunR/Bacteria882/results2")
#ggsave("B-functions.all.pdf", p_d, height = 15, width = 22)


##Figure9.Function.Nitrogen.genes####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)

rm(list=ls())

setwd("~/Documents/sunR/Bacteria882/new-picrust2")
data <- read_xlsx("function_proportion_n.xlsx", sheet=2)
data_subdata <- read_xlsx("data_subdata_n.xlsx",sheet=1)
subdata <- read_xlsx("data_subdata_n.xlsx",sheet=2)
functions <- subdata$level2
data1 <- data[data$level2 %in% functions, ]
data_subdata1 <- data_subdata[data_subdata$level2 %in% functions,]

data1$level2 <- factor(data1$level2, levels = c("Nitrogen fixation",
                                                "Amino acid metabolism","Transcription",
                                                #"Signaling molecules and interaction",
                                                "Energy metabolism",
                                                "Metabolism of terpenoids and polyketides",
                                                "Carbohydrate metabolism",
                                                #"Cell growth and death",
                                                "Nucleotide metabolism",
                                                "Xenobiotics biodegradation and metabolism",
                                                #"Folding, sorting and degradation",
                                                #"Translation",
                                                
                                                "Signal transduction",
                                                #"Information processing in viruses",
                                                "Glycan biosynthesis and metabolism",
                                                "Transport and catabolism",
                                                "Cell motility"))

data_subdata1$level2 <- factor(data_subdata1$level2, levels = c("Nitrogen fixation",
                                                                "Amino acid metabolism","Transcription",
                                                                #"Signaling molecules and interaction",
                                                                "Energy metabolism",
                                                                "Metabolism of terpenoids and polyketides",
                                                                "Carbohydrate metabolism",
                                                                #"Cell growth and death",
                                                                "Nucleotide metabolism",
                                                                "Xenobiotics biodegradation and metabolism",
                                                                #"Folding, sorting and degradation",
                                                                #"Translation",
                                                                
                                                                "Signal transduction",
                                                                #"Information processing in viruses",
                                                                "Glycan biosynthesis and metabolism",
                                                                "Transport and catabolism",
                                                                "Cell motility"))

p_d1 <- ggplot(data=data1) + 
  geom_jitter(aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =1.5) + 
  geom_text(data = data_subdata1,
            aes(x=2.0,y=anno_y1,label = rlt_sig,color=sampleType),
            parse = F, show.legend = F, size = 5,hjust = 0,
            fontface = "bold")+
  facet_wrap(~level2, scales = "free") +
  geom_smooth(aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Mean proportion (%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(colour = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=5),nrow=1),
         shape = guide_legend(title = "Compartment",nrow=1))+
  theme(plot.title = element_text(hjust = 0, size = 15, face = "bold",color="black"),
        legend.title = element_text(colour="black", size=17, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        #legend.position = "none", 
        legend.position = "top",
        legend.justification = c(1, 0),
        #legend.position.inside = c(0.87,0.1),
        #legend.box = "vertical",
        strip.text = element_text(size=11,face="bold"),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=20,face="bold",color="black")) 

p_d1

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("Figure4.pdf", p_d1, height = 12, width = 16)



##Supplementary.Figure1.Alpha.diversity####
library(vegan)
library(picante)
library(ggplot2)
library(car)
#library(ggcor)
library(ggpubr)
library(splitstackshape)
library(dplyr)
library(stringr)
library(broom)
library(agricolae)  
library(readxl)
library(patchwork)
##Bacteria-air7, root, leaf
rm(list = ls())
setwd("~/Documents/sunR/Bacteria882")

## 
group <- read_xlsx("sample_metadata.xlsx", sheet = 3)

group <- group[order(group$sample.id), ] 
group1 <- group[, c("sample.id","TP00","TP0","sampleType","sampleSite","Treatment")]
group2 <- group1 %>% filter(Treatment=="Control")
# 
#otu_all0 <- read.csv("otu_bactFlattening.csv", check.names = F, row.names = 1, header = T)
otu_all0 <- read.csv("otu_bactFlattening825.csv", check.names = F, row.names = 1, header = T)
otu_all <- otu_all0 %>% dplyr::select(all_of(group2$sample.id))

otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -702]

otu_all1$OTU.ID <- rownames(otu_all1)
otu_all1 <- otu_all1[, 1:701]  # 
otu_all1 <- t(otu_all1) #

#tree <- read.tree("amf-Newick Export.nwk")
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = "shannon", base = base)
  Simpson <- diversity(x, index = "simpson") 
  Pielou <- Shannon / log(ACE, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- "PD_whole_tree"
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}


alpha_all <- alpha(otu_all1, base = 2)
alpha_all$sample.id <- rownames(alpha_all)

##
sum(alpha_all$sample.id %in% group2$sample.id == FALSE) # 
which(alpha_all$sample.id %in% group2$sample.id == FALSE)

alpha_all0 <- merge(alpha_all, group2, by = "sample.id", all.x = TRUE) # 
#write.csv(alpha_all0, "Bacteria_alpha.csv")


####alpha.diversity.Shannon.temporal.no-N####
alpha_root <- alpha_all0[alpha_all0$sampleType == "Root", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_root)))

x <- alpha_root$TP0
y <- alpha_root$Shannon
cor.test(x, y, method = "spearman")


alpha_leaf <- alpha_all0[alpha_all0$sampleType == "Leaf", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_leaf)))

x <- alpha_leaf$TP0
y <- alpha_leaf$Shannon
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Shannon~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Shannon
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")


# anova <- aov(Shannon~TP0*sampleType, data = alpha_all0) 
# summary(anova)
library(lmerTest)
model <- lmer(Shannon ~ TP0*sampleType + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)


res <- resid(model, type = "pearson")
fit <- fitted(model)
par(mfrow=c(2,2))
plot(fit, res, xlab="Fitted", ylab="Pearson residuals"); abline(h=0, lty=2)
hist(res, breaks=30, main="Histogram of residuals")
qqnorm(res); qqline(res)                        
plot(abs(res) ~ fit, ylab="|residuals|", xlab="Fitted")  
par(mfrow=c(1,1))

shapiro.test(res)

acf(res, main="ACF of residuals")




p.sh <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Shannon, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Shannon, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Shannon index")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0.5, 10)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=2, label="Air: R = -0.334, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=1.4, label="Leaf: R = 0.194, P = 0.016", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=0.8, label="Root: R = 0.562, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") 


p.sh


p1 <- p.sh + labs(subtitle = "Time: F = 7.100, P = 0.008 \nCompartment: F = 187.389, P < 0.001 \nTime × Compartment: F = 22.794, P < 0.001")
p1 <- p1 + theme(plot.subtitle = element_text(size = 14,face = "bold",color="black")) 
p1

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.shannon.pdf",p1,height = 5, width = 5)

kw <- kruskal(alpha_all0$Shannon, alpha_all0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 430.0052, df = 2, p-value = 0
plot(kw)

#ggboxplot(genomesize0, x = "TP00", y = "value/1000000", color = "sampleType", add = "jitter")
#my_comparisons <- list(c("Air","Leaf"),c("Air","Root"),c("Leaf","Root"))
p_shannon <- ggplot(alpha_all0, aes(x = sampleType, y = Shannon)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "Shannon index") + 
  geom_text(x=2,y=2,label=expression(paste(chi^2, " = ", 430.005, ", ", "Df = ", 2, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=9.5,label="a",size=7,color="black")+
  annotate("text",x=2,y=8,label="c",size=7,color="black")+
  annotate("text",x=3,y=8,label="b",size=7,color="black")+
  ylim(0.5, 10)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_shannon

####alpha.diversity.Richness.temporal.no-N####
alpha_root <- alpha_all0[alpha_all0$sampleType == "Root", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_root)))

x <- alpha_root$TP0
y <- alpha_root$Richness
cor.test(x, y, method = "spearman")


alpha_leaf <- alpha_all0[alpha_all0$sampleType == "Leaf", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_leaf)))

x <- alpha_leaf$TP0
y <- alpha_leaf$Richness
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")

alpha_air <- alpha_all0[alpha_all0$sampleType == "Air", ]
shapiro.test(residuals(lm(Richness~TP0, alpha_air)))

x <- alpha_air$TP0
y <- alpha_air$Richness
cor.test(x, y, method = "spearman")
#cor.test(x, y, method = "pearson")


model <- lmer(Richness ~ TP0*sampleType + (1|sampleSite),
              data = alpha_all0)
summary(model)
anova(model)


res <- resid(model, type = "pearson")
fit <- fitted(model)
par(mfrow=c(2,2))
plot(fit, res, xlab="Fitted", ylab="Pearson residuals"); abline(h=0, lty=2)
hist(res, breaks=30, main="Histogram of residuals")
qqnorm(res); qqline(res)                         
plot(abs(res) ~ fit, ylab="|residuals|", xlab="Fitted") 
par(mfrow=c(1,1))

shapiro.test(res)

acf(res, main="ACF of residuals")




p.rich <- ggplot(data=alpha_all0) + 
  geom_jitter(aes(x= TP0, y = Richness, color = sampleType, shape = sampleType), 
              alpha = 0.4, size =2.1,height = 0) + 
  geom_smooth(aes(x= TP0, y= Richness, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2.6, show.legend = F)+
  scale_shape_manual(name="sampleType",values = c(15, 16, 17),guide=guide_legend(order=2)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "OTU richness")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  ylim(0,770)+
  theme_bw()+
  guides(color = guide_legend(title = "Compartment", override.aes = list(alpha = 1)),
         shape = guide_legend(title = "Compartment"))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        # legend.title = element_text(colour="black", size=12, face="bold"),
        # legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text=element_text(size=14,color = "black"),
        axis.title=element_text(size=19,face="bold",color = "black")) +
  annotate("text", x=1.5, y=120, label="Air: R = -0.396, P < 0.001", 
           size=6, hjust = 0, color="blue",fontface="bold") +
  annotate("text", x=1.5, y=70, label="Leaf: R = 0.285, P < 0.001", 
           size=6, hjust = 0, color="forestgreen",fontface="bold") +
  annotate("text", x=1.5, y=20, label="Root: R = 0.697, P < 0.001", 
           size=6, hjust = 0, color="darkred",fontface="bold") 


p.rich


p2 <- p.rich + labs(subtitle = "Time: F = 10.658, P = 0.001 \nCompartment: F = 646.375, P < 0.001 \nTime × Compartment: F = 65.233, P < 0.001")
p2 <- p2 + theme(plot.subtitle = element_text(size = 14,face = "bold",color="black")) 
p2

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.richness.pdf",p2,height = 5, width = 5)


kw <- kruskal(alpha_all0$Richness, alpha_all0$sampleType, p.adj = "bonferroni")
kw
#chi-squared = 524.9096, df = 2, p-value = 0
plot(kw)

#ggboxplot(genomesize0, x = "TP00", y = "value/1000000", color = "sampleType", add = "jitter")
#my_comparisons <- list(c("Air","Leaf"),c("Air","Root"),c("Leaf","Root"))
p_rich <- ggplot(alpha_all0, aes(x = sampleType, y = Richness)) + 
  geom_boxplot(aes(fill = sampleType),alpha = 0.5,outlier.colour = NA, 
               width = 0.4,
               position = position_dodge(width = 0.5)) +  
  geom_jitter(aes(color = sampleType),alpha = 0.5, 
              position = position_jitterdodge(dodge.width = 0.5))+
  scale_colour_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  scale_fill_manual(values= c("blue", "darkgreen","darkred"), name="Compartment") +
  labs(x = "Compartment",y = "OTU richness") + 
  geom_text(x=2,y=20,label=expression(paste(chi^2, " = ", 524.910, ", ", "Df = ", 2, ", ", "P = ", 0, sep="")),
            size=7,color="red",fontface="bold")+
  annotate("text",x=1,y=760,label="a",size=7,color="black")+
  annotate("text",x=2,y=470,label="c",size=7,color="black")+
  annotate("text",x=3,y=470,label="b",size=7,color="black")+
  ylim(0, 770)+
  theme_bw() +
  #theme(legend.position = c(0.9,0.83))
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.title = element_text(colour="black", size=12, face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "none", 
        axis.text.y = element_text(size=14,colour="black"),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title=element_text(size=19,face="bold",colour="black"))

p_rich


p <- p1+p_shannon+p2+p_rich
p

#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-alpha.shannon.richness.non.01.revised.pdf",p,height = 10, width = 10)



##Supplementary.Figure.titan.gs.rrn_n####
library(TITAN2)
library(ggplot2)
library(splitstackshape)
library(ggh4x) #facet_grid axis independent
library(colorRamps)
library(RColorBrewer)
library(ggsci)
library(patchwork)
library(readxl)
library(writexl)
rm(list = ls())
###root###
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-root-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Root"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
#otu_genome <- otu_genome[, -1]
colnames(otu_rrna)[1] <- "sample.id"
#otu_rrna <- otu_rrna[, -1]
colnames(otu_rrna)[2] <- "rrna"

library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")


genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]
rrna_up <- mean(genome_up$rrna)
genome_up0 <- na.omit(genome_up)
gs_up <- mean(genome_up0$genome_size/1000000)
gc_up <- mean(genome_up0$gc_percentage)
genome_down <- genome[genome$Titan == "z-", ]
rrna_down <- mean(genome_down$rrna)
genome_down0 <- na.omit(genome_down)
gs_down <- mean(genome_down0$genome_size/1000000)
gc_down <- mean(genome_down0$gc_percentage)



color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)


p1 <- ggplot()+
  geom_pointrange(data=dat3,
                  aes(x=order.x, y=zenv.cp, ymin=dat3$`5%`,ymax = dat3$`95%`,
                      shape=Titan,linetype=Titan,
                      color=factor(phylum.x)), size = 0.5)+
  scale_color_manual(values=color) +
  scale_shape_manual(values=c(16, 1)) + coord_flip() + ylab("Time") + xlab("Root - OTUs") + 
  scale_size_continuous(range=c(1,2))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(#size = guide_legend(title= "Abundance", nrow = 1), 
    color = guide_legend(title= "Phylum",order=2),
    shape = guide_legend(title= "Titan",order=1, nrow = 1),
    linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  #theme_set(theme_minimal() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))) +
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(colour="black", size=10, face="bold"),
    legend.position = "none",
    axis.text = element_text(size = 10,face="bold",color="black"),
    axis.title=element_text(size=16,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"))

p1


p2 <-   ggplot(dat3, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$genome_size)/1000000, 
               yend = mean(genome_down$genome_size)/1000000,
               x = 0, xend = 55, colour = "red") +
  geom_segment(y = mean(genome_up$genome_size)/1000000, 
               yend = mean(genome_up$genome_size)/1000000,
               x = 56, xend = 130, colour = "blue") +
  geom_point(data=dat3,aes(x=order.x, y=genome_size/1000000,shape=Titan,
                           color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p3 <-   ggplot(dat3, aes(order.x, gc_percentage)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = gc_percentage), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$gc_percentage), 
               yend = mean(genome_down$gc_percentage),
               x = 0, xend = 55, colour = "red") +
  geom_segment(y = mean(genome_up$gc_percentage), 
               yend = mean(genome_up$gc_percentage),
               x = 56, xend = 130, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = gc_percentage, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("GC percentage (%)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        #axis.line.x.top = element_line(color = "white"),
        #axis.line.y.right = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p3


p4 <-   ggplot(dat3, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$rrna), yend = mean(genome_down$rrna),
               x = 0, xend = 55, colour = "red") +
  geom_segment(y = mean(genome_up$rrna), yend = mean(genome_up$rrna),
               x = 56, xend = 130, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = rrna, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("16S copy number") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p4

p <- p1+p2+p4
p


###leaf###
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-leaf-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Leaf"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
colnames(otu_rrna)[1] <- "sample.id"
colnames(otu_rrna)[2] <- "rrna"


library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")


genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]
rrna_up <- mean(genome_up$rrna)
genome_up0 <- na.omit(genome_up)
gs_up <- mean(genome_up0$genome_size/1000000)
gc_up <- mean(genome_up0$gc_percentage)
genome_down <- genome[genome$Titan == "z-", ]
rrna_down <- mean(genome_down$rrna)
genome_down0 <- na.omit(genome_down)
gs_down <- mean(genome_down0$genome_size/1000000)
gc_down <- mean(genome_down0$gc_percentage)


color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)



p1 <- ggplot()+
  geom_pointrange(data=dat3,
                  aes(x=order.x, y=zenv.cp, ymin=dat3$`5%`,ymax = dat3$`95%`,
                      shape=Titan,linetype=Titan,
                      color=factor(phylum.x)), size = 0.5)+
  scale_shape_manual(values=c(16, 1)) + coord_flip() + ylab("Time") + xlab("Leaf - OTUs") + 
  scale_size_continuous(range=c(1,2))+
  scale_color_manual(values=color) +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  guides(#size = guide_legend(title= "Abundance", nrow = 1), 
    color = guide_legend(title= "Phylum",order=2),
    shape = guide_legend(title= "Titan",order=1, nrow = 1),
    linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(colour="black", size=10, face="bold"),
    axis.text = element_text(size = 10,face="bold",color="black"),
    axis.title=element_text(size=16,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"),
    legend.position = "none")

p1


p2 <-   ggplot(dat3, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$genome_size)/1000000, 
               yend = mean(genome_down$genome_size)/1000000,
               x = 0, xend = 29, colour = "red") +
  geom_segment(y = mean(genome_up$genome_size)/1000000, 
               yend = mean(genome_up$genome_size)/1000000,
               x = 30, xend = 47, colour = "blue") +
  geom_point(data=dat3,aes(x=order.x, y=genome_size/1000000,shape=Titan,
                           color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p3 <-   ggplot(dat3, aes(order.x, gc_percentage)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = gc_percentage), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$gc_percentage), 
               yend = mean(genome_down$gc_percentage),
               x = 0, xend = 29, colour = "red") +
  geom_segment(y = mean(genome_up$gc_percentage), 
               yend = mean(genome_up$gc_percentage),
               x = 30, xend = 47, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = gc_percentage, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("GC percentage (%)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p3


p4 <-   ggplot(dat3, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = mean(genome_down$rrna), yend = mean(genome_down$rrna),
               x = 0, xend = 29, colour = "red") +
  geom_segment(y = mean(genome_up$rrna), yend = mean(genome_up$rrna),
               x = 30, xend = 47, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = rrna, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("16S copy number") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p4



p <- p1+p2+p4
p


###air###
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-air-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan" 
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Air"

##genome data
otu_genome <- read_xlsx("../../otu_genometraits.xlsx",sheet=1)
otu_rrna <- read_xlsx("../../otu_genometraits.xlsx",sheet=2)
colnames(otu_genome)[2] <- "genome_size"
colnames(otu_genome)[1] <- "sample.id"
colnames(otu_rrna)[1] <- "sample.id"
colnames(otu_rrna)[2] <- "rrna"


library(dplyr)
dat2 <- left_join(dat1, otu_genome, by = "sample.id")
dat3 <- left_join(dat2, otu_rrna, by = "sample.id")

genome <- dat3[, c("sample.id", "Titan", "genome_size", "gc_percentage", "rrna")]
#write.csv(genome, "genome_average.csv")
genome_up <- genome[genome$Titan == "z+", ]

gs_up <- mean(genome_up$genome_size/1000000, na.rm=T)
gc_up <- mean(genome_up0$gc_percentage, na.rm=T)
rrna_up <- mean(genome_up$rrna, na.rm=T)

genome_down <- genome[genome$Titan == "z-", ]

gs_down <- mean(genome_down0$genome_size/1000000, na.rm=T)
gc_down <- mean(genome_down0$gc_percentage, na.rm=T)
rrna_down <- mean(genome_down$rrna, na.rm=T)



color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)


p1 <- ggplot()+
  geom_pointrange(data=dat3,
                  aes(x=order.x, y=zenv.cp, ymin=dat3$`5%`,ymax = dat3$`95%`,
                      shape=Titan,linetype=Titan,
                      color=factor(phylum.x)), size = 0.5)+
  scale_shape_manual(values=c(16, 1)) + coord_flip() + ylab("Time") + xlab("Air - OTUs") + 
  scale_size_continuous(range=c(1,2))+
  scale_color_manual(values=color) +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  # facet_wrap(~Habitat, nrow = 1,strip.position= 'top',scales ="free_y")+
  theme_bw()+
  guides(#size = guide_legend(title= "Abundance", nrow = 1), 
    color = guide_legend(title= "Phylum",order=2),
    shape = guide_legend(title= "Titan",order=1, nrow = 1),
    linetype = guide_legend(title= "Titan",order=1))+
  theme(panel.background = element_rect(fill = NA)) + 
  theme(#strip.text = element_text(size = 15,face="bold"), #facet labels
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(colour="black", size=10, face="bold"),
    axis.text = element_text(size = 10,face="bold",color="black"),
    axis.title=element_text(size=16,face="bold",color="black"),
    axis.line = element_line(linetype = "solid"),
    legend.position = "none")

p1



p2 <-   ggplot(dat3, aes(x=order.x, y=genome_size/1000000)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = genome_size/1000000), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = gs_down, 
               yend = gs_down,
               x = 0, xend = 198, colour = "red") +
  geom_segment(y = gs_up, 
               yend = gs_up,
               x = 199, xend = 325, colour = "blue") +
  geom_point(data=dat3,aes(x=order.x, y=genome_size/1000000,shape=Titan,
                           color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("Genome Size (Mb)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p2

p3 <-   ggplot(dat3, aes(order.x, gc_percentage)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = gc_percentage), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = gc_down, yend = gc_down,
               x = 0, xend = 198, colour = "red") +
  geom_segment(y = gc_up, yend = gc_up,
               x = 199, xend = 325, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = gc_percentage, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("GC percentage (%)") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p3


p4 <-   ggplot(dat3, aes(order.x, rrna)) +
  geom_segment(aes(x = order.x, xend = order.x, y = 0, yend = rrna), 
               linewidth = 0.5, color = "gray90") +
  geom_segment(y = 3.681, yend = 3.681,
               x = 0, xend = 198, colour = "red") +
  geom_segment(y = rrna_up, yend = rrna_up,
               x = 199, xend = 325, colour = "blue") +
  geom_point(data = dat3,aes(x = order.x, y = rrna, shape=Titan,
                             color=factor(phylum.x)), size = 2)+
  theme_bw() + coord_flip() + 
  scale_shape_manual(values=c(16, 1)) + ylab("16S copy number") + xlab("") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=15,face="bold",color="black"),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.line.x = element_line(linetype = "solid"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        legend.position = "none") +
  scale_color_manual(values=color)


p4

p <- p1+p2+p4
#library(aplot)
#p <- p1 %>% insert_right(p2) %>% insert_right(p4) 
p
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-titan.gs.rrn.air_n.pdf",p,height = 6, width = 8)



##TITAN_abundance####
library(TITAN2)
library(ggplot2)
library(splitstackshape)
library(ggh4x) #facet_grid axis independent
library(colorRamps)
library(RColorBrewer)
library(ggsci)
library(patchwork)
library(readxl)
library(writexl)
##otu.relative.abundance####
setwd("~/Documents/sunR/Bacteria882")
#rm(list = ls())
##1.
group0 <- read_xlsx("sample_metadata.xlsx", sheet = 3)
group <- group0[order(group0$sample.id),] 
group1 <- group[,c("sample.id","sampleType","TP00","Treatment")]
group2 <- na.omit(group1)

ID <- read_xlsx("ID_bacteria.xlsx", sheet = 1)
otu0 <- read.csv("otu_bactFlattening825.csv", head = T, row.names = 1)

otu_all <- otu0 %>% select(colnames(otu0)[colnames(otu0) %in% group2$sample.id])
otu_all$sum <- rowSums(otu_all)
zero <- which(otu_all$sum == 0)
otu_all1 <- otu_all[-zero, ]
otu_all1 <- otu_all1[, -826]

otu_rowname <- rownames(otu_all1)
rownames(ID) <- ID$names
ID <- ID[otu_rowname,]

ID[ID == "NA"] <- NA

##7.otu
ID$otu <- paste(ID$names, sep = ".")

bact.otu <- aggregate(otu_all1,by=list(ID$otu) , sum) 
rownames(bact.otu) <- bact.otu[,1]; 
bact.otu <- bact.otu[,-1] 
data6 <- data.frame(t(bact.otu))

total6 <- apply(data6, 1, sum); 
bact.relabu.otu <- data.frame(lapply(data6, function(x) {  x / total6  }) )

bact.otu0 <- bact.relabu.otu

bact.otu0 <- bact.otu0[,order(-colSums(bact.otu0))] ## 
lev <- interaction(group2$sampleType,group2$Treatment,sep = ":") ## combining factors for Barplot profiling
bact.otu.lev <- aggregate(bact.otu0, by = list(lev) , FUN = "mean") # generate the mean for each factor level

rel <- bact.otu.lev
rel1 <- setNames(as.data.frame(t(rel[,-1])), rel$Group.1)

#rm(list=ls())
setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
##air####
load("res-bact-air-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Air"


rel1_filtered <- rel1[rownames(rel1) %in% dat1$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)
rel2 <- rel1_filtered %>% 
  select("sample.id", "Air:Control") %>% 
  rename(abundance = "Air:Control")

dat1.1 <- merge(dat1,rel2,by="sample.id")


#setwd("~/Documents/sunR/Bacteria882/titan2")
##root####
load("res-bact-root-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat2 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat2$compartment = "Root"


rel1_filtered <- rel1[rownames(rel1) %in% dat2$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)
rel2 <- rel1_filtered %>% 
  select("sample.id", "Root:Control") %>% 
  rename(abundance = "Root:Control")

dat2.1 <- merge(dat2,rel2,by="sample.id")

#setwd("~/Documents/sunR/Bacteria882/titan2")
##leaf####
load("res-bact-leaf-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat3 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat3$compartment = "Leaf"

#
rel1_filtered <- rel1[rownames(rel1) %in% dat3$sample.id, ]
rel1_filtered$sample.id <- rownames(rel1_filtered)
rel2 <- rel1_filtered %>% 
  select("sample.id", "Leaf:Control") %>% 
  rename(abundance = "Leaf:Control")

dat3.1 <- merge(dat3,rel2,by="sample.id")


dat <- data.frame(rbind(dat1.1, dat2.1, dat3.1))
unique(dat$phylum)

color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)



p1 <- ggplot()+
  geom_linerange(data=dat,
                 aes(x=order.x, y=zenv.cp, ymin = dat$X5., ymax = dat$X95.,shape=Titan, 
                     linetype=Titan,color=factor(phylum)), lwd = 0.5,size =0.5)+
  geom_point(data = dat,
             aes(x = order.x, y = zenv.cp, shape = Titan,
                 color = factor(phylum), 
                 size = dat$abundance),  
             alpha = 0.8) +  
  scale_size_continuous(range = c(.1, 5)) +
  facet_grid2(~compartment,scales ="free",independent="all",)+
  scale_shape_manual(values=c(16, 1)) + 
  coord_flip() + 
  ylab("Time") + xlab("OTUs") + 
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  theme(panel.background = element_rect(fill = NA)) + 
  guides(size = guide_legend(title= "Abundance",nrow=1), 
         color=guide_legend(title= "Phylum",order=2,override.aes = list(size = 3)),
         shape=guide_legend(title= "Titan",order=1,nrow = 1,override.aes = list(size = 3)),
         linetype=guide_legend(title= "Titan",order=1))+
  theme(text = element_text(face = "bold"),
        strip.text = element_text(size = 15,face="bold",color="black"),
        legend.title = element_text(colour="black", size=15, face="bold"),
        legend.text = element_text(colour="black", size=13, face="bold"),
        plot.title = element_text(size=16, face="bold",hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        axis.line = element_line(linetype = "solid"))+
  scale_color_manual(values=color)

p1


#setwd("~/Documents/sunR/Bacteria882/n.treatment/ISMEC")
#ggsave("Bact-titan.abundance.pdf",p1,height = 6, width = 12)


##TITAN####
library(TITAN2)
library(ggplot2)
library(splitstackshape)
library(ggh4x) #facet_grid axis independent
library(colorRamps)
library(RColorBrewer)
library(ggsci)
library(patchwork)
library(readxl)
library(writexl)
##Bacteria
rm(list=ls())

setwd("~/Documents/sunR/Bacteria882/n.treatment/titan_n")
load("res-bact-air-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat1 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat1$compartment = "Air"

#setwd("~/Documents/sunR/Bacteria882/titan2")
load("res-bact-root-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat2 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
dat2$compartment = "Root"

#setwd("~/Documents/sunR/Bacteria882/titan2")
load("res-bact-leaf-titan_n.RData")

dat0 <- as.data.frame(res$sppmax)
names(dat0)[ncol(dat0)] <- "Titan"
dat0 <- dat0[which(dat0$Titan>0),]
dat0 <- dat0[order(dat0$Titan,decreasing = FALSE),]
dat0$order <- c(1:length(dat0[,"zenv.cp"]))
dat0$Titan[dat0$Titan=="1"] = "z-"
dat0$Titan[dat0$Titan=="2"] = "z+"
dat0$sample.id <- rownames(dat0)

dat3 <- merge(dat0, ID, by.x = "sample.id", by.y = "names", all.x = TRUE)
#dat2$primary_lifestyle[is.na(dat2$primary_lifestyle)] = "Others"
dat3$compartment = "Leaf"

dat <- data.frame(rbind(dat1, dat2, dat3))
unique(dat$phylum)

color <- c("Acidobacteriota" = "green", "Actinobacteriota" = "blue", "Armatimonadota" = "purple", 
           "Bacteroidota" = "black", "Chloroflexi" = "orange", "Crenarchaeota" = "lightblue", 
           "Cyanobacteria" = "red", "Deinococcota" = "darkcyan", "Desulfobacterota" = "darkred", 
           "Firmicutes" = "#00ffff", "Gemmatimonadota" = "#800080", "Myxococcota" = "darkgreen", 
           "Patescibacteria" = "hotpink", "Planctomycetota" = "deepskyblue", "Proteobacteria" = "maroon3", 
           "SAR324|clade(Marine|group|B)" = "tomato", "Verrucomicrobiota" = "navy","WPS-2"="yellow"
)

p1 <- ggplot()+
  geom_pointrange(data=dat,
                  aes(x=order.x, y=zenv.cp, ymin = dat$X5., ymax = dat$X95.,shape=Titan, 
                      linetype=Titan,color=factor(phylum)), lwd = 0.5,size =0.5)+
  facet_grid2(~compartment,scales ="free",independent="all",)+
  scale_shape_manual(values=c(16, 1)) + coord_flip() + ylab("Time") + xlab("OTUs") + 
  #ggtitle("Bacteria") +
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7),
                     labels=c("TP1","TP2","TP3","TP4","TP5","TP6","TP7"))+
  theme_bw()+
  theme(panel.background = element_rect(fill = NA)) + 
  guides(color=guide_legend(title= "Phylum",order=2),
         shape=guide_legend(title= "Titan",order=1),
         linetype=guide_legend(title= "Titan",order=1))+
  theme(text = element_text(face = "bold"),
        strip.text = element_text(size = 15,face="bold",color="black"),
        legend.title = element_text(colour="black", size=15, face="bold"),
        legend.text = element_text(colour="black", size=13, face="bold"),
        plot.title = element_text(size=16, face="bold",hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold",color="black"),
        axis.title=element_text(size=16,face="bold",color="black"),
        axis.line = element_line(linetype = "solid"))+
  scale_color_manual(values=color)

p1


#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("B-titan_n.pdf",p1,height = 6, width = 12)



##Supplementary.Figure8.Ncycle####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggh4x)
library(reshape2)
library(pals)
library(RColorBrewer)
library(psych)
library(linkET)
setwd("~/Documents/sunR/Bacteria882/new-picrust2")
#
group0 <- read_xlsx("../sample_metadata.xlsx", sheet = 3)

group2 <- group0[, c("sample.id","sampleType","TP00","TP0","Treatment")]
group2 <- na.omit(group2)
group3 <- group2 %>% filter(Treatment == "Control")


ko <- read_xlsx("ko_gene.xlsx", sheet = 1)
ko_gene <- read_xlsx("ko_gene.xlsx", sheet=2)
ko_ncycle <- read_xlsx("ko_gene.xlsx",sheet=3)

ko1 <- ko[,colnames(ko) %in% c("KO",group3$sample.id)]
ko1[, 2:702]  <- apply(ko1[, 2:702], 2, function(x) x / sum(x)) 

ko2 <- ko1[ko1$KO %in% ko_ncycle$KO,]
data1 <- merge(ko2,ko_ncycle,by="KO")
rownames(data1) <- data1$gene
data1 <- data1[,!(names(data1) %in% c("KO","gene","KO_name","EC"))]

data10 <- data.frame(t(data1))
data10$sample.id <- rownames(data10)
data11 <- merge(data10,group2,by="sample.id")

data12 <- gather(data11,key = "Ngene",value = "value", -"TP0",-"sample.id",-"sampleType",-"TP00",-"Treatment")

#data_subdata <- read_xlsx("Ngene.r.p.xlsx")
lm_rlt_Spcs <- function(df, yval, xval,ypos, subgrp = c("all")) {
  ypos_val <- ypos
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Ngene == subgrp)# change
    
  } else {
    
    df_tmp <- df
    
  }
  
  lmR_spearman <- cor.test(df_tmp[[yval]], df_tmp[[xval]], method = "spearman")
  lmR_spe_rho <- lmR_spearman$estimate %>% round(3)
  lmR_spe_pval <- lmR_spearman$p.value
  
  if(is.na(lmR_spe_pval)) {
    
    Rp_sig <- "NA"
    
  } else if(lmR_spe_pval < 0.001) {
    
    Rp_sig <- "***"
    
  } else if(lmR_spe_pval <= 0.01) {
    
    Rp_sig <- "**"
    
  } else if(lmR_spe_pval <= 0.05){
    
    Rp_sig <- "*"
    
  } else if(lmR_spe_pval <= 1){
    
    Rp_sig <- "NS"
    
  }
  
  data_sub <- data12 %>% filter(Ngene == subgrp)
  
  
  df_rlt <- data.frame(
    Ngene = subgrp,
    anno_x1 = (range(data_sub[[xval]])[2] - range(data_sub[[xval]])[1]) * 0.5 + range(data_sub[[xval]])[1],
    anno_y1 = (range(data_sub[[yval]], na.rm = T)[2] - range(data_sub[[yval]], na.rm = T)[1]) * ypos_val +
      range(data_sub[[yval]], na.rm = T)[1],
    rsig = str_c("R = ", lmR_spe_rho),
    psig = Rp_sig,
    anno_x2 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.6 + range(df_tmp[[xval]])[1],
    anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
      range(df_tmp[[yval]], na.rm = T)[1],
    pval = lmR_spe_pval
  )
  
  pval_sig <- str_c("P = ", lmR_spe_pval %>% round(3))
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, ", ", pval_sig))
  
  return(df_rlt_F)
  
}

data_air <- data12[data12$sampleType=="Air",]
data_leaf <- data12[data12$sampleType=="Leaf",]
data_root <- data12[data12$sampleType=="Root",]

level2_list <- data12$Ngene %>% unique()

data_air_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_air,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 1.1,
                                     subgrp = .))

data_air_subdata$sampleType <- "Air"

data_leaf_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_leaf,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos=1.0,
                                     subgrp = .))

data_leaf_subdata$sampleType <- "Leaf"

data_root_subdata <- 
  map_dfr(level2_list, ~ lm_rlt_Spcs(data_root,
                                     yval = "value",
                                     xval = "TP0",
                                     ypos = 0.9,
                                     subgrp = .))

data_root_subdata$sampleType <- "Root"

data_subdata1 <- rbind(data_air_subdata,data_leaf_subdata,data_root_subdata)
#write.csv(data_subdata1, "nitrogen.gene.r.p_n.csv")
data_subdata <- read_xlsx("nitrogen.gene.r.p_n.xlsx",sheet=1)




p_d <- ggplot(data=data12) + 
  geom_jitter(aes(x= TP0, y = value, color = sampleType, shape = sampleType), 
              alpha = 0.25, size =1.5) + 
  geom_text(data = data_subdata,
            aes(x=1.5,y=anno_y1,label = rlt_sig,color=sampleType),
            parse = F, show.legend = F, size = 5,hjust = 0,fontface = "bold")+
  facet_wrap(~Ngene, scales = "free") +
  geom_smooth(aes(x= TP0, y= value, group = sampleType, color = sampleType), 
              method ="lm", linewidth = 2, show.legend = F) +
  scale_shape_manual(name="sampleType",values = c(15, 16, 17)) +
  scale_colour_manual(values= c("blue", "forestgreen","darkred"), name="Time")+
  labs(x = "Time",y = "Mean proportion (%)")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))+
  theme_bw()+
  guides(colour = guide_legend(title = "Compartment", override.aes = list(alpha = 1,size=5),nrow=1),
         shape = guide_legend(title = "Compartment",nrow=1))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold",color="black"),
        legend.title = element_text(colour="black", size=20, face="bold"),
        legend.text = element_text(colour="black", size=17, face="bold"),
        #legend.position = "none", 
        legend.position = "inside",
        legend.position.inside = c(0.7,0.08),
        #legend.box = "vertical",
        strip.text = element_text(size=16,face="bold"),
        axis.text=element_text(size=10, color="black"),
        axis.title=element_text(size=25,face="bold",color="black")) 

p_d
#setwd("~/Documents/sunR/Bacteria882/n.treatment")
#ggsave("S.Figure.ncycle.pdf", p_d, height = 20, width = 22)


