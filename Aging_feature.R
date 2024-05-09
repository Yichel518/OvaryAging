setwd("F:/A_aging_lnc/")
library(readxl)
Apoptosis <- read_excel("F:/A_aging_lnc/Apoptosis.xlsx")
Autophagy <- read_excel("F:/A_aging_lnc/Autophagy.xlsx")
ROS <- read_excel("F:/A_aging_lnc/ROS.xlsx")
SASP <- read_excel("F:/A_aging_lnc/SASP.xlsx")
DNArepair <- read_excel("F:/A_aging_lnc/DNA repair.xlsx")
Fibrosis <- read_excel("F:/A_aging_lnc/Fibrosis.xlsx")
Lipid_biosynthesis <- read_excel("F:/A_aging_lnc/Lipid biosynthesis.xlsx")
Senescence <- read_excel("F:/A_aging_lnc/Senescence.xlsx")

UPR <- read_excel("F:/A_aging_lnc/UPR.xlsx")
UPR$`UPR pathway -related gene sets`

MtoH <- readRDS("./MtoH.rds")
Apoptosis <- MtoH[MtoH$HGNC.symbol%in%Apoptosis$`Apoptosis-related gene set`,]$MGI.symbol
Autophagy <- MtoH[MtoH$HGNC.symbol%in%Autophagy$`Autophagy-related gene set`,]$MGI.symbol
ROS <- MtoH[MtoH$HGNC.symbol%in%ROS$`ROS-related gene set`,]$MGI.symbol
SASP <- MtoH[MtoH$HGNC.symbol%in%SASP$`SASP-related gene set`,]$MGI.symbol
DNArepair <- MtoH[MtoH$HGNC.symbol%in%DNArepair$`DNA repair-related gene set`,]$MGI.symbol
Fibrosis <- MtoH[MtoH$HGNC.symbol%in%Fibrosis$`Fibrosis-related gene set`,]$MGI.symbol
Lipid_biosynthesis <- MtoH[MtoH$HGNC.symbol%in%Lipid_biosynthesis$`Lipid biosynthesis related genes`,]$MGI.symbol
Senescence <- MtoH[MtoH$HGNC.symbol%in%Senescence$`Senescence-related gene set`,]$MGI.symbol
UPR <- MtoH[MtoH$HGNC.symbol%in%UPR$`UPR pathway -related gene sets`[2:93],]$MGI.symbol



remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}

colnames(tpm)
Nor_tpm <- tpm[,c(7,8,9,10,11,12,16,17,18)]
colnames(Nor_tpm)

library(GSVA)
sets=list(Apoptosis,Autophagy,ROS,SASP,DNArepair,Fibrosis,Lipid_biosynthesis,Senescence)
#sets=lapply(sets, function(x) x[!is.na(x)])
sets[1]

#########calculation
#基因数据在exprMatrix中，把自己的数据赋值代入即可
exprMatrix=Nor_tpm
head(exprMatrix)
dim(exprMatrix)
gsva_matrix<- gsva(as.matrix(exprMatrix), sets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix1<- t(scale(t(gsva_matrix)))
head(gsva_matrix1)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1) 
head(nor_gsva_matrix1)
rownames(nor_gsva_matrix1) <- c("Apoptosis","Autophagy","ROS","SASP",
                                "DNArepair","Fibrosis","Lipid_biosynthesis","Senescence")

ggdata <- reshape2::melt(nor_gsva_matrix1)
ggdata$group <- substr(ggdata$Var2,1,6)

d <- ggdata %>% 
  dplyr::group_by(group) %>%
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(d, aes(value)) +
  geom_density(aes(fill = factor(group)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+facet_wrap(~Var1,scales ="free_y")


d$group <- factor(ggdata$group,levels = c("PRD_6M","PRD_9M","BRD_9M"))
ggplot(ggdata, aes(x = group, y =value,fill=group)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.8)+
  geom_boxplot(lwd=0.8,fatten = 0.8,width=0.6,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("PRD_6M","PRD_9M"),c("PRD_9M","BRD_9M")),
    map_signif_level = F,step_increase = 0.1,size = 0.8,textsize = 5,
    tip_length = 0,vjust = 0.1,test="t.test")+ylab("Gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none',plot.title = element_text(size = 16),strip.text.x =element_text(size = 16),
        strip.background = element_blank())+
  ggtitle("")+facet_wrap(~Var1,scales ="free_y",ncol = 4)
colnames(tpm3)
data <- tpm3[,c(1:3,7:12,16:18)]
n=t(scale(t(data)))
n <- as.data.frame(n)
n$SYMBOL <- rownames(n)
ggdata <- melt(n)
head(ggdata)
ggdata$group <- paste0(strsplit2(ggdata$variable,"_")[,1],"_",strsplit2(ggdata$variable,"_")[,2])
ggdata$group <- factor(ggdata$group,levels = c("PRD_6M","BRD_6M","PRD_9M","BRD_9M"))

########################## SASP #######################################
data1 <- ggdata[ggdata$SYMBOL%in%SASP,]
data1 <- data1 %>% 
   group_by(group) %>% 
   mutate(value = remove_outliers(value)) %>% 
   ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("SASP")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())


########################## Apoptosis #######################################
data1 <- ggdata[ggdata$SYMBOL%in%Apoptosis,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Apoptosis")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())


########################## Autophagy #######################################
data1 <- ggdata[ggdata$SYMBOL%in%Autophagy,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Autophagy")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## DNArepair #######################################
data1 <- ggdata[ggdata$SYMBOL%in%DNArepair,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("DNA repair")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## Fibrosis #######################################
data1 <- ggdata[ggdata$SYMBOL%in%Fibrosis,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Fibrosis")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## Lipid_biosynthesis #######################################
data1 <- ggdata[ggdata$SYMBOL%in%Lipid_biosynthesis,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Lipid storage")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## ROS #######################################
data1 <- ggdata[ggdata$SYMBOL%in%ROS,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("ROS")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## Senescence #######################################
data1 <- ggdata[ggdata$SYMBOL%in%Senescence,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Senescence")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

########################## UPR #######################################
data1 <- ggdata[ggdata$SYMBOL%in%UPR,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("UPR")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())


################################# C2 ##############################
data1 <- ggdata[ggdata$SYMBOL%in%c2$SYMBOL,]
ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())+ggtitle("Cluster1")

################################## C5 ###########################################
data1 <- ggdata[ggdata$SYMBOL%in%c5$SYMBOL,]
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()

ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())+ggtitle("Cluster2")

data1 <- ggdata[ggdata$SYMBOL%in%res[res$type=="up",]$SYMBOL,]

#去除outlier
#library(dplyr)
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()


ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("WT-upregulated Gene Set")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())
data1 <- ggdata[ggdata$SYMBOL%in%res[res$type=="down",]$SYMBOL,]
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
#去除outlier
#library(dplyr)
data1 <- data1 %>% 
  group_by(group) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()


ggplot(na.omit(data1),aes(x = group,y = log2(value+1),fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD_6M","PRD_6M"),c("BRD_9M","PRD_9M"),c("PRD_6M","PRD_9M")),
              map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("POF-upregulated Gene Set")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())
