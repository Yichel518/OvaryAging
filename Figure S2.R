setwd("F:/A_aging_lnc/")
count <- read.delim("F:/A_aging_lnc/mouse_age_newcount.txt", row.names=1, comment.char="#")
colnames(count) <- strsplit2(colnames(count),"[.]")[,7]
meta <- read.delim("F:/A_aging_lnc/metadata.txt")
colnames(count) <- meta$Sample
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)
gtf = rtracklayer::import("F:/A_MetaboAnalysis/Black_rice/weeks_22/Mus_musculus.GRCm39.104.gtf")
gtf = as.data.frame(gtf)
tra = gtf[gtf$type=="transcript",
          c("start","end","gene_id")]
glt = mutate(tra, length = end - start) %>%
  .[order(.$length),] %>% 
  filter(!duplicated(gene_id)) 
final_gene <- as.data.frame(glt[,3:4])
rownames(final_gene) <- final_gene$gene_id
final_gene <- final_gene[rownames(count),]
kb <- final_gene$length / 1000
rpk <- count/ kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.table(as.data.frame(tpm),"mm_age_tpm.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(as.data.frame(count),"mm_age_rawcount.txt", sep = '\t', col.names = NA, quote = FALSE)
ovary1 <- read.csv("F:/A_aging_lnc/GSM2906456_Ovary1_dge.txt/GSM2906456_Ovary1_dge.txt", sep="")
ovary2 <- read.csv("F:/A_aging_lnc/GSM2906457_Ovary2_dge.txt/GSM2906457_Ovary2_dge.txt", sep="")
dim(ovary1)
dim(ovary2)

rownames(ovary1)
dim(tpm)
data <- data.frame(name=rownames(ovary2))
write.table(as.data.frame(data),"name.txt", sep = '\t', col.names = NA, quote = FALSE)
mm_age_tpm$SYMBOL
mm_age_tpm$SYMBOL <- NULL
mm_age_tpm$ENTREZID <- NULL
mm_age_tpm$ENSEMBL <- NULL
mm_age_tpm$SYMBOL <- rownames(mm_age_tpm)
ggdata <- melt(mm_age_tpm)
ggdata$group <- strsplit2(ggdata$variable,"_")[,1]
ggdata$group <- factor(ggdata$group,levels = c("M3","M6","M9","M12"))
library(ggsci)
library(ggsignif)
#View(ggdata)
data <- ggdata[ggdata$SYMBOL=="Foxo3",]
ggplot(data,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3,lwd=1.2)+
  geom_boxplot(lwd=1.2,fatten = 1.2,width=0.6)+theme_classic()+
scale_fill_npg()+
  geom_signif(comparisons = list(c("M3", "M6"),c("M6","M9"),c("M9","M12")),
              map_signif_level = T,step_increase = 0.1,size = 1.2,textsize = 7,
              tip_length = 0,vjust = 0.2)+ylab("Relative expression (TPM)")+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+	
  theme(axis.title =element_text(size = 24,face = "bold"),axis.text =element_text(size = 24,face = "bold", color = 'black'))+	
  theme(axis.text = element_text(size = 24,face = "bold"),axis.text.x = element_text(face = "bold"),	
        legend.text = element_text(size=12,face = "bold"),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1.2),axis.line.y=element_line(size=1.2))+
  theme(legend.position = 'none')+ggtitle("Il17a")+ 
  theme(plot.title = element_text(size = 24, face = "bold.italic"))



mean = apply(mm_age_tpm,1,function(x){tapply(x,strsplit2(colnames(mm_age_tpm),"_")[,1],mean)}) %>% t() %>% as.data.frame()
eset <- new("ExpressionSet",exprs = as.matrix(mm_age_tpm))
# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset,min.std=0)
eset <- standardise(eset)

c <- 6
#  评估出最佳的m值
m <- mestimate(eset)
# 聚类
cl <- mfuzz(eset, c = c, m = m)
# 查看每个cluster中的基因个数
cl$size
# 提取某个cluster下的基因
cl$cluster[cl$cluster == 1]
# 查看基因和cluster之间的membership
cl$membership

library("ClusterGVis", lib.loc="C:/Users/tutu/AppData/Local/R/win-library/4.2")
colnames(mean)
ck <- clusterData(exp = mean[,c(2:4,1)],
                  cluster.method = "mfuzz",
                  cluster.num = 4)
visCluster(object = ck,
           plot.type = "line",
           textbox.size=14)

dd <- data.frame(row.names = rownames(WT_POF_Diff),`WT_vs_POF` =WT_POF_Diff$log2FoldChange)
pheatmap::pheatmap(dd,cluster_cols = F,cluster_rows = T,show_rownames = F,fontsize = 18,angle_col = 0)


tpm3$SYMBOL <- rownames(tpm3)
ggdata <- reshape2::melt(tpm3)
head(ggdata)
ggdata$group <- paste0(strsplit2(ggdata$variable,"_")[,1],"_",strsplit2(ggdata$variable,"_")[,2])
ggdata$group <- factor(ggdata$group,levels = c("PRD_6M","BRD_6M","PRD_9M","BRD_9M"))
#View(ggdata)
sym <- AnnotationDbi::select(org.Mm.eg.db,keys=mfuzz_cluster5_gene$V1,
                                   columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")

data1 <- ggdata[ggdata$SYMBOL%in%sym$SYMBOL,]
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Cluster2")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())

sym <-  AnnotationDbi::select(org.Mm.eg.db,keys=mfuzz_cluster2_gene$V1,
                                   columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")

data1 <- ggdata[ggdata$SYMBOL%in%sym$SYMBOL,]
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("Cluster1")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())




sym1 <-  AnnotationDbi::select(org.Mm.eg.db,keys=mfuzz_cluster2_gene$V1,
                              columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")

sym2 <-  AnnotationDbi::select(org.Mm.eg.db,keys=mfuzz_cluster5_gene$V1,
                               columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="ENSEMBL")

dim(mfuzz_cluster2_gene)

Aging_up <- intersect(sym1$SYMBOL,WT_POF_Diff[WT_POF_Diff$Class=="POF",]$SYMBOL)
data1 <- ggdata[ggdata$SYMBOL%in%POFup[POFup$log2FoldChange<-0.5,]$SYMBOL,]
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("POF-up")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())


data1 <- ggdata[ggdata$SYMBOL%in%WTup$SYMBOL,]
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("WT-up")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())
tpm <- tpm[,c(20:25)]
tpm$ENSEMBL <- rownames(tpm)
View(tpm)
colnames(count)
data <- count

coldata <- as.data.frame(c(rep("WT",3),rep("POF",3)))
rownames(coldata) <- colnames(data)
colnames(coldata) <- c("condition")
coldata
dds<- DESeqDataSetFromMatrix(countData =                                                                
                               data,
                             colData = coldata,
                             design = ~condition)
dds2<- DESeq(dds)
resultsNames(dds2)
res <- results(dds2,contrast=list(c("condition_WT_vs_POF")))
WT <- subset(res, padj < 0.05 & log2FoldChange > 0 )
POF <- subset(res, padj < 0.05 & log2FoldChange < 0 )
View(res)
res <- as.data.frame(res)
res$type <- ifelse(res$padj < 0.05,
                   ifelse(abs(res$log2FoldChange) >  0,
                          ifelse(res$log2FoldChange < 0 ,'down','up'),'noSig'),'noSig')


table(res$type)
res$ENSEMBL <- rownames(res)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(res),
                            columns=c("SYMBOL","ENSEMBL"),keytype="ENSEMBL")
res <- merge(as.data.frame(res),sym,all.x=T,by="ENSEMBL")
View(res)
meta <- res[res$SYMBOL%in%c("Cdkn1a","Cdkn1b","Cdkn2b","Zp3","Zp2","Zar1","Sycp2","Pcna"
),]
View(meta)
p=ggplot(na.omit(as.data.frame(res)),aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = type),size = 3,alpha=0.3) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#F39B7F','noSig'='grey','down'='#4BB5CE'),
                     # legend labels
                     label = c('up'='up (num=3456)','down'='down (num=5076)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('WT vs POF')+theme_classic()+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(plot.title = element_text(size = 24))
p
p + geom_text_repel(data = meta,aes(x = log2FoldChange,y = -log10(pvalue),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)


data1 <- ggdata[ggdata$SYMBOL%in%res[res$type=="up",]$SYMBOL,]
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
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
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("POF-upregulated Gene Set")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())
head(tpm3)
n <- t(scale(t(tpm3)))
n[n>2]=2
n[n<-2]=-2
mean_tpm3 = apply(n,1,function(x){tapply(x,Info$group,mean)}) %>% t() %>% as.data.frame()
