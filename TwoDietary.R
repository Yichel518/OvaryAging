setwd("~/scRNA_aging/output/")
all <- readRDS("./all_obj.rds")
library(Seurat)
meta <- all@meta.data[all$Diets%in%c("M6-BRD","M6-PRD","M9-BRD","M9-PRD"),]
meta$Diets <- factor(meta$Diets,levels = c("M6-PRD","M6-BRD","M9-PRD","M9-BRD"))
saveRDS(subset(all,cells=rownames(meta)),"BP_subset.rds")
rm(all)
gc()

BP <- readRDS("./BP_subset.rds")
BP$abbreviation2 <- gsub("SC","SC",BP$abbreviation)%>%gsub("TCC","SC",.)
BP$abbreviation2 <- gsub("STC","GC",BP$abbreviation2)
BP$abbreviation2[which(BP$seurat_clusters%in%c(5,16,20))]<- "EDC"

BP_meta <- BP@meta.data
unique(BP$abbreviation2)
BP <- RunUMAP(object = BP,dims = 1:40)
BP <- RunTSNE(object = BP)
BP$diet <- substr(BP$Diets,4,6)
DimPlot(object = BP, reduction = "tsne",label = F,label.size = 6,group.by = "abbreviation2")+scale_color_simpsons()
DimPlot(object = BP, reduction = "umap",label = F,label.size = 6,
        group.by = "abbreviation2")+scale_color_simpsons()+ggtitle("Total cell number: 94072")
length(colnames(BP))
DimPlot(object = BP, reduction = "umap",label = T,label.size = 6,
        group.by = "seurat_clusters")


DimPlot(object = BP, reduction = "umap",label = F,label.size = 6,
        group.by = "abbreviation2",split.by = c("age"))+scale_color_simpsons()

DimPlot(object = BP, reduction = "umap",label = F,label.size = 6,
        group.by = "abbreviation2",split.by = c("diet"))+scale_color_simpsons()


DotPlot(object = BP, features = c("Epcam","Krt19","Prom1","Aldh1a1"),group.by = "abbreviation2")


DotPlot(object = BP,group.by = "abbreviation2",
        cols = c("#4B549B","#E91C22"), 
        features = c("Amh","Cyp19a1","Inha","Cyp11a1",
                     "Col1a1","Col1a2",
                     "Epcam","Krt19","Cdh5","Pecam1","Vwf","Cebpb","Cxcr2","Csf3r","Itgam",
                     "Cd68","Clec9a","Itgax","Itgae","Cd19","Pax5","Cd79a",
                     "Nkg7","Tbx21","Ccl5","Cd4","Trac","Rorc"))+
  theme_classic()+
  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.title = element_text(size = 13) )+coord_flip()


saveRDS(BP,"./BP.rds")
BP$Diets
BP <- readRDS("./BP.rds")
BP_sub <- subset(BP,cells=BP@meta.data[BP$Diets%in%c("M6-PRD","M9-PRD","M9-BRD"),]%>%rownames(.))
rm(BP)
gc()

BP_sub$age <- substr(BP_sub$Diets,1,2)
n <- paste0(BP_sub$Diets,"-",BP_sub$abbreviation2)
names(n) <- colnames(BP_sub)
BP_sub@active.ident <-factor(n)
cor_data <- AverageExpression(BP_sub)
data <-filter_if(as.data.frame(cor_data$RNA), is.numeric, all_vars((.) != 0))
?cor
pheatmap(cor(data))






DimPlot(object = BP_sub, reduction = "umap",label = F,label.size = 6,
        group.by = "diet")+scale_color_npg()+ggtitle("Diet")

DimPlot(object = BP_sub, reduction = "umap",label = F,label.size = 6,
         group.by = "age",cols = c("#C08126","#4963B1"))+ggtitle("Age")



cellname <- paste0(BP_sub$Diets,"-",BP_sub$abbreviation2)
names(cellname) <- colnames(BP_sub)
BP_sub@active.ident <- factor(BP_sub$abbreviation2)

BP_sub@active.ident <- factor(cellname,levels = c("M6-PRD-GC","M9-PRD-GC","M9-BRD-GC",
                                                    "M6-PRD-SC","M9-PRD-SC","M9-BRD-SC",
                                                    "M6-PRD-EC","M9-PRD-EC","M9-BRD-EC",
                                                    "M6-PRD-EDC","M9-PRD-EDC","M9-BRD-EDC",
                                                    "M6-PRD-GLC","M9-PRD-GLC","M9-BRD-GLC",
                                                    "M6-PRD-M","M9-PRD-M","M9-BRD-M",
                                                    "M6-PRD-DC","M9-PRD-DC","M9-BRD-DC",
                                                    "M6-PRD-B","M9-PRD-B","M9-BRD-B",
                                                    "M6-PRD-NK","M9-PRD-NK","M9-BRD-NK",
                                                    "M6-PRD-T","M9-PRD-T","M9-BRD-T",
                                                    "M6-PRD-ILC","M9-PRD-ILC","M9-BRD-ILC"))

DotPlot(object = BP_sub,
        cols = c("#4B549B","#E91C22"),split.by = "abbreviation2",group.by = "Diets",
        features = unique(c("Amh","Cyp19a1","Inha","Cyp11a1",
                     "Col1a1","Col1a2",
                     "Epcam","Krt19","Cdh5","Pecam1","Vwf","Cebpb","Cxcr2","Csf3r","Itgam",
                     "Cd68","Clec9a","Itgax","Itgae","Cd19","Pax5","Cd79a",
                     "Nkg7","Tbx21","Ccl5","Cd4","Trac","Rorc","Cdkn1a")))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.title = element_text(size = 13) )+coord_flip()

?VlnPlot

VlnPlot(BP_sub,features = c("Cdkn1a","Cdkn2a","Cdkn2d"),split.by = "abbreviation2",group.by = "Diets")+scale_fill_simpsons()
BP_sub$Diets <- factor(BP_sub$Diets,levels = c("M6-PRD","M9-PRD","M9-BRD"))
FeaturePlot(BP_sub,features = c("Cdkn1a"),split.by = "Diets",cols = c("white","red3"),label = T)
VlnPlot(BP_sub,features = c("Cdkn1a"),group.by = "Diets",pt.size = 0)+scale_fill_d3()


BP_meta$sample <- strsplit2(rownames(BP_meta),"_")[,1]
data <- as.data.frame(table(BP_meta$sample))
data

dd <- strsplit2(data$Var1,"-")[,1:2]%>%as.data.frame()
colnames(dd) <- c("Stage","Diets")
dd$Diets <- paste0(dd$Diets,"RD")
data <- cbind(data,dd)
data$Diets <- factor(data$Diets,levels = c("BRD","PRD"))
ggplot(data, aes(fill=Diets, y=Freq, x=Stage))+
  geom_bar(position=position_dodge(),
           stat="summary",
           width=0.6,
           colour = "black",
           size=0.5)+ 
  stat_summary(fun.data = 'mean_se', 
               geom = "errorbar", 
               colour = "black",
               width = 0.2,
               position=position_dodge(0.6))+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_classic()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("#ffd401","#00b0eb"),name="")+
  labs(x="Stages",y="Cell number")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+
  annotate(geom = 'text', label="N=8196", x=0.9, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="N=8801", x=1.3, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="N=8592", x=1.9, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="N=8652", x=2.3, y=2000, angle=90,size=5)
?sc_utils
library("scProportionTest")
library(forcats, lib.loc = "/usr/local/lib/R/site-library")
library(dplyr)
library(IRanges)

BP$age <- substr(BP$Diets,1,2)
meta <- BP@meta.data[BP@meta.data$diet%in%c("PRD"),]
PRD <- subset(BP,cells=rownames(meta))
prop_test <- sc_utils(PRD)

prop_test <- permutation_test(
  prop_test, cluster_identity = "abbreviation2",
  sample_1 = "M6", sample_2 = "M9",
  sample_identity = "age"
)

permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(title  =element_text(size = 14),axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 15),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("GC","SC","EC", "EDC","GLC",
                            "M","DC","B","NK","T","ILC"))+
  ggtitle("PRD: 9m vs 6m")+xlab("Celltypes")+ylab("log2FD")


meta <- BP@meta.data[BP@meta.data$age%in%c("M9"),]
M9 <- subset(BP,cells=rownames(meta))
prop_test <- sc_utils(M9)
prop_test <- permutation_test(
  prop_test, cluster_identity = "abbreviation2",
  sample_1 = "BRD", sample_2 = "PRD",
  sample_identity = "diet"
)
permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(title  =element_text(size = 14),axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 15),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("GC","SC","EC", "EDC","GLC",
                            "M","DC","B","NK","T","ILC"))+
  ggtitle("9m: PRD vs BRD")+xlab("Celltypes")+ylab("log2FD")



permutation_plot2 <- function(
    sc_utils_obj,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    order_clusters = TRUE
) {
  
  ## Retrieve results.
  plot_data <- data.table::copy(sc_utils_obj@results$permutation)
  
  ## Mark the significant results.
  plot_data[, significance := ifelse(
    FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  )]
  
  plot_data[, significance := factor(significance, levels = c(
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  ))]
  
  ## Order the clusters by observed log2FD if requested.
  if (order_clusters) {
    plot_data[, clusters := fct_reorder(factor(clusters), dplyr::desc(obs_log2FD))]
  }
  
  ## Plot the results.
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
    theme_bw()  +
    geom_hline(yintercept = 0.58, lty = 2) +
    geom_hline(yintercept = -0.58, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("salmon", "grey")) +
    coord_flip()+theme(legend.position = "top")
  
  return(p)
}

BP$diet_cell <- paste0(BP$diet,"-",BP$abbreviation2)

cell <- c("BRD-GC","BRD-SC","BRD-EDC","BRD-EC",
          "BRD-M","BRD-DC","BRD-GLC","BRD-ILC","BRD-NK","BRD-B","BRD-T",
          "PRD-GC","PRD-SC","PRD-EDC","PRD-EC",
          "PRD-M","PRD-DC","PRD-GLC","PRD-ILC","PRD-NK","PRD-B","PRD-T")
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(BP,cells=rownames(BP@meta.data[BP$diet_cell==cell[x],]))
  saveRDS(tmp, paste0(cell[x],"_BP.rds"))
})  
  
  
  
cell <- c("GC","SC","EDC","EC","M","DC","GLC","ILC","NK","B","T")
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(BP,cells=rownames(BP@meta.data[BP$abbreviation2==cell[x],]))
  saveRDS(tmp, paste0(cell[x],"_BP.rds"))
})


meta <- BP@meta.data[BP@meta.data$age%in%c("M9"),]
M9 <- subset(BP,cells=rownames(meta))
saveRDS(M9,"./M9_BP.rds")
rm(M9)
gc()
meta <- BP@meta.data[BP@meta.data$age%in%c("M6"),]
M6 <- subset(BP,cells=rownames(meta))
saveRDS(M6,"./M6_BP.rds")
rm(M6)
gc()
rm(BP)
gc()
################################################################
########################## 差异基因 ##############################
########################## BRD ####################################
############################## GC####################
M9_BP <- readRDS("./M9_BP.rds")
BRD_GC <- readRDS("./BRD-GC_BP.rds")
BRD_GC@active.ident<- factor(BRD_GC$age)
BRD_GC_markers <- FindMarkers(BRD_GC, ident.1 = "M9", ident.2 = "M6", min.pct = 0.25)
BRD_GC_markers$group <- 0
BRD_GC_markers[BRD_GC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_GC_markers[BRD_GC_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_GC_markers$diets <- "BRD: 9m vs 6m"
data1 <- as.data.frame(table(BRD_GC_markers$group))
data1$diets <- "BRD: 9m vs 6m"
data1

PRD_GC <- readRDS("./PRD-GC_BP.rds")
PRD_GC@active.ident<-factor(PRD_GC$age)
PRD_GC_markers <- FindMarkers(PRD_GC, ident.1 = "M9", ident.2 = "M6", min.pct = 0.25)
PRD_GC_markers$group <- 0
PRD_GC_markers[PRD_GC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_GC_markers[PRD_GC_markers$avg_log2FC < -0.25,]$group <- "6M"
PRD_GC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_GC_markers$group))
data3$diets <- "PRD: 9m vs 6m"


df1=as.data.frame(M9_BP@active.ident)
idents <- paste0(M9_BP$diet,"_",M9_BP$abbreviation2)
names(idents) <- rownames(M9_BP@meta.data)
M9_BP@active.ident<-factor(idents)
unique(factor(idents))
PB_GC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_GC",
                             ident.2 = "BRD_GC", min.pct = 0.25)
PB_GC_markers$group <- 0
PB_GC_markers[PB_GC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_GC_markers[PB_GC_markers$avg_log2FC < -0.25,]$group <- "BRD"
PB_GC_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_GC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "GC"
d1 <- data
head(d1)
#Var1 Freq         diets celltype
#1    0 2295 BRD: 9m vs 6m       GC
#2   6M  804 BRD: 9m vs 6m       GC
#3   9M  652 BRD: 9m vs 6m       GC
#4    0 1957 PRD: 9m vs 6m       GC
#5   6M 1399 PRD: 9m vs 6m       GC
#6   9M 1506 PRD: 9m vs 6m       GC
GC_age_up <- intersect(PRD_GC_markers[PRD_GC_markers$group=="9M",]%>%rownames(.),
          PB_GC_markers[PB_GC_markers$group=="PRD",]%>%rownames(.))
length(GC_age_up) # 965
GC_age_down <-intersect(PRD_GC_markers[PRD_GC_markers$group=="6M",]%>%rownames(.),
          PB_GC_markers[PB_GC_markers$group=="BRD",]%>%rownames(.))
length(GC_age_down) # 657
length(PRD_GC_markers%>%rownames(.))
length(BRD_GC_markers%>%rownames(.))

length(BRD_GC_markers[abs(BRD_GC_markers$avg_log2FC > 0.25),]$group)
length(PRD_GC_markers[abs(PRD_GC_markers$avg_log2FC > 0.25),]$group)

########################## SC ####################################
rm(BRD_GC,PRD_GC)
gc()
BRD_SC <- readRDS("./BRD-SC_BP.rds")
df1=as.data.frame(BRD_SC@meta.data)
View(df1)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_SC@meta.data)
BRD_SC@active.ident<-factor(idents)
BRD_SC_markers <- FindMarkers(BRD_SC, ident.1 = "M9_SC", ident.2 = "M6_SC", min.pct = 0.25)
BRD_SC_markers$group <- 0
BRD_SC_markers[BRD_SC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_SC_markers[BRD_SC_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_SC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_SC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_SC <- readRDS("./PRD-SC_BP.rds")
df1=as.data.frame(PRD_SC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_SC@meta.data)
PRD_SC@active.ident<-factor(idents)
PRD_SC_markers <- FindMarkers(PRD_SC, ident.1 = "M9_SC", ident.2 = "M6_SC", min.pct = 0.25)
PRD_SC_markers$group <- 0
PRD_SC_markers[PRD_SC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_SC_markers[PRD_SC_markers$avg_log2FC < -0.25,]$group <- "6M"
PRD_SC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_SC_markers$group))
data3$diets <- "PRD: 9m vs 6m"

PB_SC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_SC",
                             ident.2 = "BRD_SC", min.pct = 0.25)
PB_SC_markers$group <- 0
PB_SC_markers[PB_SC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_SC_markers[PB_SC_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_SC_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_SC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "SC"
d5 <- data
d5
#Var1 Freq          diets celltype
#1    0 2094  BRD: 9m vs 6m       SC
#2   6M  982  BRD: 9m vs 6m       SC
#3   9M  402  BRD: 9m vs 6m       SC
#4    0 1878  PRD: 9m vs 6m       SC
#5   6M  811  PRD: 9m vs 6m       SC
#6   9M 1860  PRD: 9m vs 6m       SC
#7    0 1884 9m: PRD vs BRD       SC
#8  BRD  655 9m: PRD vs BRD       SC
#9  PRD 1974 9m: PRD vs BRD       SC
SC_age_up <- intersect(PRD_SC_markers[PRD_SC_markers$group=="9M",]%>%rownames(.),
                       PB_SC_markers[PB_SC_markers$group=="PRD",]%>%rownames(.))

SC_age_down <-intersect(PRD_SC_markers[PRD_SC_markers$group=="6M",]%>%rownames(.),
                        PB_SC_markers[PB_SC_markers$group=="BRD",]%>%rownames(.))


length(PRD_SC_markers[PRD_SC_markers$group=="9M",]%>%rownames(.))

length(SC_age_up) # 1422
length(SC_age_down) # 456
length(BRD_SC_markers[abs(BRD_SC_markers$avg_log2FC>0.25),]$group)
length(PRD_SC_markers[abs(PRD_SC_markers$avg_log2FC>0.25),]$group)
rm(BRD_SC,PRD_SC)
gc()
########################## B ####################################
BRD_B <- readRDS("./BRD-B_BP.rds")
df1=as.data.frame(BRD_B@active.ident)
idents <- gsub("B cells","B",BRD_B$group_age)
names(idents) <- rownames(BRD_B@meta.data)
BRD_B@active.ident<-factor(idents)
BRD_B_markers <- FindMarkers(BRD_B, ident.1 = "M9_B", ident.2 = "M6_B", min.pct = 0.25)
BRD_B_markers$group <- 0
BRD_B_markers[BRD_B_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_B_markers[BRD_B_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_B_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_B_markers$group))
data1$diets <- "BRD: 9m vs 6m"

PRD_B <- readRDS("./PRD-B_BP.rds")
df1=as.data.frame(PRD_B@active.ident)
idents <- gsub("B cells","B",PRD_B$group_age)
names(idents) <- rownames(PRD_B@meta.data)
PRD_B@active.ident<-factor(idents)
PRD_B_markers <- FindMarkers(PRD_B, ident.1 = "M9_B", ident.2 = "M6_B", min.pct = 0.25)
PRD_B_markers$group <- 0
PRD_B_markers[PRD_B_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_B_markers[PRD_B_markers$avg_log2FC < -0.25,]$group <- "6M"
PRD_B_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_B_markers$group))
data3$diets <- "PRD: 9m vs 6m"

PB_B_markers <- FindMarkers(M9_BP, ident.1 = "PRD_B",
                             ident.2 = "BRD_B", min.pct = 0.25)
PB_B_markers$group <- 0
PB_B_markers[PB_B_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_B_markers[PB_B_markers$avg_log2FC < -0.25,]$group <- "BRD"
PB_B_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_B_markers$group))
data2$diets <- "9m: PRD vs BRD"
data <- rbind(data1,data3)%>%rbind(.,data2)

data$celltype <- "B"
d2 <- data

B_age_up <- intersect(PRD_B_markers[PRD_B_markers$group=="9M",]%>%rownames(.),
                       PB_B_markers[PB_B_markers$group=="PRD",]%>%rownames(.))

B_age_down <-intersect(PRD_B_markers[PRD_B_markers$group=="6M",]%>%rownames(.),
                        PB_B_markers[PB_B_markers$group=="BRD",]%>%rownames(.))

########################## DC ####################################
BRD_DC <- readRDS("./BRD-DC_BP.rds")
df1=as.data.frame(BRD_DC@active.ident)
idents <- gsub("Dendritic cells","DC",BRD_DC$group_age)
names(idents) <- rownames(BRD_DC@meta.data)
BRD_DC@active.ident<-factor(idents)
BRD_DC_markers <- FindMarkers(BRD_DC, ident.1 = "M9_DC", ident.2 = "M6_DC", min.pct = 0.25)
BRD_DC_markers$group <- 0
BRD_DC_markers[BRD_DC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_DC_markers[BRD_DC_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_DC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_DC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_DC <- readRDS("./PRD-DC_BP.rds")
df1=as.data.frame(PRD_DC@active.ident)
idents <- gsub("Dendritic cells","DC",PRD_DC$group_age)
names(idents) <- rownames(PRD_DC@meta.data)
PRD_DC@active.ident<-factor(idents)
PRD_DC_markers <- FindMarkers(PRD_DC, ident.1 = "M9_DC", ident.2 = "M6_DC", min.pct = 0.25)
PRD_DC_markers$group <- 0
PRD_DC_markers[PRD_DC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_DC_markers[PRD_DC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_DC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_DC_markers$group))
data3$diets <- "PRD: 9m vs 6m"

PB_DC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_DC",
                             ident.2 = "BRD_DC", min.pct = 0.25)
PB_DC_markers$group <- 0
PB_DC_markers[PB_DC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_DC_markers[PB_DC_markers$avg_log2FC < -0.25,]$group <- "BRD"
PB_DC_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_DC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "DC"
d3 <- data
DC_age_up <- intersect(PRD_DC_markers[PRD_DC_markers$group=="9M",]%>%rownames(.),
                       PB_DC_markers[PB_DC_markers$group=="PRD",]%>%rownames(.))

DC_age_down <-intersect(PRD_DC_markers[PRD_DC_markers$group=="6M",]%>%rownames(.),
                        PB_DC_markers[PB_DC_markers$group=="BRD",]%>%rownames(.))

########################## EDC ####################################
BRD_EDC <- readRDS("./BRD-EDC_BP.rds")
df1=as.data.frame(BRD_EDC@active.ident)
idents <- gsub("Endothelial cells","EDC",BRD_EDC$group_age)
names(idents) <- rownames(BRD_EDC@meta.data)
BRD_EDC@active.ident<-factor(idents)
BRD_EDC_markers <- FindMarkers(BRD_EDC, ident.1 = "M9_EDC", ident.2 = "M6_EDC", min.pct = 0.25)
BRD_EDC_markers$group <- 0
BRD_EDC_markers[BRD_EDC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_EDC_markers[BRD_EDC_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_EDC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_EDC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_EDC <- readRDS("./PRD-EDC_BP.rds")
df1=as.data.frame(PRD_EDC@active.ident)
idents <- gsub("Endothelial cells","EDC",PRD_EDC$group_age)
names(idents) <- rownames(PRD_EDC@meta.data)
PRD_EDC@active.ident<-factor(idents)
PRD_EDC_markers <- FindMarkers(PRD_EDC, ident.1 = "M9_EDC", ident.2 = "M6_EDC", min.pct = 0.25)
PRD_EDC_markers$group <- 0
PRD_EDC_markers[PRD_EDC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_EDC_markers[PRD_EDC_markers$avg_log2FC < -0.25,]$group <- "6M"
PRD_EDC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_EDC_markers$group))
data3$diets <- "PRD: 9m vs 6m"

PB_EDC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_EDC",
                             ident.2 = "BRD_EDC", min.pct = 0.25)
PB_EDC_markers$group <- 0
PB_EDC_markers[PB_EDC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_EDC_markers[PB_EDC_markers$avg_log2FC < -0.25,]$group <- "BRD"
PB_EDC_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_EDC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "EDC"
d4 <- data
EDC_age_up <- intersect(PRD_EDC_markers[PRD_EDC_markers$avg_log2FC>0,]%>%rownames(.),
                       PB_EDC_markers[PB_EDC_markers$avg_log2FC>0,]%>%rownames(.))

EDC_age_down <-intersect(PRD_EDC_markers[PRD_EDC_markers$avg_log2FC<0,]%>%rownames(.),
                        PB_EDC_markers[PB_EDC_markers$avg_log2FC<0,]%>%rownames(.))


#PRD_EDC_markers[PRD_EDC_markers$avg_log2FC > 0.25,]%>%rownames()

ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=EDC_age_up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")

KEGG_EDC_resistance <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=EDC_age_down,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_EDC_rescue <- enrichKEGG(gene =ncbi$ENTREZID,
                                  organism = 'mmu',
                                  pvalueCutoff =0.05)

KEGG_EDC_rescue@result
KEGG_EDC_rescue@resu
write.csv(KEGG_EDC_resistance@result,"./KEGG_EDC_resistance.csv")
write.csv(KEGG_EDC_rescue@result,"./KEGG_EDC_rescue.csv")









########################## EC ####################################
BRD_EC <- readRDS("./BRD-EC_BP.rds")
df1=as.data.frame(BRD_EC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_EC@meta.data)
BRD_EC@active.ident<-factor(idents)
BRD_EC_markers <- FindMarkers(BRD_EC, ident.1 = "M9_EC", ident.2 = "M6_EC", min.pct = 0.25)
BRD_EC_markers$group <- 0
BRD_EC_markers[BRD_EC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_EC_markers[BRD_EC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_EC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_EC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_EC <- readRDS("./PRD-EC_BP.rds")
df1=as.data.frame(PRD_EC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_EC@meta.data)
PRD_EC@active.ident<-factor(idents)
PRD_EC_markers <- FindMarkers(PRD_EC, ident.1 = "M9_EC", ident.2 = "M6_EC", min.pct = 0.25)
PRD_EC_markers$group <- 0
PRD_EC_markers[PRD_EC_markers$avg_log2FC >0,]$group <- "9M"
PRD_EC_markers[PRD_EC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_EC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_EC_markers$group))
data3$diets <- "PRD: 9m vs 6m"

PB_EC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_EC",
                             ident.2 = "BRD_EC", min.pct = 0.25)
PB_EC_markers$group <- 0
PB_EC_markers[PB_EC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_EC_markers[PB_EC_markers$avg_log2FC < -0.25,]$group <- "BRD"
PB_EC_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_EC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "EC"
d6 <- data
EC_age_up <- intersect(PRD_EC_markers[PRD_EC_markers$group=="9M",]%>%rownames(.),
                       PB_EC_markers[PB_EC_markers$group=="PRD",]%>%rownames(.))

EC_age_down <-intersect(PRD_EC_markers[PRD_EC_markers$group=="6M",]%>%rownames(.),
                        PB_EC_markers[PB_EC_markers$group=="BRD",]%>%rownames(.))

########################## T ####################################
BRD_T <- readRDS("./BRD-T_BP.rds")
df1=as.data.frame(BRD_T@active.ident)
idents <- gsub("T cells","T",BRD_T$group_age)
names(idents) <- rownames(BRD_T@meta.data)
BRD_T@active.ident<-factor(idents)
BRD_T_markers <- FindMarkers(BRD_T, ident.1 = "M9_T", ident.2 = "M6_T", min.pct = 0.25)
BRD_T_markers$group <- 0
BRD_T_markers[BRD_T_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_T_markers[BRD_T_markers$avg_log2FC < -0.25,]$group <- "6M"
BRD_T_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_T_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_T <- readRDS("./PRD-T_BP.rds")
df1=as.data.frame(PRD_T@active.ident)
idents <- gsub("T cells","T",PRD_T$group_age)
names(idents) <- rownames(PRD_T@meta.data)
PRD_T@active.ident<-factor(idents)
PRD_T_markers <- FindMarkers(PRD_T, ident.1 = "M9_T", ident.2 = "M6_T", min.pct = 0.25)
PRD_T_markers$group <- 0
PRD_T_markers[PRD_T_markers$avg_log2FC >0,]$group <- "9M"
PRD_T_markers[PRD_T_markers$avg_log2FC < -0.25,]$group <- "6M"
PRD_T_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_T_markers$group))
data3$diets <- "PRD: 9m vs 6m"
PB_T_markers <- FindMarkers(M9_BP, ident.1 = "PRD_T",
                             ident.2 = "BRD_T", min.pct = 0.25)
PB_T_markers$group <- 0
PB_T_markers[PB_T_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_T_markers[PB_T_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_T_markers$diets <- "9m: PRD vs BRD"
data2 <- as.data.frame(table(PB_T_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "T"
d7 <- data
T_age_up <- intersect(PRD_T_markers[PRD_T_markers$group=="9M",]%>%rownames(.),
                       PB_T_markers[PB_T_markers$group=="PRD",]%>%rownames(.))

T_age_down <-intersect(PRD_T_markers[PRD_T_markers$group=="6M",]%>%rownames(.),
                        PB_T_markers[PB_T_markers$group=="BRD",]%>%rownames(.))
################### 合并 ################################
data <- rbind(d1,d2)%>%rbind(.,d3)%>%
  rbind(.,d4)%>%rbind(.,d5)%>%rbind(.,d6)%>%rbind(.,d7)
head(data)
data$diets <- factor(data$diets,
                     levels = c("BRD: 9m vs 6m",
                                "PRD: 9m vs 6m",
                                "9m: PRD vs BRD"))
data$celltype <- factor(data$celltype,levels = c("GC","EDC","DC","B","SC","EC","T"))

#"#F39EA2","#EC2124"

data$stage<- data$Var1
data$condition <- data$stage
ggplot(data[data$stage%in%c("BRD","PRD","9M","6M"),], aes(x =celltype,y = Freq,fill = condition))+
  geom_bar(stat="identity",position = "stack",color="black")+
  facet_wrap(~diets)+theme_classic()+
  scale_fill_manual(values = c("#C4E0ED","#4E7DB8","#8A9197","#D2AF81"))+
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18,color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust=1),	
        legend.text = element_text(size=12),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(color = 'black'),axis.line.y=element_line(color = 'black'))+
  ylab("N of DEGs")+
  theme(strip.text.x = element_text(size=14))

unique(c(GC_age_up,SC_age_up,B_age_up,DC_age_up,EDC_age_up,EC_age_up,T_age_up))%>%length()
intersect(GC_age_up,SC_age_up)%>%intersect(.,B_age_up)%>%intersect(.,DC_age_up)%>%
  intersect(.,EDC_age_up)%>%intersect(.,EC_age_up)%>%intersect(.,T_age_up)

unique(c(GC_age_down,SC_age_down,B_age_down,DC_age_down,EDC_age_down,EC_age_down,T_age_down))%>%length()

intersect(GC_age_down,SC_age_down)%>%intersect(.,B_age_down)%>%intersect(.,DC_age_down)%>%
  intersect(.,EDC_age_down)%>%intersect(.,EC_age_down)%>%intersect(.,T_age_down)



M9_BP <- readRDS("./M9_BP.rds")
M9_BP@active.ident <- factor(gsub("-","_",M9_BP$diet_cell))
########################## GLC ####################################
BRD_GLC <- readRDS("./BRD-GLC_BP.rds")
df1=as.data.frame(BRD_GLC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_GLC@meta.data)
BRD_GLC@active.ident<-factor(idents)
BRD_GLC_markers <- FindMarkers(BRD_GLC, ident.1 = "M9_GLC", ident.2 = "M6_GLC", min.pct = 0.25)
BRD_GLC_markers$group <- 0
BRD_GLC_markers[BRD_GLC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_GLC_markers[BRD_GLC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_GLC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_GLC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_GLC <- readRDS("./PRD-GLC_BP.rds")
df1=as.data.frame(PRD_GLC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_GLC@meta.data)
PRD_GLC@active.ident<-factor(idents)
PRD_GLC_markers <- FindMarkers(PRD_GLC, ident.1 = "M9_GLC", ident.2 = "M6_GLC", min.pct = 0.25)
PRD_GLC_markers$group <- 0
PRD_GLC_markers[PRD_GLC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_GLC_markers[PRD_GLC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_GLC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_GLC_markers$group))
data3$diets <- "PRD: 9m vs 6m"
M9_BP@active.ident
PB_GLC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_GLC",
                             ident.2 = "BRD_GLC", min.pct = 0.25)
PB_GLC_markers$group <- 0
PB_GLC_markers[PB_GLC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_GLC_markers[PB_GLC_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_GLC_markers$diets <- "9m: PRD vs BRD"

data2 <- as.data.frame(table(PB_GLC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "GLC"
d1 <- data
GLC_resistance <- intersect(PRD_GLC_markers[PRD_GLC_markers$group=="9M",]%>%rownames(.),
                       PB_GLC_markers[PB_GLC_markers$group=="PRD",]%>%rownames(.))

GLC_rescue <-intersect(PRD_GLC_markers[PRD_GLC_markers$group=="6M",]%>%rownames(.),
                        PB_GLC_markers[PB_GLC_markers$group=="BRD",]%>%rownames(.))



########################## ILC ####################################
BRD_ILC <- readRDS("./BRD-ILC_BP.rds")
df1=as.data.frame(BRD_ILC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_ILC@meta.data)
BRD_ILC@active.ident<-factor(idents)
BRD_ILC_markers <- FindMarkers(BRD_ILC, ident.1 = "M9_ILC", ident.2 = "M6_ILC", min.pct = 0.25)
BRD_ILC_markers$group <- 0
BRD_ILC_markers[BRD_ILC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_ILC_markers[BRD_ILC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_ILC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_ILC_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_ILC <- readRDS("./PRD-ILC_BP.rds")
df1=as.data.frame(PRD_ILC@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_ILC@meta.data)
PRD_ILC@active.ident<-factor(idents)
PRD_ILC_markers <- FindMarkers(PRD_ILC, ident.1 = "M9_ILC", ident.2 = "M6_ILC", min.pct = 0.25)
PRD_ILC_markers$group <- 0
PRD_ILC_markers[PRD_ILC_markers$avg_log2FC >0,]$group <- "9M"
PRD_ILC_markers[PRD_ILC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_ILC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_ILC_markers$group))
data3$diets <- "PRD: 9m vs 6m"
M9_BP@active.ident
PB_ILC_markers <- FindMarkers(M9_BP, ident.1 = "PRD_ILC",
                              ident.2 = "BRD_ILC", min.pct = 0.25)
PB_ILC_markers$group <- 0
PB_ILC_markers[PB_ILC_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_ILC_markers[PB_ILC_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_ILC_markers$diets <- "9m: PRD vs BRD"

data2 <- as.data.frame(table(PB_ILC_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "ILC"
d2 <- data
ILC_resistance <- intersect(PRD_ILC_markers[PRD_ILC_markers$group=="9M",]%>%rownames(.),
                            PB_ILC_markers[PB_ILC_markers$group=="PRD",]%>%rownames(.))

ILC_rescue <-intersect(PRD_ILC_markers[PRD_ILC_markers$group=="6M",]%>%rownames(.),
                       PB_ILC_markers[PB_ILC_markers$group=="BRD",]%>%rownames(.))


########################## M ####################################

BRD_M <- readRDS("./BRD-M_BP.rds")
df1=as.data.frame(BRD_M@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_M@meta.data)
BRD_M@active.ident<-factor(idents)
BRD_M_markers <- FindMarkers(BRD_M, ident.1 = "M9_M", ident.2 = "M6_M", min.pct = 0.25)
BRD_M_markers$group <- 0
BRD_M_markers[BRD_M_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_M_markers[BRD_M_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_M_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_M_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_M <- readRDS("./PRD-M_BP.rds")
df1=as.data.frame(PRD_M@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_M@meta.data)
PRD_M@active.ident<-factor(idents)
PRD_M_markers <- FindMarkers(PRD_M, ident.1 = "M9_M", ident.2 = "M6_M", min.pct = 0.25)
PRD_M_markers$group <- 0
PRD_M_markers[PRD_M_markers$avg_log2FC >0,]$group <- "9M"
PRD_M_markers[PRD_M_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_M_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_M_markers$group))
data3$diets <- "PRD: 9m vs 6m"
M9_BP@active.ident
PB_M_markers <- FindMarkers(M9_BP, ident.1 = "PRD_M",
                              ident.2 = "BRD_M", min.pct = 0.25)
PB_M_markers$group <- 0
PB_M_markers[PB_M_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_M_markers[PB_M_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_M_markers$diets <- "9m: PRD vs BRD"

data2 <- as.data.frame(table(PB_M_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "M"
d3 <- data
M_resistance <- intersect(PRD_M_markers[PRD_M_markers$group=="9M",]%>%rownames(.),
                            PB_M_markers[PB_M_markers$group=="PRD",]%>%rownames(.))

M_rescue <-intersect(PRD_M_markers[PRD_M_markers$group=="6M",]%>%rownames(.),
                       PB_M_markers[PB_M_markers$group=="BRD",]%>%rownames(.))




########################## NK ####################################

BRD_NK <- readRDS("./BRD-NK_BP.rds")
df1=as.data.frame(BRD_NK@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(BRD_NK@meta.data)
BRD_NK@active.ident<-factor(idents)
BRD_NK_markers <- FindMarkers(BRD_NK, ident.1 = "M9_NK", ident.2 = "M6_NK", min.pct = 0.25)
BRD_NK_markers$group <- 0
BRD_NK_markers[BRD_NK_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_NK_markers[BRD_NK_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_NK_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_NK_markers$group))
data1$diets <- "BRD: 9m vs 6m"


PRD_NK <- readRDS("./PRD-NK_BP.rds")
df1=as.data.frame(PRD_NK@meta.data)
df1$cell <- paste0(df1$age,"_",df1$abbreviation2)
idents <- df1$cell 
names(idents) <- rownames(PRD_NK@meta.data)
PRD_NK@active.ident<-factor(idents)
PRD_NK_markers <- FindMarkers(PRD_NK, ident.1 = "M9_NK", ident.2 = "M6_NK", min.pct = 0.25)
PRD_NK_markers$group <- 0
PRD_NK_markers[PRD_NK_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_NK_markers[PRD_NK_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_NK_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_NK_markers$group))
data3$diets <- "PRD: 9m vs 6m"
M9_BP@active.ident
PB_NK_markers <- FindMarkers(M9_BP, ident.1 = "PRD_NK",
                              ident.2 = "BRD_NK", min.pct = 0.25)
PB_NK_markers$group <- 0
PB_NK_markers[PB_NK_markers$avg_log2FC > 0.25,]$group <- "PRD"
PB_NK_markers[PB_NK_markers$avg_log2FC< -0.25,]$group <- "BRD"
PB_NK_markers$diets <- "9m: PRD vs BRD"

data2 <- as.data.frame(table(PB_NK_markers$group))
data2$diets <- "9m: PRD vs BRD"

data <- rbind(data1,data3)%>%rbind(.,data2)
data$celltype <- "NK"
d4 <- data
NK_resistance <- intersect(PRD_NK_markers[PRD_NK_markers$group=="9M",]%>%rownames(.),
                            PB_NK_markers[PB_NK_markers$group=="PRD",]%>%rownames(.))

NK_rescue <-intersect(PRD_NK_markers[PRD_NK_markers$group=="6M",]%>%rownames(.),
                       PB_NK_markers[PB_NK_markers$group=="BRD",]%>%rownames(.))

data <- rbind(d1,d2)%>%rbind(.,d3)%>%
  rbind(.,d4)
data$diets <- factor(data$diets,
                     levels = c("BRD: 9m vs 6m",
                                "PRD: 9m vs 6m",
                                "9m: PRD vs BRD"))
data$celltype <- factor(data$celltype,levels = c("GLC","ILC","M","NK"))

#"#F39EA2","#EC2124"

data$stage<- data$Var1
data$condition <- data$stage
data <- data[data$Var1%in%c("6M","9M","BRD","PRD"),]
ggplot(data, aes(x =celltype,y = Freq,fill = condition))+
  geom_bar(stat="identity",position = "stack",color="black")+
  facet_wrap(~diets)+theme_classic()+
  scale_fill_manual(values = c("#C4E0ED","#4E7DB8","#8A9197","#D2AF81"))+
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18,color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust=1),	
        legend.text = element_text(size=12),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(color = 'black'),axis.line.y=element_line(color = 'black'))+
  ylab("N of DEGs")+
  theme(strip.text.x = element_text(size=14))


aging <- unique(c(PRD_NK_markers[PRD_NK_markers$group=="9M",]%>%rownames(.),
               PRD_GLC_markers[PRD_GLC_markers$group=="9M",]%>%rownames(.),
               PRD_ILC_markers[PRD_ILC_markers$group=="9M",]%>%rownames(.),
               PRD_M_markers[PRD_M_markers$group=="9M",]%>%rownames(.)))
length(aging)



young <- unique(c(PRD_NK_markers[PRD_NK_markers$group=="6M",]%>%rownames(.),
                  PRD_GLC_markers[PRD_GLC_markers$group=="6M",]%>%rownames(.),
                  PRD_ILC_markers[PRD_ILC_markers$group=="6M",]%>%rownames(.),
                  PRD_M_markers[PRD_M_markers$group=="6M",]%>%rownames(.)))
length(young)


PRD <- unique(c(PB_NK_markers[PB_NK_markers$group=="PRD",]%>%rownames(.),
                  PB_GLC_markers[PB_GLC_markers$group=="PRD",]%>%rownames(.),
                  PB_ILC_markers[PB_ILC_markers$group=="PRD",]%>%rownames(.),
                  PB_M_markers[PB_M_markers$group=="PRD",]%>%rownames(.)))
length(PRD)

BRD <- unique(c(PB_NK_markers[PB_NK_markers$group=="BRD",]%>%rownames(.),
                PB_GLC_markers[PB_GLC_markers$group=="BRD",]%>%rownames(.),
                PB_ILC_markers[PB_ILC_markers$group=="BRD",]%>%rownames(.),
                PB_M_markers[PB_M_markers$group=="BRD",]%>%rownames(.)))
length(BRD)

resistance <- unique(c(GLC_resistance,ILC_resistance,M_resistance,NK_resistance))
length(resistance)
length(GLC_resistance)
length(GLC_rescue)
length(PRD_GLC_markers[PRD_GLC_markers$group=="9M",]%>%rownames(.))
length(PRD_GLC_markers[PRD_GLC_markers$group=="6M",]%>%rownames(.))

length(ILC_resistance)
length(ILC_rescue)
length(PRD_ILC_markers[PRD_ILC_markers$group=="9M",]%>%rownames(.))
length(PRD_ILC_markers[PRD_ILC_markers$group=="6M",]%>%rownames(.))


length(M_resistance)
length(M_rescue)
length(PRD_M_markers[PRD_M_markers$group=="9M",]%>%rownames(.))
length(PRD_M_markers[PRD_M_markers$group=="6M",]%>%rownames(.))

length(NK_resistance)
length(NK_rescue)
length(PRD_NK_markers[PRD_NK_markers$group=="9M",]%>%rownames(.))
length(PRD_NK_markers[PRD_NK_markers$group=="6M",]%>%rownames(.))

rescue <- unique(c(GLC_rescue,ILC_rescue,M_rescue,NK_rescue))
length(rescue)



r1 <- c(GC_age_up,EDC_age_up,DC_age_up,B_age_up,SC_age_up,EC_age_up,T_age_up)
r2 <- c(GC_age_down,EDC_age_down,DC_age_down,B_age_down,SC_age_down,EC_age_down,T_age_down)

r3 <- c(GLC_resistance,ILC_resistance,M_resistance,NK_resistance)
r4 <- c(GLC_rescue,ILC_rescue,M_rescue,NK_rescue)
BP <- readRDS("./BP.rds")
BP_sub <- subset(BP,cells=BP@meta.data[BP$Diets%in%c("M6-PRD","M9-PRD","M9-BRD"),]%>%rownames(.))
rm(BP)
gc()

meta <-BP_sub@meta.data
n <- t(scale(t(as.data.frame(BP_sub@assays$RNA@data))))%>% na.omit(.)

avg <- apply(as.data.frame(n),1,function(x){tapply(x,meta$Diets,mean)}) %>% t() %>% as.data.frame()

avg1 <- avg[r1,c("GC","EDC","DC","B","SC","EC","T")]
avg2 <- avgt[r2,c("GC","EDC","DC","B","SC","EC","T")]

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))

h1 <- Heatmap(t(avg1),name="scale",row_dend_reorder = F,row_title = "Antral GC",
              col=col_fun,column_title = "Resistance",show_column_names  = T)
h1
dim(mean_Ant_s2)
h2 <- Heatmap(t(avg2),name="scale",row_dend_reorder = T,row_title = "Antral GC",
              col=col_fun,column_title = "Rescue",show_column_names  = T)
h1+h2






######################  Aging Pathway 1 ###########################
up <- c(GC_age_up,EDC_age_up,DC_age_up,B_age_up,
        SC_age_up,EC_age_up,T_age_up)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_up <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)

up <- c(GC_age_up,EDC_age_up,DC_age_up,B_age_up,
        SC_age_up,EC_age_up,T_age_up)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_up <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)
write.csv(KEGG_up@result,"./KEGG_age_up.csv")
down <- c(GC_age_down,EDC_age_down,DC_age_down,B_age_down,
        SC_age_down,EC_age_down,T_age_down)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=down,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_down <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)
write.csv(KEGG_down@result,"./KEGG_age_down.csv")

KEGG_age_down <- read.csv("./KEGG_age_down.csv")
s <- KEGG_age_down[KEGG_age_down$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s$Description
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD rescue") +
  geom_text(aes(x = rep(0.002,12), #用新增的重复数组控制文本标签起始位置
  label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)

KEGG_age_up <- read.csv("./KEGG_age_up.csv")
s <- KEGG_age_up[KEGG_age_up$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s$Description
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller( direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD resistance") +
  geom_text(aes(x = rep(0.002,12), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)

######################  Aging Pathway 2 ###########################

up <- c(GLC_resistance,ILC_resistance,M_resistance,NK_resistance)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_up <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)
write.csv(KEGG_up@result,"./KEGG_resistance.csv")
down <- c(GLC_rescue,ILC_rescue,M_rescue,NK_rescue)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=down,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_down <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)
write.csv(KEGG_down@result,"./KEGG_rescue.csv")




KEGG_rescue <- read.csv("./KEGG_rescue.csv")
s <- KEGG_rescue[KEGG_rescue$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s$Description
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD rescue") +
  geom_text(aes(x = rep(0.002,12), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)

KEGG_resistance <- read.csv("./KEGG_resistance.csv")
s <- KEGG_resistance[KEGG_resistance$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s$Description
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller( direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD resistance") +
  geom_text(aes(x = rep(0.002,12), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)


PRD_GC_markers[PRD_GC_markers$avg_log2FC>0,]

ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_GC_markers[PRD_GC_markers$avg_log2FC<0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC <- enrichKEGG(gene =ncbi$ENTREZID,
                       organism = 'mmu',
                       pvalueCutoff =0.05)



ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_DC_markers[PRD_DC_markers$avg_log2FC>0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_DC <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)



ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_EDC_markers[PRD_EDC_markers$avg_log2FC>0.25,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_EDC <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_B_markers[PRD_B_markers$avg_log2FC>0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_B <- enrichKEGG(gene =ncbi$ENTREZID,
                       organism = 'mmu',
                       pvalueCutoff =0.05)


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_SC_markers[PRD_SC_markers$avg_log2FC>0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_SC <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_EC_markers[PRD_EC_markers$avg_log2FC>0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_EC <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=rownames(PRD_T_markers[PRD_T_markers$avg_log2FC>0,]),
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_T <- enrichKEGG(gene =ncbi$ENTREZID,
                       organism = 'mmu',
                       pvalueCutoff =0.05)

write.csv(KEGG_GC@result, 'KEGG_GC_BP.csv')
write.csv(KEGG_B@result, 'KEGG_B_BP.csv')
write.csv(KEGG_DC@result, 'KEGG_DC_BP.csv')
write.csv(KEGG_SC@result, 'KEGG_SC_BP.csv')
write.csv(KEGG_EDC@result, 'KEGG_EDC_BP.csv')
write.csv(KEGG_T@result, 'KEGG_T_BP.csv')
write.csv(KEGG_EC@result, 'KEGG_EC_BP.csv')

d1 <- KEGG_GC@result
d1$Class <- "GC"
d2 <- KEGG_EDC@result
d2$Class <- "EDC"
d3 <- KEGG_DC@result
d3$Class <- "DC"
d4 <- KEGG_B@result
d4$Class <- "B"

d5 <- KEGG_SC@result
d5$Class <- "SC"
d6 <- KEGG_EC@result
d6$Class <- "EC"
d7 <- KEGG_T@result
d7$Class <- "T"

all <- rbind(d1,d2)%>%rbind(.,d3)%>%rbind(.,d4)%>%rbind(.,d5)%>%rbind(.,d6)%>%rbind(.,d7)

all$Description <- strsplit2(all$Description," - ")[,1]%>%
  gsub("Chemical carcinogenesis - reactive oxygen species","Reactive oxygen species",.)


View(d1)
KEGG_select_BP <- read.csv("~/scRNA_aging/output/KEGG_select_BP.csv")
KEGG_select_BP$pvalue
KEGG_select_BP$logp <--log10(KEGG_select_BP$pvalue)
ggplot(KEGG_select_BP,aes(x=Class,y=Description))+
  geom_tile(aes(fill=logp))+
  theme(panel.background = element_blank())+scale_x_discrete(limits=c("GC","EDC","DC","B","SC","EC","T"))+xlab("")+
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 14,color = 'black'))+	
  theme(axis.text = element_text(size = 14),axis.text.x = element_text(angle = 45,hjust=1),	
        legend.text = element_text(size=13),legend.position = "right",
        legend.title = element_text(size=13),
        axis.line.x=element_line(),axis.line.y=element_line())+
  scale_fill_gradient2(low="white",high = "#E91C22")+labs(fill="-log10(pvalue)")+scale_y_discrete(limits=unique(KEGG_select_BP$Description))

KEGG_select_BP <- all[all$Description%in%unique(KEGG_select_BP$Description),]
KEGG_select_BP$RichFactor=as.numeric(KEGG_select_BP$Count)/ as.numeric(sub("/\\d+", "",KEGG_select_BP$BgRatio))
KEGG_select_BP$Class
ggplot(KEGG_select_BP,aes(Class,Description)) + geom_point(aes(color=pvalue,size=RichFactor))+						
  scale_colour_gradient(low="#FF603F",high= "#6997ED")+						
  labs(size="RichFactor",x="",y="",title="")+					
  theme(axis.ticks = element_blank())+
  scale_y_discrete(limits=rev(unique(KEGG_select_BP$Description)))+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust=1),	
        legend.text = element_text(size=12),legend.title  = element_text(size=12),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  scale_x_discrete(limits=c("GC","EDC","DC","B","SC","EC","T"))

a=KEGG_select_BP[KEGG_select_BP$Description%in%c("Cellular senescence","Reactive oxygen species"),]$geneID
b=KEGG_select_BP[KEGG_select_BP$Description%in%c("Cellular senescence"),]$geneID
b


strsplit2(a,"/")%>%c()


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=strsplit2(a,"/")%>%c(),
                             columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")
ncbi
na.omit(ncbi)
BP <- readRDS("./BP.rds")










cell=BP@meta.data[BP@meta.data$Diets%in%c("M6-PRD"),]%>%rownames()
M6 <- subset(BP,cells=cell)
FeaturePlot(M6,features = "Foxp1")
length(colnames(M6))
#M6_select <- subset(M6,cells = sample(colnames(M6),40000))


cellname <- M6$abbreviation2
names(cellname) <- colnames(M6)
M6@active.ident <- factor(cellname)
M6_markers <- FindAllMarkers(object =M6,only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)


write.csv(M6_markers,"M6_markers.csv")
saveRDS(M6,"M6.rds")

factor(M6$abbreviation2)



######################### Cell Proportion ############################
BP <- readRDS("./BP_sub.rds")
meta <- BP@meta.data[BP@meta.data$Diets%in%c("M6-PRD","M9-PRD","M9-BRD"),]
levels(factor(meta$abbreviation2))
rm(BP)
gc()
View(meta)
meta$orig.ident
meta$abbreviation2
data <- meta%>%dplyr::group_by(orig.ident)%>%
  summarise(n=table(abbreviation2))%>%as.data.frame()
#View(data)
data$celltype <- rep(levels(factor(meta$abbreviation2)),8)
#View(data)
data1 <- data%>% group_by(orig.ident)%>%summarise(freq=n/sum(n))%>%as.data.frame()
head(data1)
data1$celltype <- rep(levels(factor(meta$abbreviation2)),8)
??spread
data2 <- spread(data1,key="orig.ident",value = "freq")

data3 <- apply(data2[,2:9],1,function(x){tapply(x,c(rep("M6-PRD",3),rep("M9-BRD",3),rep("M9-PRD",2)),mean)}) %>% t() %>% as.data.frame()
data3$celltype <- data2$celltype
write.csv(data3,"cell_prop.csv")

data <- melt(data3)
head(data)
data$variable <- factor(data$variable,levels = c("M6-PRD","M9-PRD","M9-BRD")) 
data
ggplot(data, aes( x = variable,y=100 * value,fill = celltype,
                 stratum =  celltype, alluvium =  celltype))+
  geom_stratum(width = 0.8, color='white')+
  geom_alluvium(alpha = 0.5,
                width = 0.8,
                curve_type = "linear")+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust=1),	
        legend.text = element_text(size=12),legend.title  = element_text(size=12),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+xlab("")+
  ylab("Cell Proportion (%)")+scale_fill_simpsons()

BP






install.packages("ggalluvial")
library(ggalluvial)
aging_gene_sets <- read.delim("~/scRNA_aging/output/aging_gene_sets.txt")

#View(aging_gene_sets)

head(aging_gene_sets)
gene1 <- as.list(aging_gene_sets[,1])[1:34]
gene2 <- as.list(aging_gene_sets[,2])
gene3 <- as.list(aging_gene_sets[,3])[1:40]
gene4 <- as.list(aging_gene_sets[,4])[1:27][-c(20,8)]

gene5 <- as.list(unique(c("Gdf9","Bmp15","Foxl2","Gja1","Fshr","Lhcgr",
                   "Lhcgr","Esr2","Fshr","Inha","Amh",
                   "Star","Cyp17a1","Cyp19a1")))



############################# 八大衰老通路 #########################
library(readxl)
Apoptosis <- read_excel("../Aging_genes/Apoptosis.xlsx")
Autophagy <- read_excel("../Aging_genes/Autophagy.xlsx")
ROS <- read_excel("../Aging_genes/ROS.xlsx")
SASP <- read_excel("../Aging_genes/SASP.xlsx")
DNArepair <- read_excel("../Aging_genes/DNA repair.xlsx")
Fibrosis <- read_excel("../Aging_genes/Fibrosis.xlsx")
Lipid_biosynthesis <- read_excel("../Aging_genes/Lipid biosynthesis.xlsx")
Senescence <- read_excel("../Aging_genes/Senescence.xlsx")
UPR <- read_excel("../Aging_genes/UPR.xlsx")



MtoH <- readRDS("./MtoH.rds")

Apoptosis <- MtoH[MtoH$HGNC.symbol%in%Apoptosis$`Apoptosis-related gene set`,]$MGI.symbol
Autophagy <- MtoH[MtoH$HGNC.symbol%in%Autophagy$`Autophagy-related gene set`,]$MGI.symbol
ROS <- MtoH[MtoH$HGNC.symbol%in%ROS$`ROS-related gene set`,]$MGI.symbol
SASP <- MtoH[MtoH$HGNC.symbol%in%SASP$`SASP-related gene set`,]$MGI.symbol
DNArepair <- MtoH[MtoH$HGNC.symbol%in%DNArepair$`DNA repair-related gene set`,]$MGI.symbol
Fibrosis <- MtoH[MtoH$HGNC.symbol%in%Fibrosis$`Fibrosis-related gene set`,]$MGI.symbol
Lipid_biosynthesis <- MtoH[MtoH$HGNC.symbol%in%Lipid_biosynthesis$`Lipid biosynthesis related genes`,]$MGI.symbol
Senescence <- MtoH[MtoH$HGNC.symbol%in%Senescence$`Senescence-related gene set`,]$MGI.symbol
UPR <- MtoH[MtoH$HGNC.symbol%in%UPR$`UPR pathway -related gene sets`,]$MGI.symbol




################################## Seurat Scores ###########################################
BP <- AddModuleScore(object = BP,features = gene1,ctrl = 100,name = "Cell_Cycle")
BP <- AddModuleScore(object = BP,features = gene2,ctrl = 100,name = "DNA_Repair")
BP <- AddModuleScore(object = BP,features = gene3,ctrl = 100,name = "SASP")
BP <- AddModuleScore(object = BP,features = gene4,ctrl = 100,name = "Inflammation")

BP <- AddModuleScore(object = BP,features = Apoptosis,ctrl = 100,name = "Apoptosis")
BP <- AddModuleScore(object = BP,features = Autophagy,ctrl = 100,name = "Autophagy")
BP <- AddModuleScore(object = BP,features = ROS,ctrl = 100,name = "ROS")
BP <- AddModuleScore(object = BP,features = SASP,ctrl = 100,name = "SASP")
BP <- AddModuleScore(object = BP,features = DNArepair,ctrl = 100,name = "DNArepair")
BP <- AddModuleScore(object = BP,features = Fibrosis,ctrl = 100,name = "Fibrosis")
BP <- AddModuleScore(object = BP,features = Lipid_biosynthesis,ctrl = 100,name = "Lipid_biosynthesis")
BP <- AddModuleScore(object = BP,features = Senescence,ctrl = 100,name = "Senescence")
BP <- AddModuleScore(object = BP,features = UPR,ctrl = 100,name = "UPR")


BP$Group <- paste0(substr(BP$Diets,4,6),"-",BP$abbreviation2)
BP$Diets
data1 <- FetchData(BP,vars = c("Group","Diets","Cell_Cycle1"))
data2 <- FetchData(BP,vars = c("Group","Diets","DNA_Repair1"))
data3 <- FetchData(BP,vars = c("Group","Diets","SASP1"))
data4 <- FetchData(BP,vars = c("Group","Diets","Inflammation1"))

data5 <- FetchData(BP,vars = c("Group","Diets","Apoptosis1"))
data6 <- FetchData(BP,vars = c("Group","Diets","Autophagy1"))
data7 <- FetchData(BP,vars = c("Group","Diets","ROS1"))
#data8 <- FetchData(BP,vars = c("Group","Diets","SASP1"))
data9 <- FetchData(BP,vars = c("Group","Diets","DNArepair1"))
data10 <- FetchData(BP,vars = c("Group","Diets","Fibrosis1"))
data11 <- FetchData(BP,vars = c("Group","Diets","Lipid_biosynthesis1"))
data12 <- FetchData(BP,vars = c("Group","Diets","Senescence1"))
data13 <- FetchData(BP,vars = c("Group","Diets","UPR1"))

head(data1)
##################### GC1 衰老途径 #############################
GC1 <- AddModuleScore(object = GC1,features = gene1,ctrl = 100,name = "Cell_Cycle")
GC1 <- AddModuleScore(object = GC1,features = gene2,ctrl = 100,name = "DNA_Repair")
GC1 <- AddModuleScore(object = GC1,features = gene3,ctrl = 100,name = "SASP")
GC1 <- AddModuleScore(object = GC1,features = gene4,ctrl = 100,name = "Inflammation")

GC1 <- AddModuleScore(object = GC1,features = Apoptosis,ctrl = 100,name = "Apoptosis")
GC1 <- AddModuleScore(object = GC1,features = Autophagy,ctrl = 100,name = "Autophagy")
GC1 <- AddModuleScore(object = GC1,features = ROS,ctrl = 100,name = "ROS")
GC1 <- AddModuleScore(object = GC1,features = SASP,ctrl = 100,name = "SASP")
GC1 <- AddModuleScore(object = GC1,features = DNArepair,ctrl = 100,name = "DNArepair")
GC1 <- AddModuleScore(object = GC1,features = Fibrosis,ctrl = 100,name = "Fibrosis")
GC1 <- AddModuleScore(object = GC1,features = Lipid_biosynthesis,ctrl = 100,name = "Lipid_biosynthesis")
GC1 <- AddModuleScore(object = GC1,features = Senescence,ctrl = 100,name = "Senescence")
GC1 <- AddModuleScore(object = GC1,features = UPR,ctrl = 100,name = "UPR")
GC1 <- AddModuleScore(object = GC1,features = gene5,ctrl = 100,name = "Functions")
GC1$Funtions5

##################### SC1 衰老途径 #############################
SC1 <- AddModuleScore(object = SC1,features = gene1,ctrl = 100,name = "Cell_Cycle")
SC1 <- AddModuleScore(object = SC1,features = gene2,ctrl = 100,name = "DNA_Repair")
SC1 <- AddModuleScore(object = SC1,features = gene3,ctrl = 100,name = "SASP")
SC1 <- AddModuleScore(object = SC1,features = gene4,ctrl = 100,name = "Inflammation")

SC1 <- AddModuleScore(object = SC1,features = Apoptosis,ctrl = 100,name = "Apoptosis")
SC1 <- AddModuleScore(object = SC1,features = Autophagy,ctrl = 100,name = "Autophagy")
SC1 <- AddModuleScore(object = SC1,features = ROS,ctrl = 100,name = "ROS")
SC1 <- AddModuleScore(object = SC1,features = SASP,ctrl = 100,name = "SASP")
SC1 <- AddModuleScore(object = SC1,features = DNArepair,ctrl = 100,name = "DNArepair")
SC1 <- AddModuleScore(object = SC1,features = Fibrosis,ctrl = 100,name = "Fibrosis")
SC1 <- AddModuleScore(object = SC1,features = Lipid_biosynthesis,ctrl = 100,name = "Lipid_biosynthesis")
SC1 <- AddModuleScore(object = SC1,features = Senescence,ctrl = 100,name = "Senescence")
SC1 <- AddModuleScore(object = SC1,features = UPR,ctrl = 100,name = "UPR")
SC1 <- AddModuleScore(object = SC1,features = gene5,ctrl = 100,name = "Functions")

SC1$Group <- paste0(paste0(substr(colnames(SC1),1,4),"RD"),"-",SC1$abbreviation2)
SC1$Diets <- paste0(substr(colnames(SC1),1,4),"RD")
data1 <- FetchData(SC1,vars = c("Group","Diets","Cell_Cycle1"))
data2 <- FetchData(SC1,vars = c("Group","Diets","DNA_Repair1"))
data3 <- FetchData(SC1,vars = c("Group","Diets","SASP1"))
data4 <- FetchData(SC1,vars = c("Group","Diets","Inflammation1"))

data5 <- FetchData(SC1,vars = c("Group","Diets","Apoptosis1"))
data6 <- FetchData(SC1,vars = c("Group","Diets","Autophagy1"))
data7 <- FetchData(SC1,vars = c("Group","Diets","ROS1"))
data8 <- FetchData(SC1,vars = c("Group","Diets","SASP1"))
data9 <- FetchData(SC1,vars = c("Group","Diets","DNArepair1"))
data10 <- FetchData(SC1,vars = c("Group","Diets","Fibrosis1"))
data11 <- FetchData(SC1,vars = c("Group","Diets","Lipid_biosynthesis1"))
data12 <- FetchData(SC1,vars = c("Group","Diets","Senescence1"))
data13 <- FetchData(SC1,vars = c("Group","Diets","UPR1"))

data14 <- FetchData(SC1,vars = c("Group","Diets","Functions5"))

head(data1)
##############################################################################
d1 <- data.frame(Diets=data1$Diets,Group=data1$Group,
                   Celltype=strsplit2(data1$Group,"-")[,2],Scores=data1$Cell_Cycle1)

d1$Class <- "Cell Cycle"

d2 <- data.frame(Diets=data2$Diets,Group=data2$Group,
                 Celltype=strsplit2(data2$Group,"-")[,2],Scores=data2$DNA_Repair1)
d2$Class <- "DNA Repair"

d3 <- data.frame(Diets=data3$Diets,Group=data3$Group,
                 Celltype=strsplit2(data3$Group,"-")[,2],Scores=data3$SASP1)
d3$Class <- "SASP"

d4 <- data.frame(Diets=data4$Diets,Group=data4$Group,
                 Celltype=strsplit2(data4$Group,"-")[,2],Scores=data4$Inflammation1)
d4$Class <- "Inflammation"

data <- rbind(d1,d2)%>%rbind(.,d3)%>%rbind(.,d4)

################ 
d5 <- data.frame(Diets=data5$Diets,Group=data5$Group,
                 Celltype=strsplit2(data5$Group,"-")[,2],Scores=data5$Apoptosis1)
d5$Class <- "Apoptosis"
d6 <- data.frame(Diets=data6$Diets,Group=data6$Group,
                 Celltype=strsplit2(data6$Group,"-")[,2],Scores=data6$Autophagy1)
d6$Class <- "Autophagy"
d7 <- data.frame(Diets=data7$Diets,Group=data7$Group,
                 Celltype=strsplit2(data7$Group,"-")[,2],Scores=data7$ROS1)
d7$Class <- "ROS"
d8 <- data.frame(Diets=data8$Diets,Group=data8$Group,
                 Celltype=strsplit2(data8$Group,"-")[,2],Scores=data8$SASP1)
d8$Class <- "SASP"

d9 <- data.frame(Diets=data9$Diets,Group=data9$Group,
                 Celltype=strsplit2(data9$Group,"-")[,2],Scores=data9$DNArepair1)
d9$Class <- "DNA repair"
d10 <- data.frame(Diets=data10$Diets,Group=data10$Group,
                 Celltype=strsplit2(data10$Group,"-")[,2],Scores=data10$Fibrosis1)
d10$Class <- "Fibrosis"
d11 <- data.frame(Diets=data11$Diets,Group=data11$Group,
                  Celltype=strsplit2(data11$Group,"-")[,2],Scores=data11$Lipid_biosynthesis1)
d11$Class <- "Lipid biosynthesis"
d12 <- data.frame(Diets=data12$Diets,Group=data12$Group,
                  Celltype=strsplit2(data12$Group,"-")[,2],Scores=data12$Senescence1)
d12$Class <- "Senescence"
d13 <- data.frame(Diets=data13$Diets,Group=data13$Group,
                  Celltype=strsplit2(data13$Group,"-")[,2],Scores=data13$UPR1)
d13$Class <- "UPR"


d14 <- data.frame(Diets=data14$Diets,Group=data14$Group,
                  Celltype=strsplit2(data14$Group,"-")[,2],Scores=data14$Functions5)
d14$Class <- "Functions"


data <- rbind(d5,d6)%>%rbind(.,d7)%>%rbind(.,d8)%>%rbind(.,d9)%>%rbind(.,d10)%>%
  rbind(.,d11)%>%rbind(.,d12)%>%rbind(.,d13)





remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
d1$Group <- paste0(d1$Diets,"-",d1$Celltype)
d <- d1 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
#d <- d[d$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]
ggplot(d1, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)

ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
mean(ddd[ddd$Diets=="M6-PRD",]$Scores)
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)
#wilcox.test(ddd[ddd$Diets=="M6-PRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Cell cycle gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.49,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+annotate("text", label = "p = 0.0008909", x = 1, y = 2.5, size =8, colour = "black",face = "italic")


####SASP
d3$Group <- paste0(d3$Diets,"-",d3$Celltype)
d <- d3 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
d$Celltype <- factor(d$Celltype,levels = c("GC","EDC","DC","B"))

library(ggrepel)
#提取umap坐标数据
umap <- data.frame(BP@meta.data, BP@reductions$umap@cell.embeddings)
head(umap)
umap$abbreviation2
#1）计算每个celltype的median坐标位置

umap <- umap[umap$Diets%in%c("M9-BRD","M9-PRD"),]
head(umap)
cell_type_med <- umap %>%
  group_by(abbreviation2) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
#2）使用ggrepel包geom_label_repel 添加注释
umap$SASP  <- umap$SASP1
umap$DNA_Repair <- umap$DNA_Repair1
ggplot(umap, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = SASP),size=0.2) + 
  theme_classic()+facet_wrap(~Diets)+
  scale_color_continuous(high="red",low="white")+
  geom_point(aes(colour  = DNA_Repair),size=0.2) 

#dd <- d3[d3$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]
M9 <- readRDS("./M9_BP.rds")
M9$abbreviation2
DotPlot(object = M9,group.by = "diet",
        cols = c("#4B549B","#E91C22"), 
        features = unlist(gene2))+
  theme_classic()+
#  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
#                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.position = "top",
        legend.title = element_text(size = 13))+xlab("SASP gene set")

M9$abbreviation2
?DotPlot
DotPlot(object = M9,group.by = "abbreviation2",split.by = "diet",
        cols = c("#4B549B","#E91C22"), 
        features = unlist(gene3))+
  theme_classic()+
  #  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
  #                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.position = "top",
        legend.title = element_text(size = 13))+NoLegend()


ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("SASP 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7,scales = "free_y")
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets
ddd$Diets <- factor(ddd$Diets,levels = c("M9-BRD","M9-PRD"))
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)
ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "SASP gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.49,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.52,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 0.0008909", x = 1, y = 2.5, size =8, colour = "black",face = "italic")


####DNA repair
d2$Group <- paste0(d2$Diets,"-",d2$Celltype)
d <- d2 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
#d <- d[d$Celltype%in%c("SC","EC","T"),]
#d$Celltype <- factor(d$Celltype,levels = c("GC","EDC","DC","B","SC","EC","T"))
#dd <- d2[d2$Celltype%in%c("SC","EC","T"),]

ggplot(dd, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("DNA repair 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7,scales = "free_y")

ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)
ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "DNA repair gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.54,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.48,linetype="dashed",color="#00BFC4")+annotate("text", label = "p < 2.2e-16", x = 1, y =3, size =8, colour = "black",face = "italic")

####Inflammation
d4$Group <- paste0(d4$Diets,"-",d4$Celltype)
d <- d4 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
d <- d[d$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]
d$Celltype <- factor(d$Celltype,levels = c("GC","EDC","DC","B","SC","EC","T"))
dd <- d4[d4$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]

ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Inflammation
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7,scales = "free_y")
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)
ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Inflammation gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.73,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.84,linetype="dashed",color="#00BFC4")+annotate("text", label = "p = 0.0008909", x = 1, y = 2.5, size =8, colour = "black",face = "italic")

dat<- data.frame(BP@meta.data, 
                 BP@reductions$umap@cell.embeddings,
                 BP@active.ident)
colnames(dat)[ncol(dat)] = "seurat_annotation"
class_avg <- dat %>%
  group_by(seurat_annotation) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

library(ggpubr)
ggplot(dat, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = Cluster1)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = seurat_annotation),
                            data = class_avg,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + theme_bw()

############################# Apoptosis  ###################3
library(dplyr)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
d5$Group <- paste0(d5$Diets,"-",d5$Celltype)
d <- d5 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Apoptosis gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.49,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+annotate("text", label = "p = 0.0008909", x = 1, y = 2.5, size =8, colour = "black",face = "italic")


############################# Autophagy  ###################3
library(dplyr)
d6$Group <- paste0(d6$Diets,"-",d6$Celltype)
d <- d6%>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))#+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets
ddd$Diets <- factor(ddd$Diets,levels = c("M9-BRD","M9-PRD"))
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)


wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)


median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Autophagy gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.61))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.55,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 0.003397", x = 1, y = 4, size =8, colour = "black",face = "italic")

############################# ROS  ###################3
library(dplyr)
d7$Group <- paste0(d7$Diets,"-",d7$Celltype)
d <- d7 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "ROS gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.60,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.47,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 6.007e-15", x = 1, y = 2.5, size =8, colour = "black",face = "italic")


############################# SASP  ###################3
library(dplyr)
d8$Group <- paste0(d8$Diets,"-",d8$Celltype)
d <- d8  %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup() 
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd8446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets
ddd$Diets <- factor(ddd$Diets,levels = c("M9-BRD","M9-PRD"))
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "SASP gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.54,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.56,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 0.04715", x = 1, y = 2.5, size =8, colour = "black",face = "italic")

############################# DNA repair ###################3
library(dplyr)
d9$Group <- paste0(d9$Diets,"-",d9$Celltype)
d <- d9 %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd9446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "DNA repair gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.55,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 3.141e-07", x = 1, y = 2.5, size =8, colour = "black",face = "italic")

############################# Fibrosis###################3
library(dplyr)
head(d10)
d10$Group <- paste0(d10$Diets,"-",d10$Celltype)
d <- d10  %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC","B"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd10446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,scales = "free_y",ncol =6)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Fibrosis gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.82,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.90,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 0.0001503", x = 1, y = 2.5, size =8, colour = "black",face = "italic")


############################# Lipid ###################
library(dplyr)
d11$Group <- paste0(d11$Diets,"-",d11$Celltype)
d <- d11  %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd11446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets <- factor(ddd$Diets,levels = c("M9-BRD","M9-PRD"))

mean(ddd[ddd$Diets=="M9-PRD",]$Scores)
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Lipid storage gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.44,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.55,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 3.745e-07", x = 1, y = 2.5, size =8, colour = "black",face = "italic")
#F9837B


############################# Senescense ###################
library(dplyr)
d12$Group <- paste0(d12$Diets,"-",d12$Celltype)
d <- d12  %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC","B","SC","EC","T"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd11446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Fibrosis gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.55,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 3.141e-07", x = 1, y = 2.5, size =8, colour = "black",face = "italic")

############################# UPR ###################
library(dplyr)
d13$Group <- paste0(d13$Diets,"-",d13$Celltype)
d <- d13  %>% 
  dplyr::group_by(Group) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd11446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets

mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Fibrosis gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.55,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept = 0.53,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 3.141e-07", x = 1, y = 2.5, size =8, colour = "black",face = "italic")

M9
DotPlot(object = M9,group.by = "diet",
        cols = c("#4B549B","#E91C22"), 
        features = unique(Lipid_biosynthesis))+
  theme_classic()+
  #  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
  #                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.position = "top",
        legend.title = element_text(size = 13))+xlab("SASP gene set")





############################# Functions ###################
library(dplyr)
d14$Group <- paste0(d14$Diets,"-",d14$Celltype)
d <- d14  %>% 
  dplyr::group_by(Diets) %>%
  mutate(value = remove_outliers(Scores)) %>% 
  ungroup()  
#d <- d[d$Celltype%in%c("GC","EDC","DC"),]
ggplot(d, aes(x = Diets, y =Scores,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#Fd11446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Cell cycle 
gene-set score")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16))+facet_wrap(~Celltype,ncol = 7)
#d <- subset(d,Scores>0)
ddd <- d[d$Diets%in%c("M9-BRD","M9-PRD"),]
ddd <- subset(ddd,Scores>0)
ddd$Diets
ddd$Diets <- factor(ddd$Diets,levels = c("M9-BRD","M9-PRD"))
mean(ddd[ddd$Diets=="M9-BRD",]$Scores)
mean(ddd[ddd$Diets=="M9-PRD",]$Scores)

wilcox.test(ddd[ddd$Diets=="M9-BRD",]$Scores,ddd[ddd$Diets=="M9-PRD",]$Scores)

median(ddd[ddd$Diets=="M9-BRD",]$Scores)
median(ddd[ddd$Diets=="M9-PRD",]$Scores)

ggplot(ddd, aes(Scores)) +
  geom_density(aes(fill = factor(Diets)), alpha =0.9) +
  labs(
    title = "",
    subtitle = "",
    caption = "",
    x = "Steroid Hormone Synthesis/
Folliculogenesis gene-set score",
    y="Density",
    fill = "Class"
  )+theme_classic()+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),
        legend.text = element_text(size=14),legend.position = "right",legend.title =  element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+geom_vline(xintercept = 0.76,linetype="dashed",color="#F8766D")+
  geom_vline(xintercept =0.59 ,linetype="dashed",color="#00BFC4")+
  annotate("text", label = "p = 0.0002117", x = 1, y = 2.5, size =8, colour = "black",face = "italic")






############################### Cell interection #######################
###############GC.OO############
sigmean_Y=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/Y_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
sigmean_Y <- sigmean_Y[sigmean_Y$secreted=="True",]
LGC_Y <- sigmean_Y[sigmean_Y$receptor_a=="False",]
a <- na.omit(LGC_Y[,c('gene_a',"GC.OO")])
colnames(a) <- c('ligand',"level")
a$class <- "GC.OO"
a
LGC_Y <- sigmean_Y[sigmean_Y$receptor_b=="False",]
b <- na.omit(LGC_Y[,c('gene_a',"OO.GC")])
colnames(b) <- c('ligand',"level")
b$class <- "OO.GC"
head(b)
a
GCOO_Y <- rbind(a,b)
dim(GCOO_Y)
GCOO_Y

sigmean_O=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/O_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
sigmean_O <- sigmean_O[sigmean_O$secreted=="True",]
LGC_O <- sigmean_O[sigmean_O$receptor_a=="False",]
a <- na.omit(LGC_O[,c('gene_a',"GC.OO")])
colnames(a) <- c('ligand',"level")
a$class <- "GC.OO"

LGC_O <- sigmean_O[sigmean_O$receptor_b=="False",]
b <- na.omit(LGC_O[,c('gene_a',"OO.GC")])
colnames(b) <- c('ligand',"level")
b$class <- "OO.GC"
head(b)
a
GCOO_O <- rbind(a,b)
dim(GCOO_O)
GCOO_O
setdiff(GCOO_Y$ligand,GCOO_O$ligand)

list <- setdiff(GCOO_O$ligand,GCOO_Y$ligand)%>%c()
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")

BP_GCOO <- clusterProfiler::simplify(enrichGO(gene=list,keyType = "SYMBOL",
                                         OrgDb= org.Hs.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

View(BP_GCOO@result)
KEGG_GCOO <- enrichKEGG(gene =sym$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff =0.05)
#View(KEGG_GCOO@result)
write.table(KEGG_GCOO@result, 'KEGG_GC_OO_pairs.txt', sep='\t', quote=F)

KEGG_GC_OO_pairs <- read.delim("./KEGG_GC_OO_pairs.txt")

s <- KEGG_GC_OO_pairs[KEGG_GC_OO_pairs$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
dim(s)
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for GC-OO LRs") +
  geom_text(aes(x = rep(0.002,7), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)


MtoH <- readRDS("MtoH.rds")
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
#strsplit2(s$geneID,"/")%>%c()
#sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL
list_M <- MtoH[MtoH$HGNC.symbol%in%list,]$MGI.symbol
#list_M <- MtoH[MtoH$HGNC.symbol%in%sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL,]$MGI.symbol

#View(s)

GC <- readRDS("./GC_BP.rds")
GC$abbreviation
data <- as.data.frame(GC@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M6-BRD","M9-PRD","M9-BRD"))
head(ggdata)
#data <- ggdata[ggdata$SYMBOL%in%c("Pdgfb","Col4a2","Ncam1","Fgfr2","Kitl","Efna4", "Tyro3","Axl","Nectin1","Fgfr1","Spp1","Vegfa","Efna1","Notch2"),]
data <- ggdata[ggdata$SYMBOL%in%list_M,]
# DoHeatmap(GC1,features =list_M,group.by = "Diets",slot = "scale.data")+NoLegend()
DotPlot(GC,features = list_M,group.by = "Diets",cols = c("#555397","#E42228"))+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  ggtitle("Ligands: GC to OO")+scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))

#data <- ggdata[ggdata$SYMBOL%in%list,]
head(data)
d <- data %>% 
  dplyr::group_by(group) %>%
  mutate(value = remove_outliers(value)) %>% 
  ungroup()


ggplot(data, aes(x = group, y =log(value+1),fill=group)) +
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+ 
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#8A9197","#D2AF81"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 5,
    tip_length = 0,vjust = 0.01,test = "t.test")+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 18))+
  ggtitle("Ligands: GC to OO")+scale_x_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))

mean(data[data$group=="M9-BRD",]$value)
mean(data[data$group=="M6-PRD",]$value)


mean_GC <- apply(as.data.frame(GC@assays$RNA@data),1,function(x){tapply(x,paste0(GC$Diets,"-","GC"),mean)}) %>% t() %>% as.data.frame()
head(mean_GC)
dim(mean_GC)
colnames(mean_GC)=strsplit2(colnames(mean_GC),"_")[,1]

data1 <- mean_GC[rownames(mean_GC)%in%c("Cdkn2a","Cdkn2d","Cdkn1b","Cdkn1a"),]
data1

GC$Diets
SC$g <- paste0(SC$Diets,"-","SC")

mean_SC <- AverageExpression(SC,group.by = "g")

head(mean_SC)
dim(mean_SC)
colnames(mean_SC)=strsplit2(colnames(mean_SC),"_")[,1]

data2 <- mean_SC$SCT[rownames(mean_SC$SCT)%in%c("Cdkn2a","Cdkn2d","Cdkn1b","Cdkn1a"),]
data3 <- cbind(data1,data2)
n <- t(scale(t(data2)))
#n[n>1]=1
#n[n<1]=-1
rownames(data)
pheatmap::pheatmap(t(n),fontsize = 14,border_color = "black",
                   cluster_cols = F,cluster_rows = F)





###############SC.OO############
sigmean_Y=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/Y_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
#View(sigmean_Y)
a <- na.omit(sigmean_Y[,c('gene_a',"gene_b","SC.OO")])
colnames(a) <- c('gene_a',"gene_b","level")
a$class <- "SC.OO"
b <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.SC")])
colnames(b) <- c('gene_a',"gene_b","level")
b$class <- "OO.SC"
dim(b)
SCOO_Y <- rbind(a,b)
dim(SCOO_Y)
SCOO_Y$pairs <- paste0(SCOO_Y$gene_a,"_",SCOO_Y$gene_b)

sigmean_O=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/O_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
c <- na.omit(sigmean_O[,c('gene_a',"gene_b","SC.OO")])
colnames(c) <- c('gene_a',"gene_b","level")
c$class <- "SC.OO"
dim(c)
d <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.SC")])
colnames(d) <- c('gene_a',"gene_b","level")
d$class <- "OO.SC"
dim(d)

SCOO_O <- rbind(c,d)
dim(SCOO_O)

SCOO_O$pairs <- paste0(SCOO_O$gene_a,"_",SCOO_O$gene_b)
list <- strsplit2(setdiff(SCOO_Y$pairs,SCOO_O$pairs),"_")%>%c()
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
BP_SCOO <- clusterProfiler::simplify(enrichGO(gene=list,keyType = "SYMBOL",
                                         OrgDb= org.Hs.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

View(BP_SCOO@result)
KEGG_SCOO <- enrichKEGG(gene =sym$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff =0.05)
#View(KEGG_SCOO@result)
write.table(KEGG_SCOO@result, 'KEGG_SC_OO_pairs.txt', sep='\t', quote=F)
KEGG_SC_OO_pairs <- read.delim("./KEGG_SC_OO_pairs.txt", row.names=1)
s <- KEGG_SC_OO_pairs[KEGG_SC_OO_pairs$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for SC-OO LRs") +
  geom_text(aes(x = rep(0.002,8), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)
s[s$Description=="Efferocytosis",]


MtoH <- readRDS("MtoH.rds")
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
strsplit2(s$geneID,"/")%>%c()
sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL

list_M_SC <- MtoH[MtoH$HGNC.symbol%in%strsplit2(setdiff(SCOO_Y$pairs,paste0(SCOO_O$gene_a,"_",SCOO_Y$gene_b)),"_")[,1],]$MGI.symbol
#list_M <- MtoH[MtoH$HGNC.symbol%in%sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL,]$MGI.symbol

View(s)

SC1$Diets
SC <- readRDS("./SC_BP.rds")
SC$abbreviation2
data <- as.data.frame(SC@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M6-BRD","M9-PRD","M9-BRD"))

head(ggdata)
data <- ggdata[ggdata$SYMBOL%in%list_M,]
data$group <- paste0(substr(data$variable,1,4),"RD")
# DoHeatmap(SC1,features =list_M,group.by = "Diets",slot = "scale.data")+NoLegend()
DotPlot(SC,features = list_M,group.by = "Diets",cols = c("white","#E42228"))+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+ggtitle("Ligands from SC to OO")
data$group
d <- data %>% 
  dplyr::group_by(group) %>%
  mutate(value = remove_outliers(value)) %>% 
  ungroup()
d <- subset(data,value>0)
ggplot(na.omit(d), aes(x = group, y =value,fill=group))+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+ 
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-PRD","M6-BRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 5,
    tip_length = 0,vjust = 0.1)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))+ggtitle("Ligands from SC to OO")

###############M.OO############

sigmean_Y=read.table("/home/tuyx/MRNA_aging/cellphonedb-data-4.1.0/Y_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
#View(sigmean_Y)
a <- na.omit(sigmean_Y[,c('gene_a',"gene_b","M.OO")])
colnames(a) <- c('gene_a',"gene_b","level")
a$class <- "M.OO"
b <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.M")])
colnames(b) <- c('gene_a',"gene_b","level")
b$class <- "OO.M"
dim(b)
MOO_Y <- rbind(a,b)
dim(MOO_Y)
MOO_Y$pairs <- paste0(MOO_Y$gene_a,"_",MOO_Y$gene_b)

sigmean_O=read.table("/home/tuyx/MRNA_aging/cellphonedb-data-4.1.0/O_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
c <- na.omit(sigmean_O[,c('gene_a',"gene_b","M.OO")])
colnames(c) <- c('gene_a',"gene_b","level")
c$class <- "M.OO"
dim(c)
d <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.M")])
colnames(d) <- c('gene_a',"gene_b","level")
d$class <- "OO.M"
dim(d)

MOO_O <- rbind(c,d)
dim(MOO_O)

MOO_O$pairs <- paste0(MOO_O$gene_a,"_",MOO_O$gene_b)
list <- strsplit2(setdiff(MOO_Y$pairs,MOO_O$pairs),"_")%>%c()
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
#BP_MOO <- clusterProfiler::simplify(enrichGO(gene=list,keyType = "SYMBOL",
#                                              OrgDb= org.Hs.eg.db, ont="BP",
#                                              pAdjustMethod="BH",pvalueCutoff=0.05,
#                                              qvalueCutoff=0.05))

#View(BP_MOO@result)
KEGG_MOO <- enrichKEGG(gene =sym$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff =0.05)
#View(KEGG_MOO@result)
write.table(KEGG_MOO@result, 'KEGG_M_OO_pairs.txt', sep='\t', quote=F)

KEGG_M_OO_pairs <- read.delim("./KEGG_M_OO_pairs.txt", row.names=1)
s <- KEGG_M_OO_pairs[KEGG_M_OO_pairs$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for M-OO LRs") +
  geom_text(aes(x = rep(0.002,8), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)
s[s$Description=="Efferocytosis",]
######################### GC subcell ##############################
gc()
GC <- readRDS("./GC_BP.rds")
GC <- FindNeighbors(GC, dims = 1:30)
GC <- FindClusters(GC, resolution = 0.4)
GC <- RunUMAP(GC, dims = 1:30)
GC <- RunTSNE(GC, dims = 1:30)
ImmGen.se=ImmGenData() #(鼠)
Mouse.se=MouseRNAseqData() #(鼠)

GC_SingleR <- GetAssayData(GC, slot="data") ##获取标准化矩阵
GC.hesc <- SingleR(test = GC_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
GC.hesc2 <- SingleR(test = GC_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
GC@meta.data$SingleR_Mouse <- GC.hesc$labels
GC@meta.data$SingleR_ImmGen <- GC.hesc2$labels
DimPlot(GC, reduction = "umap", pt.size = 1, label = T)

DimPlot(GC, reduction = "umap", pt.size = 1, label = T,group.by = "seurat_clusters")

# 3 5 Mitotic GC
# 2 6 Preantral GC
# 4 7 Atretic GC
# 1 Antral GC
DotPlot(GC,features = c(#"Cyp11a1",#Luteal 7
  "Top2a","Cenpa","Pcna",# Mitotic GCs 5
  "Inhbb","Fst","Gja1",# Antral GCs 34
  "Igfbp5","Gatm","Col18a1", # Preantral GCs 1
  "Inhba","Bx4","Nr5a2",
  "Foxl2","Cyp19a1","Fshr",
  "Itih5","Cald1","Pik3ip1" # 'Atretic GC' 027
),cols = c("#584EA1","#ED1C26"),group.by = "seurat_clusters")+theme(axis.text.x = element_text(angle = 45,hjust =1))

DotPlot(GC,features = c(#"Cyp11a1",#Luteal 7
                        "Top2a","Cenpa","Pcna",# Mitotic GCs 5
                        "Inhbb","Fst","Gja1",# Antral GCs 34
                        "Igfbp5","Gatm","Col18a1", # Preantral GCs 1
                        "Inhba","Bx4","Nr5a2",
                        "Foxl2","Cyp19a1","Fshr",
                        "Itih5","Cald1","Pik3ip1" # 'Atretic GC' 027
                        ),cols = c("#584EA1","#ED1C26"))+theme(axis.text.x = element_text(angle = 45,hjust =1))

DotPlot(GC,features = c("Top2a","Cenpa","Pcna",
                        "Inhbb","Fst","Gja1",
                        "Igfbp5","Gatm","Col18a1",
                        "Inhba","Bx4","Nr5a2",
                        "Foxl2","Cyp19a1","Fshr",
                        "Itih5","Cald1","Pik3ip1"),cols = c("white","#584EA1"))+theme(axis.text.x = element_text(angle = 45,hjust =1))


GC1 <- subset(GC,idents=c(0:5,7))
new.cluster.ids <- c('Atretic GC','Preantral GC','Atretic GC','Antral GC','Antral GC','Mitotic GC','Atretic GC')
names(x=new.cluster.ids)=levels(x = GC1)
GC1 <- RenameIdents(object =GC1, new.cluster.ids)
GC1 <- readRDS("./GC1_BP.rds")
DotPlot(GC1,features = c("Igfbp5","Gatm","Col18a1", # Preantral GCs 1
                         "Top2a","Cenpa","Pcna",# Mitotic GCs 5
                         "Inhbb","Fst","Gja1",# Antral GCs 34
                         
                         "Inhba","Bx4","Nr5a2",
                         "Foxl2","Cyp19a1","Fshr",
                         "Itih5","Cald1","Pik3ip1"),
        cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))+
  scale_y_discrete(limits=rev(c("Preantral GC","Mitotic GC","Antral GC","Atretic GC")))

#INSL3, APOE, GSTA1, APOA1, FDX1 and CYP17A1 
DotPlot(GC1,features = c("Insl3","Apoe","Gsta1","Apoa1","Fdx1","Cyp17a1",
                         "Amh","Fst","Hsd17b1","Serpine2","Prkar2b",
                         "Dcn","Lgals1","Lgals3","Ptgfr"),
        cols = c("#584EA1","#ED1C26"),group.by = "celltype")+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))

gc <- MtoH[MtoH$HGNC.symbol%in%c("P3H4","CCDC36",
                                 "FMN2","HFM1","M1AP","MEI4","MEIOB","SPANXA2","TERF1","TEX15",
                                 "DMC1","DMRTC2","ING2","LEMD2","PTTG3P","SPATA22","SYCE3",
                                 "ATM","CCDC156","CCNE1","ESPL1","FBXO5","MEI1","MSH5","PLK1",
                                 "RAD51D","SLX4","STRA8","SYCE1L","SYCE2","UBR2","WEE2",
                                 "BRCA2","BUB1B","CENPC","EME1","ERCC4","FANCD2","FANCM","MASTL",
                                 "MEIKIN","MND1","MSH2","RAD21"),]$MGI.symbol

DotPlot(GC1,features = gc)+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))



# INSL3, APOE, GSTA1, APOA1, FDX1 and CYP17A1  Hormone synthesis
# AMH, FST, HSD17B1, SERPINE2 and PRKAR2B Follicular development
# DCN, LGALS1 and LGALS3
GC1 <- RunTSNE(GC1, dims = 1:20)
GC1 <- RunUMAP(GC1, dims = 1:30)

#saveRDS(GC1,"./GC1_BP.rds")

GC1 <- readRDS("./GC1_BP.rds")
GC1@active.ident <- factor(GC1@active.ident,
                           levels = c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))

DimPlot(GC1, reduction = "tsne", pt.size = 1, label = T,group.by = "celltype")+scale_color_brewer()
DimPlot(GC1, reduction = "umap", pt.size = 1, label = T,group.by = "celltype")+scale_color_brewer()

saveRDS(GC1,"./GC1_BP.rds")
GC1$

GC1$stage <- substr(GC1$Diets,1,2)
df1=as.data.frame(GC1@active.ident)
GC1@meta.data$celltype=df1$`GC1@active.ident`
GC1@meta.data$stage_cell <- paste0(GC1$stage,"-",GC1@meta.data$celltype)

df1=as.data.frame(GC1@active.ident)
idents <- paste0(GC1$stage,"-",df1$`GC1@active.ident`)
names(idents) <- rownames(GC1@meta.data)
GC1@active.ident<-factor(idents)
GC1$group <- paste0(substr(GC1$Diets,1,6),"-",strsplit2(GC1@active.ident,"-")[,2])
GC1@meta.data$group

Mito<- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M6-BRD-Mitotic GC",group.by = "group", min.pct = 0.25)
Pre<- FindMarkers(GC1, ident.1 = "M9-BRD-Preantral GC", ident.2 = "M6-BRD-Preantral GC",group.by = "group", min.pct = 0.25)
Ant<- FindMarkers(GC1, ident.1 = "M9-BRD-Antral GC", ident.2 = "M6-BRD-Antral GC",group.by = "group", min.pct = 0.25)
Atr<- FindMarkers(GC1, ident.1 = "M9-BRD-Atretic GC", ident.2 = "M6-BRD-Atretic GC",group.by = "group", min.pct = 0.25)

data1=data.frame(cell=c("Mit GC","Pre GC","Ant GC","Atr GC"),diff=c(dim(Mito)[1],
                                                                    dim(Pre)[1],
                                                                    dim(Ant)[1],
                                                                    dim(Atr)[1]))
data1$diets <- "BRD"
Ant2<- FindMarkers(GC1, ident.1 = "M9-PRD-Antral GC", ident.2 = "M9-BRD-Antral GC",group.by = "group", min.pct = 0.25)

################## PRD ##########################
Mito3<-FindMarkers(GC1, ident.1 = "M9-PRD-Mitotic GC", ident.2 = "M6-PRD-Mitotic GC",group.by = "group", min.pct = 0.25)
Pre3<- FindMarkers(GC1, ident.1 = "M9-PRD-Preantral GC", ident.2 = "M6-PRD-Preantral GC",group.by = "group", min.pct = 0.25)
Ant3<- FindMarkers(GC1, ident.1 = "M9-PRD-Antral GC", ident.2 = "M6-PRD-Antral GC",group.by = "group", min.pct = 0.25)
Atr3<- FindMarkers(GC1, ident.1 = "M9-PRD-Atretic GC", ident.2 = "M6-PRD-Atretic GC",group.by = "group", min.pct = 0.25)
dim(Ant2)
dim(Ant3)
Ant_resistance <- intersect(rownames(Ant2[Ant2$avg_log2FC>0,]),rownames(Ant3[Ant3$avg_log2FC>0,]))
Ant_rescue <- intersect(rownames(Ant2[Ant2$avg_log2FC<0,]),rownames(Ant3[Ant3$avg_log2FC<0,]))
which(Ant_resistance=="Foxo3")
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Ant_resistance,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Ant_resistance <- enrichKEGG(gene =ncbi$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Ant_rescue,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Ant_rescue <- enrichKEGG(gene =ncbi$ENTREZID,
                                  organism = 'mmu',
                                  pvalueCutoff =0.05)

write.csv(KEGG_Ant_rescue@result,"./KEGG_Ant_rescue.csv")
write.csv(KEGG_Ant_resistance@result,"./KEGG_Ant_resistance.csv")
KEGG_Ant_rescue <- read.csv("./KEGG_Ant_rescue.csv")
KEGG_Ant_resistance <- read.csv("./KEGG_Ant_resistance.csv")
KEGG_Ant_rescue <- KEGG_Ant_rescue[KEGG_Ant_rescue$select==1,]%>%na.omit()
KEGG_Ant_resistance <- KEGG_Ant_resistance [KEGG_Ant_resistance$select==1,]%>%na.omit()

KEGG_Ant_rescue$cell <- "Rescue"
KEGG_Ant_resistance$cell <- "Resistance"
KEGG_Ant_select <- rbind(KEGG_Ant_rescue,KEGG_Ant_resistance)
b <- lapply(str_split(KEGG_Ant_select$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
KEGG_Ant_select$GeneRatio <- b
ggplot(KEGG_Ant_select,aes(cell,Description,size = GeneRatio))+
  geom_point(shape=21,aes(fill= pvalue),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,8))+
  theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle=45, hjust=1),
        axis.text.y  = element_text(color="black",size=16),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black")
  ) +scale_x_discrete(limits=c("Resistance","Rescue"))+
  scale_y_discrete(limits= unique(KEGG_Ant_select$Description)) 



a=KEGG_Ant_resistance[KEGG_Ant_resistance$Description%in%c("DNA replication","Mismatch repair"),]$geneID%>%str_split("/")%>%unlist()


b=KEGG_Ant_rescue$geneID%>%str_split("/")%>%unlist()
b=KEGG_Ant_rescue[KEGG_Ant_rescue$Description%in%c("FoxO signaling pathway"),]$geneID%>%str_split("/")%>%unlist()

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=unique(a),columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")
ncbi2<- AnnotationDbi::select(org.Mm.eg.db,keys=unique(b),columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")

ncbi$SYMBOL
GC1@active.ident
#Ant <- subset(GC1,idents=c("M6-Antral GC","M9-Antral GC"))
n <- t(scale(t(as.data.frame(Ant@assays$RNA@data))))
n <- na.omit(n)
mean_Ant <- apply(as.data.frame(n),1,function(x){tapply(x,Ant$Diets,mean)}) %>% t() %>% as.data.frame()

mean_Ant_s <- mean_Ant[ncbi$SYMBOL,c(3,4,6)]
mean_Ant_s2 <- mean_Ant[ncbi2$SYMBOL,c(3,4,6)]

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))
dim(mean_Ant_s)
h1 <- Heatmap(t(mean_Ant_s),name="scale",row_dend_reorder = F,row_title = "Antral GC",
               col=col_fun,column_title = "Resistance",show_column_names  = T)
h1
dim(mean_Ant_s2)
h2 <- Heatmap(t(mean_Ant_s2),name="scale",row_dend_reorder = T,row_title = "Antral GC",
              col=col_fun,column_title = "Rescue",show_column_names  = T)
h1+h2








data3=data.frame(cell=c("Mit GC","Pre GC","Ant GC","Atr GC"),diff=c(dim(Mito3)[1],
                                                                    dim(Pre3)[1],
                                                                    dim(Ant3)[1],
                                                                    dim(Atr3)[1]))
data3$diets <- "PRD"

data <- rbind(data1,data3)
data$diets <- factor(data$diets,levels = c("BRD","WRD","PRD"))
data$cell <- factor(data$cell,levels = c("Mit GC", "Pre GC", "Ant GC", "Atr GC"))
ggplot(data, aes(x =cell,y = diff))+
  geom_bar(stat="identity",fill="#6BAED6",color="black",width = 0.8)+
  facet_wrap(~diets)+
  theme_classic()+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle =45,hjust=1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(color = 'black'),axis.line.y=element_line(color = 'black'))+theme(legend.position = 'none')+
  ylab("Differentially expressed
genes (DEGs)")+
  theme(strip.text.x = element_text(size=14))

An_GC

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Mito3[Mito3$avg_log2FC>0,]%>%rownames(),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_mito_O <- enrichKEGG(gene =ncbi$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)
KEGG_mito_O@result <- subset(KEGG_mito_O@result,pvalue<0.05)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Pre3[Pre3$avg_log2FC>0,]%>%rownames(),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_pre_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
KEGG_pre_O@result <- subset(KEGG_pre_O@result,pvalue<0.05)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Ant3[Ant3$avg_log2FC>0,]%>%rownames(),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Ant_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
KEGG_Ant_O@result <- subset(KEGG_Ant_O@result,pvalue<0.05)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=Atr3[Atr3$avg_log2FC>0,]%>%rownames(),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Atr_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
KEGG_Atr_O@result <- subset(KEGG_Atr_O@result,pvalue<0.05)
View(KEGG_Ant_O@result)
write.csv(KEGG_Ant_O@result,"./KEGG2_Ant_O.csv")
write.csv(KEGG_mito_O@result,"./KEGG2_Mito_O.csv")
write.csv(KEGG_pre_O@result,"./KEGG2_Pre_O.csv")
write.csv(KEGG_Atr_O@result,"./KEGG2_Atr_O.csv")

KEGG_Mito_O <- read.csv("~/scRNA_aging/output/KEGG2_Mito_O.csv", row.names=1)
KEGG_Pre_O <- read.csv("~/scRNA_aging/output/KEGG2_Pre_O.csv", row.names=1)
KEGG_Ant_O <- read.csv("~/scRNA_aging/output/KEGG2_Ant_O.csv", row.names=1)
KEGG_Atr_O <- read.csv("~/scRNA_aging/output/KEGG2_Atr_O.csv", row.names=1)

KEGG_Mito_O <- KEGG_Mito_O[KEGG_Mito_O$select==1,]%>%na.omit()
KEGG_Pre_O <- KEGG_Pre_O[KEGG_Pre_O$select==1,]%>%na.omit()
KEGG_Ant_O <- KEGG_Ant_O[KEGG_Ant_O$select==1,]%>%na.omit()
KEGG_Atr_O <- KEGG_Atr_O[KEGG_Atr_O$select==1,]%>%na.omit()

KEGG_Mito_O$class <- "Mito GC"
KEGG_Pre_O$class <- "Pre GC"
KEGG_Ant_O$class <- "Ant GC"
KEGG_Atr_O$class <- "Atr GC"

KEGG_select <- rbind(KEGG_Mito_O,KEGG_Pre_O)%>%rbind(.,KEGG_Ant_O)%>%rbind(.,KEGG_Atr_O)
KEGG_select
b <- lapply(str_split(KEGG_select$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
KEGG_select$GeneRatio <- b
ggplot(KEGG_select,aes(class,Description,size = GeneRatio))+
  geom_point(shape=21,aes(fill= pvalue),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,8))+
  theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle=45, hjust=1),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black")
  ) +scale_x_discrete(limits=c("Mito GC","Pre GC","Ant GC","Atr GC"))+
  scale_y_discrete(limits= unique(KEGG_select$Description)) 

nn <- KEGG_select[KEGG_select$Description=="Apoptosis",]$geneID%>%strsplit2(.,"/")%>%c()

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=nn,columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")
data <- as.data.frame(GC1@assays$RNA@data)
data$gene <- rownames(data)
ggdata <- melt(data)

ggdata$diet <- paste0(substr(ggdata$variable,1,4),"RD")

unique(ggdata$diet)
ggdata <- ggdata[ggdata$gene%in%ncbi$SYMBOL,]

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
data <- ggdata %>% 
  group_by(diet) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()

ggplot(ggdata,aes(x = diet,y = log2(value+1),fill=diet))+geom_violin(aes(fill = diet)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M9-BRD","M9-PRD"),c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0.1,test="t.test")+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16),
        strip.text.x = element_text(size=16),strip.background = element_blank())+ggtitle("Jund")


###################### subGC cell proportion ####################
library("scProportionTest")
GC1@meta.data$celltype=df1$`GC1@active.ident`


library("scProportionTest")
#cell <- GC1@meta.data[grep(7,GC1$seurat_clusters,invert = T),]%>%rownames(.)
#GC2 <- subset(GC1,cells=cell)

prop_test <- sc_utils(GC1)
prop_test2 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6-PRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test2,log2FD_threshold = log2(1.5),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))+ggtitle("PRD: 9m vs 6m")


prop_test3 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test3,log2FD_threshold = log2(1.5),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))+ggtitle("9m: PRD vs BRD")
GC1$stage_cell
CS<- FindMarkers(GC1, ident.1 = "M6-Preantral GC", ident.2 = "M6-Antral GC",group.by = "stage_cell", min.pct = 0.25)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=CS%>%rownames(),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_CS_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
View(KEGG_CS_O@result)
KEGG_CS_O@result <- subset(KEGG_CS_O@result,pvalue<0.05)
BP_CS<- clusterProfiler::simplify(enrichGO(gene=ncbi$SYMBOL,keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

View(BP_CS@result)

int <- intersect(ncbi$SYMBOL,c("Fshr","Lhr","Amh","Cyp19a1",
                        "Inha","Inhba","Bmp15","Gdf9",
                        "Ptgs2","Nr5a1","Nr5a2","Star","Hsd3b1",
                        "Hsd17b1","Hsd17b7"))

data <- as.data.frame(GC1@assays$RNA@data)
data$gene <- rownames(data)
ggdata <- melt(data)
ggdata$variable
ggdata$diet <- paste0(substr(ggdata$variable,1,4),"RD")
unique(ggdata$diet)
ggdata <- ggdata[ggdata$gene%in%c("Fshr","Lhr","Amh","Cyp19a1",
                                  "Inha","Inhba","Bmp15","Gdf9",
                                  "Ptgs2","Nr5a1","Nr5a2","Star","Hsd3b1",
                                  "Hsd17b1","Hsd17b7"),]

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}
data <- ggdata %>% 
  group_by(diet) %>% 
  mutate(value = remove_outliers(value)) %>% 
  ungroup()

ggplot(ggdata,aes(x = diet,y = log2(value+1),fill=diet))+geom_violin(aes(fill = diet)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M9-BRD","M9-PRD"),c("M6-BRD","M6-PRD"),c("M6-PRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0.1,test="t.test")+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16),
        strip.text.x = element_text(size=16),strip.background = element_blank())+ggtitle("Jund")

mean_GC <- apply(as.data.frame(GC1@assays$RNA@data),1,function(x){tapply(x,GC1$Diets,mean)}) %>% t() %>% as.data.frame()
head(mean_GC)
dim(mean_GC)
colnames(mean_GC)=strsplit2(colnames(mean_GC),"_")[,1]
data_GC <- as.data.frame(GC1@assays$SCT@data)[
  rownames(GC1)%in%c("Fshr","Inha","Nr5a1","Cyp19a1"),]
data_GC$SYMBOL <- rownames(data_GC)
ggdata <- reshape2::melt(data_GC)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M9-PRD","M9-BRD"))
ggdata$celltype <- "GC1"
head(ggdata)
ggplot(na.omit(ggdata), aes(x = group, y =value,fill=group)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.05,size = 0.5,textsize = 3,
              tip_length = 0,vjust = 0.1)+ylab("Cumulus Expansion Factor Expression")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ 
  theme(plot.title = element_text(size = 14),
        strip.text.x =element_text(size = 14),strip.background = element_blank())+
  facet_wrap(~SYMBOL,ncol = 2,scales = "free_y")












data1 <- mean_GC[rownames(mean_GC)%in%c("Fshr","Lhr","Amh","Cyp19a1",
                                       "Inha","Inhba","Bmp15","Gdf9",
                                       "Ptgs2","Nr5a1","Nr5a2","Star","Hsd3b1",
                                       "Hsd17b1","Hsd17b7"),c(1,3,4,6)]

data <- mean_GC[rownames(mean_GC)%in%c("Kitl","Amh","Cx43","Fshr","Inha","Gdf9","Bmp15","Kit","Foxl2",
                                       "Notch1","Egfr","Fgfr2","Igf1","Igf1r","Vcam1","Adm","Hbegf",
                                       "Cxcr4","S1pr1","Nlrp5"),c(1,3,4,6)]

data
n <- t(scale(t(data)))
n
#n[n>1]=1
#n[n<1]=-1
rownames(data)
pheatmap::pheatmap(t(n),fontsize = 14,border_color = "black",
                   cluster_cols = T,cluster_rows  = F)

########################## GC1 Monocle2 ###########################
devtools::load_all("~/monocle")

library(monocle)
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
GC1 <- readRDS("./GC1_BP.rds")
GC1$celltype <- GC1@active.ident

GC2 <- subset(GC1,idents=c("Mitotic GC", "Antral GC"))
#saveRDS(GC2,"./GC2.rds")

data <- as(as.matrix(GC1@assays$RNA@data), 'sparseMatrix')
pd <- GC1@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame",
          data=pd)
fd <- new("AnnotatedDataFrame",
          data=fData)
cd <- newCellDataSet(data,
                     phenoData = pd,
                     featureData =fd)#, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
rm(data)
gc()
sc_cds <- estimateSizeFactors(cd)
sc_cds <- estimateDispersions(sc_cds)
?dispersionTable
disp_table <- dispersionTable(sc_cds)
# first method for genes selected
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
sc_cds <- setOrderingFilter(sc_cds, unsup_clustering_genes$gene_id)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds$celltype <- pd$celltype
plot_cell_trajectory(sc_cds, color_by ="celltype",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

saveRDS(sc_cds,"GC21_sc_cds.rds")

diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = T)

plot_cell_trajectory(sc_cds, color_by ="celltype",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")

plot_cell_trajectory(sc_cds, color_by ="Pseudotime",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")

plot_cell_trajectory(sc_cds, color_by ="State",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")

saveRDS(sc_cds,"GC22_sc_cds.rds")
plot_cell_trajectory(sc_cds, color_by ="State",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")

sub <- data.frame(State=sc_cds$State)
rownames(sub) <- colnames(sc_cds)
sub$Cell <-  colnames(sc_cds)
sub1 <- sub[sub$State%in%c(1,2,3),]
GC2 <- subset(GC1,cells=sub1$Cell)
rm(GC1)
gc()

GC2$celltype <- GC2@active.ident

#GC2 <- subset(GC1,idents=c("Mitotic GC", "Antral GC"))
#saveRDS(GC2,"./GC2.rds")

data <- as(as.matrix(GC2@assays$RNA@data), 'sparseMatrix')
pd <- GC2@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame",
          data=pd)
fd <- new("AnnotatedDataFrame",
          data=fData)
cd <- newCellDataSet(data,
                     phenoData = pd,
                     featureData =fd)#, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
rm(data)
gc()
sc_cds <- estimateSizeFactors(cd)
sc_cds <- estimateDispersions(sc_cds)
?dispersionTable
disp_table <- dispersionTable(sc_cds)
# first method for genes selected
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
sc_cds <- setOrderingFilter(sc_cds, unsup_clustering_genes$gene_id)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds$celltype <- pd$celltype
plot_cell_trajectory(sc_cds, color_by ="celltype",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))



















################### Monocle 3 ##################
library(monocle3)
GC1$celltype <- GC1@active.ident
data <- as(as.matrix(GC1@assays$RNA@data), 'sparseMatrix')
pd <- GC1@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame",
          data=pd)
fd <- new("AnnotatedDataFrame",
          data=fData)
cd <- newCellDataSet(data,
                     phenoData = pd,
                     featureData =fd)#, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
rm(data)
gc()
sc_cds <- estimateSizeFactors(cd)
sc_cds <- estimateDispersions(sc_cds)

cds <- learn_graph(sc_cds)     #是的，就这么简单，用 learn_graph即可
#画图  我们可以根据meta信息来画出相关的图片
head(colData(cds))
rm(GC1)
gc()
plot_cells(sc_cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)                                    #这个图就可以看出细胞直接的分化轨迹了
#黑色的线显示的是graph的结构。
#数字带白色圆圈表示不同的结局，也就是叶子。
#数字带黑色圆圈代表分叉点，从这个点开始，细胞可以有多个结局。
#这些数字可以通过label_leaves和label_branch_points参数设置。







######################### GC1 细胞周期分析###############
library(Seurat)
str(cc.genes)
GC1 <- readRDS("./GC1_BP.rds")
GC1$celltype <- GC1@active.ident
#GC1 <- subset(x = GC1, subset = nFeature_RNA > 200 & nFeature_RNA<6000 & percent.mito<15)

MtoH <- readRDS("./MtoH.rds")
cc.genes <- readRDS("./cc.genes.rds")
s.genes <- MtoH[MtoH$HGNC.symbol%in%cc.genes$s.genes,]$MGI.symbol
g2m.genes <- MtoH[MtoH$HGNC.symbol%in%cc.genes$s.genes,]$MGI.symbol


GC1 <- CellCycleScoring(GC1,s.features = s.genes, 
                           g2m.features = g2m.genes, 
                           set.ident = TRUE)

GC1$Phase
?RidgePlot
GC1$celltype
cycle <- GC1@meta.data%>%dplyr::group_by(GC1@meta.data$celltype,GC1@meta.data$Diets)%>%
  dplyr::summarise(table(Phase))

head(cycle)
cycle$Phase <- names(cycle$`table(Phase)`)

colnames(cycle) <- c("Celltype","Group","Cell_Number","Phase")
cycle <- cycle[cycle$Group%in%c("M6-PRD","M9-PRD","M9-BRD"),]
cycle$Group <- factor(cycle$Group,levels = c("M6-PRD","M9-PRD","M9-BRD"))
cycle$Phase <- factor(cycle$Phase,levels = c("G1","S","G2M"))
cycle$Celltype <- factor(cycle$Celltype,
                                 levels = c("Preantral GC","Mitotic GC",
                                            "Antral GC","Atretic GC"))
ggplot(cycle) +				
  geom_bar(aes(x = Group,y=Cell_Number,fill = Phase),stat = "identity",position = "fill",colour= 'black') +						
  scale_fill_manual(values  = c("#6686c6","#a0c1db","gray","#E65F92","#c77364","#ce8f5c",						
                                "#7bac80","#75a5c0","#b5181a","#b72d9e",						
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f","#39cf40",						
                                "#3e6926","#0a0883","#49ddd0","#e0f8f7",						
                                "#651d91","#9d3eaf","#b9a3d0","#5b5a5b","#8f8f8f"))+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5))+guides(color=guide_legend(title = NULL))+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=14),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+	
  xlab("")+ylab("% of cell phase")+
  facet_wrap(~Celltype,ncol = 4,scales = "free_y")+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))+ 
  theme(
        strip.text.x =element_text(size = 14),strip.background = element_blank())						








GC1@meta.data$G2M.Score

GC1@meta.data$Diets <- factor(GC1@meta.data$Diets,levels = c("M6-PRD","M9-PRD","M9-BRD"))
GC1@meta.data$celltype <- factor(GC1@meta.data$celltype,
                                 levels = c("Mitotic GC","Preantral GC",
                                            "Antral GC","Atretic GC"))
ggplot(na.omit(GC1@meta.data[GC1$celltype%in%c( "Antral GC"),]), aes(x = Diets, y =S.Score,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.05,size = 0.5,textsize = 3,
              tip_length = 0,vjust = 0.1)+ylab("S Phase Scores")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  ggtitle("")+ 
  theme(plot.title = element_text(size = 14),
        strip.text.x =element_text(size = 14),strip.background = element_blank())+
  facet_wrap(~celltype,ncol = 4,scales = "free_y")

GC1$Foxo3 <- GC1@assays$SCT@data["Foxo3",]
GC1$Stat3 <- GC1@assays$SCT@data["Stat3",]
GC1$Jund <- GC1@assays$SCT@data["Jund",]
GC1$Foxp1 <- GC1@assays$SCT@data["Foxp1",]
GC1$Cdkn2a <- GC1@assays$SCT@data["Cdkn2a",]
GC1$Cdkn2d <- GC1@assays$SCT@data["Cdkn2d",]

ggplot(na.omit(GC1@meta.data[GC1$celltype%in%c( "Antral GC"),]), aes(x = Diets, y =Foxp1,fill=Diets)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.05,size = 0.5,textsize = 3,
              tip_length = 0,vjust = 0.1)+ylab("Foxp1 Expression Level")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  ggtitle("")+ 
  theme(plot.title = element_text(size = 14),
        strip.text.x =element_text(size = 14),strip.background = element_blank())+NoLegend()


ggplot(na.omit(GC1@meta.data[GC1$celltype%in%c( "Antral GC"),]), aes(x = Diets, y =log2(Cdkn2d+1),fill=Diets)) +
  stat_boxplot( width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.05,size = 0.5,textsize = 3,
              tip_length = 0,vjust = 0.1)+ylab("Cdkn2a Expression Level")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  ggtitle("")+ 
  theme(plot.title = element_text(size = 14),
        strip.text.x =element_text(size = 14),strip.background = element_blank())+NoLegend()


GC1@active.ident <- GC1$old.ident
Ant <- subset(GC1,idents="Antral GC")
cell <- Ant@meta.data[Ant$Diets%in%c("M6-PRD","M9-PRD","M9-BRD"),]
Ant <- subset(GC1,cells=rownames(cell))

s.genes
RidgePlot(Ant, group.by = "Diets",
          features = c("Pcna","Casp8ap2","Mcm6"), 
          cols = pal_npg("nrc", alpha = 0.7)(3),
          ncol = 4)


?group_by

RidgePlot(GC1,  
          cols = pal_npg("nrc", alpha = 0.7)(3),
          ncol = 2)

Ant_tf <- read.delim("~/scRNA_aging/TFs/output_Ant_recue/Step2_regulonTargetsInfo.tsv")
Ant_tf <- Ant_tf[Ant_tf$highConfAnnot=="TRUE",]
unique(Ant_tf$TF)
View(Ant_tf)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=c(Ant_tf$gene,Ant_tf$TF),
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Ant_tf <- enrichKEGG(gene =sym$ENTREZID,
                            organism = 'mmu',
                            pvalueCutoff =0.05)
View(KEGG_Ant_tf@result)
write.table(KEGG_Ant_tf@result, 'KEGG_Ant_tf.txt', sep='\t', quote=F)
KEGG_Ant_tf <- read.delim("./KEGG_Ant_tf.txt", row.names=1)

s <- KEGG_Ant_tf[KEGG_Ant_tf$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
s$Description
dim(s)
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd",direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "") +
  geom_text(aes(x = rep(0.002,9), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=14))+
  scale_y_discrete(limits=s$Description)







######################### SC subcell ##############################
SC()
SC <- readRDS("./SC_BP.rds")
SC <- FindNeighbors(SC, dims = 1:30)
SC <- FindClusters(SC, resolution = 0.4)
SC <- RunUMAP(SC, dims = 1:30)
SC <- RunTSNE(SC, dims = 1:30)
#ImmGen.se=ImmGenData() #(鼠)
#Mouse.se=MouseRNAseqData() #(鼠)
#SC_SingleR <- GetAssayData(SC, slot="data") ##获取标准化矩阵
#SC.hesc <- SingleR(test = SC_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
#SC.hesc2 <- SingleR(test = SC_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
#SC@meta.data$SingleR_Mouse <- SC.hesc$labels
#SC@meta.data$SingleR_ImmGen <- SC.hesc2$labels
DimPlot(SC, reduction = "umap", pt.size = 1, label = T)
DimPlot(SC, reduction = "umap", pt.size = 1, label = T,group.by = "seurat_clusters")

DotPlot(SC,features = c("Inhba","Bex4","Nr5a2",# GC 2  12
  "Dcn","Mgp","Lum",# Fibroblast-like 0
  "Cyp11a1","Mgarp","Cyp17a1", # Theca 5 11
  "Enpep","Ptch1","Hhip",# Early Theca 6 7 10
  "Rgs5","Notch3","Ebf1",# Pericytes 3
  "Cnn1","Myh11","Actg2" # Smooth Muscle 4
  ),cols = c("#584EA1","#ED1C26"),group.by = "seurat_clusters")+theme(axis.text.x = element_text(angle = 45,hjust =1))

DotPlot(SC,features = c("Inhba","Bex4","Nr5a2",# GC 2  12
                        "Dcn","Mgp","Lum",# Fibroblast-like 0 9
                        "Cyp11a1","Mgarp","Cyp17a1", # Theca 5 11
                        "Enpep","Ptch1","Hhip",# Early Theca 6  10
                        "Rgs5","Notch3","Ebf1",# Pericytes 3
                        "Cnn1","Myh11","Actg2", # Smooth Muscle 4
                        "Cd74"
),group.by = "seurat_clusters")+theme(axis.text.x = element_text(angle = 45,hjust =1))
SC <- subset(SC,idents=c(0,2,3,4,5,6,9,10,11,12))

DotPlot(SC,features = c("Inhba","Bex4","Nr5a2",# GC 2  12
                        "Dcn","Mgp","Lum",# Fibroblast-like 0 9
                        "Cyp11a1","Mgarp","Cyp17a1", # Theca 5 11
                        "Enpep","Ptch1","Hhip",# Early Theca 6  10
                        "Rgs5","Notch3","Ebf1",# Pericytes 3
                        "Cnn1","Myh11","Actg2"),cols = c("white","#584EA1"))+theme(axis.text.x = element_text(angle = 45,hjust =1))

new.cluster.ids <- c('Fibroblast-like','GC','Pericytes','Smooth Muscle','Theca',
                     'Early Theca','Fibroblast-like',"Early Theca","Theca","GC")
names(x=new.cluster.ids)=levels(x = SC)
SC1 <- RenameIdents(object =SC, new.cluster.ids)
SC1 <- subset(SC1,idents=c("Fibroblast-like","Theca",
                           "Early Theca","Pericytes","Smooth Muscle"))
DotPlot(SC1,features = c(#"Inhba","Bex4","Nr5a2",# GC 2 8 12
                         "Dcn","Mgp","Lum",# Fibroblast-like 0 9
                         "Cyp11a1","Mgarp","Cyp17a1", # Theca 5 11
                         "Enpep","Ptch1","Hhip",# Early Theca 6  10
                         "Rgs5","Notch3","Ebf1",# Pericytes 3
                         "Cnn1","Myh11","Actg2"),cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))+
  scale_y_discrete(limits=rev(c("Fibroblast-like","Theca",
                                "Early Theca","Pericytes","Smooth Muscle")))


#SC1 <- RunTSNE(SC1, dims = 1:22)
SC1 <- RunUMAP(SC1, dims = 1:30)
saveRDS(SC1,"./SC1_ <- .rds")
SC1 <- readRDS("./SC1_BP.rds")
SC1@active.ident <- factor(SC1@active.ident,
                           levels = c("Fibroblast-like","Theca",
                                      "Early Theca","Pericytes","Smooth Muscle"))
DimPlot(SC1, reduction = "umap", pt.size = 1, label = F)+scale_color_futurama()
#saveRDS(SC1,"./SC1_BP.rds")

SC1$celltype <- factor(SC1@active.ident)
prop_test <- sc_utils(SC1)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6-PRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(title  =element_text(size = 14),axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 15),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Fibroblast-like","Theca",
                            "Early Theca","Pericytes","Smooth Muscle"))+
  ggtitle("PRD: 9m vs 6m")+xlab("Celltypes")+ylab("log2FD")
SC1$diet
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "BRD", sample_2 = "PRD",
  sample_identity = "diet"
)
permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(title  =element_text(size = 14),axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 15),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Fibroblast-like","Theca",
                            "Early Theca","Pericytes","Smooth Muscle"))+
  ggtitle("9m: PRD vs BRD")+xlab("Celltypes")+ylab("log2FD")

permutation_plot2 <- function(
    sc_utils_obj,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1),
    order_clusters = TRUE
) {
  
  ## Retrieve results.
  plot_data <- data.table::copy(sc_utils_obj@results$permutation)
  
  ## Mark the significant results.
  plot_data[, significance := ifelse(
    FDR < FDR_threshold & abs(obs_log2FD) > log2FD_threshold,
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  )]
  
  plot_data[, significance := factor(significance, levels = c(
    paste("FDR <", FDR_threshold, "& abs(Log2FD) >", round(log2FD_threshold, 2)),
    "n.s."
  ))]
  
  ## Order the clusters by observed log2FD if requested.
  if (order_clusters) {
    plot_data[, clusters := fct_reorder(factor(clusters), dplyr::desc(obs_log2FD))]
  }
  
  ## Plot the results.
  p <- ggplot(plot_data, aes(x = clusters, y = obs_log2FD)) +
    geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
    theme_bw()  +
    geom_hline(yintercept = 0.58, lty = 2) +
    geom_hline(yintercept = -0.58, lty = 2) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("salmon", "grey")) +
    coord_flip()+theme(legend.position = "top")
  
  return(p)
}
SC1 <- readRDS("./SC1_BP.rds")
data <- as.data.frame(SC1@assays$SCT@data)[
  rownames(SC1)%in%c("Col1a1","Col1a2","Col3a1","Col4a1"),]

data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
rm(SC1,data)
gc()
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M9-PRD","M9-BRD"))
ggdata$celltype <- "SC1"
head(ggdata)
ggplot(na.omit(ggdata), aes(x = group, y =value,fill=group)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
    tip_length = 0,vjust = 0.1)+ylab("Collegen Expression")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ 
  theme(plot.title = element_text(size = 16),strip.text.x =element_text(size = 16),strip.background =element_blank())+
  facet_wrap(~SYMBOL,ncol = 7,scales = "free_y")

Fib <- subset(SC1,idents=c("Fibroblast-like"))
Fib$Diets
rm(SC1)
gc()
saveRDS(Fib,"./Fib.rds")
Fib <- readRDS("./Fib.rds")
Fib$Diets
Fib_69<- FindMarkers(Fib, ident.1 = "M6-PRD", ident.2 = "M9-PRD",
                     group.by = "Diets")
head(Fib_69)
?FindMarkers
Fib_bp<- FindMarkers(Fib, ident.1 = "M9-BRD", ident.2 = "M9-PRD",
                     group.by = "Diets")
View(Fib_bp)
mean_Fib <- apply(as.data.frame(Fib@assays$RNA@data),1,function(x){tapply(x,Fib$Diets,median)}) %>% t() %>% as.data.frame()
Fib_rescue <- read.delim("~/scRNA_aging/TFs/output_Fib_rescue/Step2_regulonTargetsInfo.tsv")
head(Fib_rescue)
Fib_rescue <- Fib_rescue[Fib_rescue$highConfAnnot=="TRUE",]
Fib_resistance <- read.delim("~/scRNA_aging/TFs/output_Fib_resistance//Step2_regulonTargetsInfo.tsv")
Fib_resistance <- Fib_resistance[Fib_resistance$highConfAnnot=="TRUE",]

Fib_bp_markers$SYMBOL <- rownames(Fib_bp_markers)


Fib_bp_markers$type <- ifelse(Fib_bp_markers$p_val_adj < 0.05,
                   ifelse(abs(Fib_bp_markers$avg_log2FC) > 0.25 ,
                          ifelse(Fib_bp_markers$avg_log2FC < -0.25 ,'down','up'),'noSig'),'noSig')


table(Fib_bp_markers$type)
View(Fib_bp_markers)
meta <- Fib_bp_markers[Fib_bp_markers$SYMBOL%in%c("Col3a1","Col5a2","Timp1","Jund"),]

p=ggplot(na.omit(as.data.frame(Fib_bp_markers)),aes(x = avg_log2FC,y = -log10(p_val_adj))) +
  geom_point(aes(color = type),size = 3,alpha=0.3) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#F39B7F','noSig'='grey','down'='#4BB5CE'),
                     # legend labels
                     label = c('up'='up (num=422)','down'='down (num=2864)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('9m: PRD vs BRD')+theme_classic()+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(plot.title = element_text(size = 24))
p
p + geom_text_repel(data = meta,aes(x = avg_log2FC,y = -log10(p_val_adj),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)




Fib_69_markers$SYMBOL <- rownames(Fib_69_markers)


Fib_69_markers$type <- ifelse(Fib_69_markers$p_val_adj < 0.05,
                              ifelse(abs(Fib_69_markers$avg_log2FC) > 0.25 ,
                                     ifelse(Fib_69_markers$avg_log2FC < -0.25 ,'down','up'),'noSig'),'noSig')


table(Fib_69_markers$type)
View(Fib_69_markers)
meta <- Fib_69_markers[Fib_69_markers$SYMBOL%in%c("Col3a1","Col1a1","Col1a2","Col5a1","Timp1"),]

p=ggplot(na.omit(as.data.frame(Fib_69_markers)),aes(x = avg_log2FC,y = -log10(p_val_adj))) +
  geom_point(aes(color = type),size = 3,alpha=0.3) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#F39B7F','noSig'='grey','down'='#4BB5CE'),
                     # legend labels
                     label = c('up'='up (num=587)','down'='down (num=2889)')) +
  # 阈值分割线
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  ggtitle('PRD: 9m vs 6m')+theme_classic()+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(),	
        legend.text = element_text(size=14),legend.position = "top",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(plot.title = element_text(size = 24))
p
p + geom_text_repel(data = meta,aes(x = avg_log2FC,y = -log10(p_val_adj),label = meta$SYMBOL),
                    force=20,color="black",size=6,point.padding = 0.5,hjust = 0.5,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
                    segment.color="black",segment.size=0.2,nudge_y=1)





#mean_Fib_s <- mean_Fib[c(a,b),c(3,4,6)]

n <- t(scale(t(as.data.frame(Fib@assays$RNA@data))))
n <- na.omit(n)
mean_Fib <- apply(as.data.frame(n),1,function(x){tapply(x,Fib$Diets,mean)}) %>% t() %>% as.data.frame()
head(mean_Fib)
#p <- pheatmap::pheatmap(na.omit(n),show_rownames = F,cluster_cols = F,fontsize = 14,
#                        cluster_rows = T,angle_col = "0")
#p

mean_Fib_s <- mean_Fib[c(unique(c(Fib_resistance$TF,Fib_resistance$gene)),unique(c(Fib_rescue$TF,Fib_rescue$gene))),c(3,4,6)]

library(ComplexHeatmap)
rownames(mean_Fib_s)
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))

h1 <- Heatmap(t(mean_Fib_s[intersect(rownames(mean_Fib_s),unique(c(Fib_resistance$TF,Fib_resistance$gene))),]),name="scale",row_dend_reorder = T,row_title = "Fibroblast-like",
              col=col_fun,column_title = "Resistance",show_column_names  = T)
h2 <- Heatmap(t(mean_Fib_s[intersect(rownames(mean_Fib_s),unique(c(Fib_rescue$TF,Fib_rescue$gene))[-c(4,5,15,18,26,31)]),]),name="scale",row_dend_reorder = T,row_title = "Fibroblast-like",
              col=col_fun,column_title = "Rescue",show_column_names  = T)
h1+h2
a
b

DotPlot(object = Fib,group.by = "Diets",
        cols = c("#4B549B","#E91C22"), 
        features = intersect(rownames(mean_Fib_s),unique(c(Fib_resistance$TF,Fib_resistance$gene))))+
  theme_classic()+
  #  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
  #                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.position = "right",
        legend.title = element_text(size = 13))+xlab("")+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))+ylab("")

DotPlot(object = Fib,group.by = "Diets",
        cols = c("#4B549B","#E91C22"), 
        features = intersect(rownames(mean_Fib_s),unique(c(Fib_rescue$TF,Fib_rescue$gene))))+
  theme_classic()+
  #  scale_y_discrete(limits=c("GC","SC","EC", "EDC","GLC",
  #                            "M","DC","B","NK","T","ILC"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25,hjust = 1),
        axis.text = element_text(size=14,colour = "black"),
        legend.text =element_text(size=12),legend.position = "right",
        legend.title = element_text(size = 13))+xlab("")+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))



left_anno = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#0ECBC5","#FFDD93")),
                       labels =c("resistance","rescue"),
                       labels_gp = gpar(col = "black", fontsize = 14)))

col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))
n=t(scale(t(mean_Fib_s)))
?Heatmap
p1 <- Heatmap(na.omit(n),col=col_fun,border = "black",
              show_row_names = FALSE,show_column_names = F,
              cluster_columns = FALSE,cluster_rows = T,
              row_split = c("resistance","rescue"),
              left_annotation = left_anno,
              row_title = NULL,column_title = NULL,
              heatmap_legend_param = list(
              title = "scale",
              title_position = "leftcenter-rot"))

p1




View(Fib_69)
953
3834
Fib_69 <- Fib_69_markers 
Fib_bp <-  Fib_bp_markers

#Fib_69_markers <- Fib_69
#Fib_bp_markers <- Fib_bp
a=intersect(rownames(Fib_69[which(Fib_69$avg_log2FC>0.25),]),
            rownames(Fib_bp[which(Fib_bp$avg_log2FC>0.25),]))
a
b=intersect(rownames(Fib_69[which(Fib_69$avg_log2FC< -0.25),]),
            rownames(Fib_bp[which(Fib_bp$avg_log2FC< -0.25),]))
b
length(a)
View(as.data.frame(b))

intersect(a,b)
sym<- AnnotationDbi::select(org.Mm.eg.db,keys=a,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_rescue <- enrichKEGG(gene =sym$ENTREZID,
                       organism = 'mmu',
                       pvalueCutoff =0.05)

write.table(KEGG_rescue@result, 'KEGG_rescue.txt', sep='\t', quote=F)

sym<- AnnotationDbi::select(org.Mm.eg.db,keys=b,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_restrict <- enrichKEGG(gene =sym$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)

write.table(KEGG_restrict@result, 'KEGG_restrict.txt', sep='\t', quote=F)
KEGG_restrict <- read.delim("./KEGG_restrict.txt", row.names=1)

s <- KEGG_restrict[KEGG_restrict$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
s$Description
dim(s)
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD-resistance") +
  geom_text(aes(x = rep(0.002,6), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)



KEGG_rescue <- read.delim("./KEGG_rescue.txt", row.names=1)

s <- KEGG_rescue[KEGG_rescue$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
s$Description
dim(s)
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd",direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD-rescue") +
  geom_text(aes(x = rep(0.002,7), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)





rm(BRD_STC,BRD_GC,BRD_DC,BRD_EDC,BRD_T,BRD_SCTC,PRD_STC,PRD_GC,PRD_DC,PRD_EDC,PRD_T,PRD_SCTC)
gc()
BP  <- readRDS("./BP_subset.rds")
BP_GC <- subset(BP,cells=rownames(BP@meta.data[BP$abbreviation=="GC",]))
rm(BP)
gc()


data <- as.data.frame(BP_GC@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M6-BRD","M9-PRD","M9-BRD"))

head(ggdata)
aging_gene_sets <- read.delim("~/scRNA_aging/output/aging_gene_sets.txt")

data <- ggdata[ggdata$SYMBOL%in%aging,]

ggplot(na.omit(data), aes(x = group, y =log(value+1),fill=group)) +
  geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("M6-PRD","M6-BRD"),c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
    map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 5,
    tip_length = 0,vjust = 0.1)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))

mean_GC <- apply(as.data.frame(BP_GC@assays$RNA@data),1,function(x){tapply(x,BP_GC$group_age_diet,mean)}) %>% t() %>% as.data.frame()
head(mean_GC)
dim(mean_GC)
colnames(mean_GC)=strsplit2(colnames(mean_GC),"_")[,1]

data <- mean_GC[rownames(mean_GC)%in%c("Cdkn2a","Cdkn2d","Cdkn1b","Cdkn1a"),]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1
rownames(data)
pheatmap::pheatmap(n,fontsize = 14,border_color = "black",cluster_cols = F)




####################### cellphone ############################

###################### BRD ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_BRD/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))

head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.01))
#View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
statdf$total_number=0
#View(statdf)
for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}

rankname=sort(unique(statdf$indexa)) 
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)
#View(statdf) # 226
unique(statdf$indexb)
p1 <- statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,230))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)
p1
###################### PRD ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_PRD/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))

head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.001))
#View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")


statdf$total_number=0
#View(statdf)
for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}

rankname=sort(unique(statdf$indexa)) 
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)
#View(statdf) # 226
unique(statdf$indexb)
p2 <- statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,230))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)
p2

###################### 6M ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M6/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))

head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.001))
#View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")


statdf$total_number=0
#View(statdf)
for (i in 1:dim(statdf)[1]) {
  tmp_indexb=statdf[i,"indexb"]
  tmp_indexa=statdf[i,"indexa"]
  if (tmp_indexa == tmp_indexb) {
    statdf[i,"total_number"] = statdf[i,"number"]
  } else {
    statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
      statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
  }
}

rankname=sort(unique(statdf$indexa)) 
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)
#View(statdf) # 226
unique(statdf$indexb)
p3 <- statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,230))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)
p3


M6 <- readRDS("./M6_BP.rds")
meta <- M6@meta.data[M6@meta.data$abbreviation2%in%c("STC","GC"),]
STGC <- subset(M6,cells=rownames(meta))
head(STGC@meta.data)

STGC_marker <- FindAllMarkers(object =STGC,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
STGC_top10  <- STGC_marker  %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)
STGC$abbreviation2
View(STGC_marker)
DoHeatmap(STGC,features = STGC_top20$gene)
DoHeatmap(STGC,
          features = as.character(unique(STGC_top10$gene)),
          group.by = "abbreviation2",
          group.colors = c("#370335","#FD7446"))+
  scale_fill_gradientn(colors = c("white","white","firebrick3"))


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=STGC_marker[STGC_marker$cluster=="Granulosa cells",]$gene,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)
View(KEGG_GC@result)
BP_GC<- clusterProfiler::simplify(enrichGO(gene=ncbi$SYMBOL,keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))
View(BP_GC@result)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=STGC_marker[STGC_marker$cluster=="Stem cells",]$gene,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")

BP_STC<- clusterProfiler::simplify(enrichGO(gene=ncbi$SYMBOL,keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_STC@result)
BP1 <- BP_GC@result[BP_GC@result$Description%in%c("sex differentiation","Wnt signaling pathway",
                                                  "connective tissue development","hormone-mediated signaling pathway",
                                                  "ovulation cycle","response to estradiol",
                                                  "cellular response to fatty acid"),]
BP1$class <- "GC"

BP2 <- BP_STC@result[BP_STC@result$Description%in%c("RNA splicing","oxidative phosphorylation",
                                                    "regulation of mitotic cell cycle",
                                                    "DNA replication","female gamete generation"),]
BP2$class <- "STC"
BP <- rbind(BP1,BP2)
b <- lapply(str_split(BP$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
BP$GeneRatio <- b
ggplot(BP,aes(class,Description,size = GeneRatio))+
  geom_point(shape=21,aes(fill= pvalue),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,8))+
  theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black")
  ) +scale_x_discrete(limits=c("GC","STC"))+
  scale_y_discrete(limits= unique(BP$Description)) 



df <- readRDS("./df_monocle.rds")
head(df)
M9<- readRDS("./M9_BP.rds")
data <- M9@assays$RNA@data %>% as.data.frame(.)%>%
  .[rownames(.)%in%df$wide.res$gene,]


M9@meta.data[M9@meta.data$abbreviation2=="GC",]
ggdata <- melt(data)
head(ggdata)
ggdata <- ggdata[intersect(M9@meta.data[M9@meta.data$abbreviation2=="GC",]%>%rownames(),ggdata$variable),]
head(ggdata)
ggdata$diets <- paste0(strsplit2(ggdata$variable,"-")[,2],"RD")

ggplot(ggdata, aes(x = diets, y =log(value+1),fill=diets)) +
  geom_violin(aes(fill = diets)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD","PRD")),
              map_signif_level = T,step_increase = 0,size = 0.5,textsize = 5,
              tip_length = 0,vjust = 0)+ylab("")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.8),axis.line.y=element_line(size=0.8))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))+ggtitle("GC")


unique(M9$abbreviation)
M9 <- readRDS("./M9_BP.rds")
M9$abbreviation2 <- gsub("SC","SCTC",M9$abbreviation)%>%gsub("TCC","SCTC",.)
unique(M9$abbreviation2)
unique(M9$abbreviation2)
M9$abbreviation2 <- factor(M9$abbreviation2,levels = c("STC",
                                                       "GC","SCTC","EC","EDC",
                                                       "GLC","M","DC","B","NK","T","ILC"))
unique(M9$abbreviation2)
M9@meta.data$age <- substr(M9$group_age,1,2)
M9$Group <- paste0(M9$age,"-",M9$abbreviation2)
M9@active.ident <- M9$abbreviation2
STGC <- subset(M9,idents=c("STC","GC"))
data <- as.data.frame(STGC@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M9-BRD","M9-PRD"))
head(ggdata)
df$wide.res$SYMBOL <- df$wide.res$gene
data <- merge(df$wide.res,ggdata,all.x=T,by="SYMBOL")%>%na.omit(.)
data1 <- data[data$value>0.1,]
data1$group <- substr(data1$group,4,6)
data1$group <- factor(data1$group,levels = c("BRD","PRD"))
data1$cluster <- paste0("Cluster",data1$cluster)
data1$variable
STGC@meta.data$variable <- rownames(STGC@meta.data)
data1 <- merge(data1,STGC@meta.data,all.x=T,by="variable")
data2 <-data1[data1$abbreviation=="STC",] 
head(data2)

cluster <- data2[data2$cluster=="Cluster1",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data2[data2$cluster=="Cluster2",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data2[data2$cluster=="Cluster3",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data2[data2$cluster=="Cluster4",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
data

ggplot(data2,aes(x = group,y = value,fill=group))+geom_violin(aes(fill = group)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD","PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0.1,test="t.test")+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 16),
        strip.text.x = element_text(size=16),strip.background = element_blank())+
  facet_wrap(~cluster,ncol=2,scales = "free_y")


######################### Immune ##############################
BP <- readRDS("./BP_subset.rds")
unique(BP@active.ident)
IM <- subset(BP,idents=c("B cells","T cells","Natural killer cells",
                         "Innate lymphoid cells","Macrophages",
                         "Dendritic cells","Granulocytes"))
rm(BP)
gc()
ElbowPlot(IM, ndims=50, reduction="pca") 
IM <- FindNeighbors(IM, dims = 1:30)
IM <- FindClusters(IM, resolution = 0.6)
IM <- RunUMAP(IM, dims = 1:30)
#IM <- RunTSNE(IM, dims = 1:30)
DimPlot(IM, reduction = "umap", pt.size = 1, label = T)
# CD300e+ Mono 16
# ILC2 15
# T17LC 4 9
# NTL/Mono 5
# CD4 T 13
# DC 10 11 17
# NK 3
# B 7
# T1LC 1
# TMacro 0 2 6
# CD45 8
DimPlot(IM, reduction = "umap", pt.size = 1, label = T,group.by = "seurat_clusters")

DotPlot(IM,features = c("Ptprc","C1qb","C1qa","Apoe","Mafb","Fcgr1","Axl",
                        "Cd28","Tbx21","Ifng","Igkc","Cd79a","Cd79b",
                        "Ncr1","Klrb1b","Klrb1c","Cd209a","Flt3","Clec9a",
                        "Lef1","Tcf7","Ccr7","F10","Clec4a","Hp",
                        "Rorc","Il18r1","Zbtb16","Cenpe","Top2a","Ccna2",
                        "H2-Eb1","H2-Aa","Cd74","Il1rl1","Areg","Gata3",
                        "Pparg","Cd300e","Fcgr4"),group.by = "seurat_clusters")+theme(axis.text.x = element_text(angle = 45,hjust =1))
IM <- subset(IM,idents=c(0,1,2,3,4,5,6,7,8,9,10,11,13,15,16,17))
DotPlot(IM,features = c("Ptprc","C1qb","C1qa","Apoe","Mafb","Fcgr1","Axl",
                        "Cd28","Tbx21","Ifng","Igkc","Cd79a","Cd79b",
                        "Ncr1","Klrb1b","Klrb1c","Cd209a","Flt3","Clec9a",
                        "Lef1","Tcf7","Ccr7","F10","Clec4a","Hp",
                        "Rorc","Il18r1","Zbtb16","Cenpe","Top2a","Ccna2",
                        "H2-Eb1","H2-Aa","Cd74","Il1rl1","Areg","Gata3",
                        "Pparg","Cd300e","Fcgr4"),cols = c("white","#584EA1"))+theme(axis.text.x = element_text(angle = 45,hjust =1))
# CD300e+ Mono 16
# ILC2 15
# T17LC 4 9
# NTL/Mono 5
# CD4 T 13
# DC 10 11 17
# NK 3
# B 7
# T1LC 1
# TMacro 0 2 6
# CD45 8
new.cluster.ids <- c("M","T1LC","M","NK","T17LC","Mo","M","B","IML","T17LC",
                     "DC","DC","CD4T","ILC2","CD300e+Mo","DC")
names(x=new.cluster.ids)=levels(x = IM)

IM1 <- RenameIdents(object =IM, new.cluster.ids)

DotPlot(IM1,features = c("Ptprc","C1qb","C1qa","Apoe","Mafb","Fcgr1","Axl",
                         "Cd28","Tbx21","Ifng","Igkc","Cd79a","Cd79b",
                         "Ncr1","Klrb1b","Klrb1c","Cd209a","Flt3","Clec9a",
                         "Lef1","Tcf7","Ccr7","F10","Clec4e","Hp",
                         "Rorc","Il18r1","Zbtb16",
                         "H2-Eb1","H2-Aa","Cd74","Il1rl1","Areg","Gata3",
                         "Pparg","Cd300e","Fcgr4"),cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 14))+
  scale_y_discrete(limits=c("M","T1LC","B","NK","DC" ,"CD4T" ,"Mo","T17LC",
                            "IML","ILC2","CD300e+Mo"))


IM1 <- RunUMAP(IM1, dims = 1:30)
saveRDS(IM1,"./IM1_BP.rds")

IM1@active.ident <- factor(IM1@active.ident,
                           levels = c("M","T1LC","B","NK","DC" ,"CD4T" ,"Mo","T17LC",
                                      "IML","ILC2","CD300e+Mo"))
DimPlot(IM1, reduction = "umap", pt.size = 1, label = F)+
  scale_color_futurama()
saveRDS(IM1,"./IM1_BP.rds")


df1=as.data.frame(IM1@active.ident)
IM1@meta.data$celltype=df1$`IM1@active.ident`
prop_test <- sc_utils(IM1)
prop_test2 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6-PRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)

permutation_plot2(prop_test2,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+ggtitle("PRD: 9m vs 6m")


prop_test3 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test3,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+ggtitle("9m: PRD vs BRD")

?FindMarkers
IM1$D_name <- paste0(IM1$Diets,"-",IM1$celltype)
idents <- gsub("CD300e+Mo","Mo",IM1$D_name)
names(idents) <- rownames(IM1@meta.data)
IM1@active.ident<-factor(idents)

PRD_M_markers <- FindMarkers(IM1, ident.1 = "M9-PRD-M", ident.2 = "M6-PRD-M", min.pct = 0.25)
PRD_M_markers$group <- 0
PRD_M_markers[PRD_M_markers$avg_log2FC > 0,]$group <- "9M"
PRD_M_markers[PRD_M_markers$avg_log2FC< 0,]$group <- "6M"


PB_M_markers <- FindMarkers(IM1, ident.1 = "M9-PRD-M", ident.2 = "M9-BRD-M", min.pct = 0.25)
PB_M_markers$group <- 0
PB_M_markers[PB_M_markers$avg_log2FC > 0,]$group <- "PRD"
PB_M_markers[PB_M_markers$avg_log2FC< 0,]$group <- "BRD"

unique(IM1$D_name)

PRD_Mo_markers <- FindMarkers(IM1, ident.1 = "M9-PRD-Mo", ident.2 = "M6-PRD-Mo", min.pct = 0.25)
PRD_Mo_markers$group <- 0
PRD_Mo_markers[PRD_Mo_markers$avg_log2FC > 0,]$group <- "9M"
PRD_Mo_markers[PRD_Mo_markers$avg_log2FC< 0,]$group <- "6M"

PB_Mo_markers <- FindMarkers(IM1, ident.1 = "M9-PRD-Mo", ident.2 = "M9-BRD-Mo", min.pct = 0.25)
PB_Mo_markers$group <- 0
PB_Mo_markers[PB_Mo_markers$avg_log2FC > 0,]$group <- "PRD"
PB_Mo_markers[PB_Mo_markers$avg_log2FC< 0,]$group <- "BRD"


M_age_up <- intersect(PRD_M_markers[PRD_M_markers$group=="9M",]%>%rownames(.),
                       PB_M_markers[PB_M_markers$group=="PRD",]%>%rownames(.))

M_age_down <-intersect(PRD_M_markers[PRD_M_markers$group=="6M",]%>%rownames(.),
                        PB_M_markers[PB_M_markers$group=="BRD",]%>%rownames(.))

Mo_age_up <- intersect(PRD_Mo_markers[PRD_Mo_markers$group=="9M",]%>%rownames(.),
                      PB_Mo_markers[PB_Mo_markers$group=="PRD",]%>%rownames(.))

Mo_age_down <-intersect(PRD_Mo_markers[PRD_Mo_markers$group=="6M",]%>%rownames(.),
                       PB_Mo_markers[PB_Mo_markers$group=="BRD",]%>%rownames(.))


ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=M_age_up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_M_resistance <- enrichKEGG(gene =ncbi$ENTREZID,
                      organism = 'mmu',
                      pvalueCutoff =0.05)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=M_age_down,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_M_rescue <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=Mo_age_up,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Mo_resistance <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,
                             keys=Mo_age_down,
                             columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Mo_rescue <- enrichKEGG(gene =ncbi$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)
write.csv(KEGG_M_resistance@result,"./KEGG_M_resistance.csv")
write.csv(KEGG_M_rescue@result,"./KEGG_M_rescue.csv")
write.csv(KEGG_Mo_resistance@result,"./KEGG_Mo_resistance.csv")
write.csv(KEGG_Mo_rescue@result,"./KEGG_Mo_rescue.csv")


BP_M_resistance<- clusterProfiler::simplify(enrichGO(gene=M_age_up,keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

BP_M_rescue<- clusterProfiler::simplify(enrichGO(gene=M_age_down,keyType = "SYMBOL",
                                                     OrgDb= org.Mm.eg.db, ont="BP",
                                                     pAdjustMethod="BH",pvalueCutoff=0.05,
                                                     qvalueCutoff=0.05))

BP_Mo_resistance<- clusterProfiler::simplify(enrichGO(gene=Mo_age_up,keyType = "SYMBOL",
                                                     OrgDb= org.Mm.eg.db, ont="BP",
                                                     pAdjustMethod="BH",pvalueCutoff=0.05,
                                                     qvalueCutoff=0.05))

BP_Mo_rescue<- clusterProfiler::simplify(enrichGO(gene=Mo_age_down,keyType = "SYMBOL",
                                                 OrgDb= org.Mm.eg.db, ont="BP",
                                                 pAdjustMethod="BH",pvalueCutoff=0.05,
                                                 qvalueCutoff=0.05))

write.csv(BP_M_resistance@result,"./BP_M_resistance.csv")
write.csv(BP_M_rescue@result,"./BP_M_rescue.csv")
write.csv(BP_Mo_resistance@result,"./BP_Mo_resistance.csv")
write.csv(BP_Mo_rescue@result,"./BP_Mo_rescue.csv")
unique(IM1$D_name)
#"Ccr2","Csf2ra","Csf1r","Ighm","Clec4a2","Clec4a3" 



BP_M_resistance <- read.csv("./BP_M_resistance.csv")
BP_M_rescue <- read.csv("./BP_M_rescue.csv")
s <- BP_M_resistance[BP_M_resistance$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
s$Description
dim(s)

ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD-resistance") +
  geom_text(aes(x = rep(0.002,12), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)

s <- BP_M_rescue[BP_M_rescue$select==1,]%>%na.omit()
s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b
s <- s[order(s$GeneRatio),]
s$Description
dim(s)
ggplot(data = s,aes(x = GeneRatio, y = Description)) +
  geom_bar(aes(fill = -log10(pvalue)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd",direction = 1) +
  labs(x = "RichFactor", y = "Pathway", title = "KEGG enrichment for BRD-rescue") +
  geom_text(aes(x = rep(0.002,8), #用新增的重复数组控制文本标签起始位置
                label= Description),hjust= 0)+theme_classic() +
  theme(axis.text.y = element_blank(),axis.ticks.y  = element_blank(),
        axis.text= element_text(size=12))+
  scale_y_discrete(limits=s$Description)
ElbowPlot(IM2, ndims=50, reduction="pca")

IM2 <- FindNeighbors(IM2, dims = 1:40)
IM2 <- FindClusters(IM2, resolution = 0.3)
IM2 <- RunUMAP(IM2, dims = 1:40)

DimPlot(IM2, reduction = "umap", pt.size = 1, 
        label = T,group.by = "seurat_clusters")+scale_color_jama()

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 
library(DoubletFinder)
?paramSweep
IM2 <- NormalizeData(object =IM2, normalization.method = "LogNormalize", scale.factor = 1e4)
sweep.data <- paramSweep(IM2,PCs=1:40,sct=T)
sweep.stats <- summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)
homotypic.prop=modelHomotypic(IM2@meta.data$RNA_snn_res.0.3)
nExp_poi=round(0.075*length(IM2$orig.ident))
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
IM3=doubletFinder(IM2, PCs = 1:40, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),
                      nExp = nExp_poi.adj, reuse.pANN = FALSE,sct = T)
dim(IM3@meta.data)
IM3@meta.data$DF_hi.lo<- IM3@meta.data[,21]   #注意，这里是按照之前教程的顺序没有添加其他列，所以Doublet排在第8列，如果在其他列，请更改相关代码

DimPlot(IM3, group.by ="DF_hi.lo",cols=c("black","gold","red"),reduction = "umap")

#绘图，结果用金色代表Singlet，黑色代表Doublet。添加红色是为了防止出现为添加相关标签的NA值准备，通常结果不会存在红色的点，如果存在说明结果有误。
Doublet<-table(IM3@meta.data$DF_hi.lo=="Doublet")
Idents(IM3) <- "DF_hi.lo"    #激活相关ID
IM3<-subset(x = IM3, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet
DimPlot(IM3, reduction = "tsne", pt.size = 1, 
        label = T,group.by = "seurat_clusters")+scale_color_jama()

IM3 <- FindNeighbors(IM3, dims = 1:40)
IM3 <- FindClusters(IM3, resolution = 0.6)
IM3 <- RunUMAP(IM3, dims = 1:40)

FeaturePlot(IM3,reduction = "umap",features = c("Retnla","Ly6c1","Tnf","Il1b","Cd86","Cd80","Ccl2",
                             "Arg1","Il10","Mrc1","Cd163","Mgl2"))

DotPlot(IM3,features = c("Ly6c1","Tnf","Il1b","Il6","Cd86","Cd80","Ccl2","Il12b",
                             "Arg1","Il10","Mrc1","Cd163","Cd209a"),group.by = "seurat_clusters")

FeaturePlot(IM3,reduction = "tsne",features = c("Il18",
                                                "Il6","Il1b","Il10","Tnf"))



DotPlot(IM3,features = c("Ccl2","Ccl3","Ccl5","Ccl7","Ccl12"),
        group.by = "Diets")

#### 吞噬作用
DotPlot(IM3,features = c("Cd68","Cd36","Fcgr1","Mrc1","Cd14"),
        group.by = "Diets",cols = c("#584EA1","#ED1C26"))+
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))+coord_flip()+ggtitle("Phagocytosis")
### Cell proliferation 
DotPlot(IM3,features = c("Pcna","Mki67"),
        group.by = "Diets",cols = c("#584EA1","#ED1C26"))+
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))+coord_flip()+ggtitle("Phagocytosis")

IM1 <- readRDS("./IM1_BP.rds")
IM1$celltype <- IM1@active.ident
cells <- IM1@meta.data[IM1$celltype%in%c("Mo","M"),]%>%rownames(.)
IM2 <- subset(IM1,cells=cells)
IM2$celltype <- factor(IM2$celltype,levels = c("Mo","M"))
IM2 <- RunUMAP(IM2, dims = 1:40)
DimPlot(IM2, reduction = "umap", pt.size = 1, 
        label = T)+scale_color_jama()

FeaturePlot(IM2,reduction = "umap",features = c("Retnla","Ly6c1","Tnf","Il1b","Cd86","Cd80","Ccl2",
                                                "Arg1","Il10","Mrc1","Cd163","Mgl2"))
FeaturePlot(IM2,reduction = "umap",features = c("Cd86","Arg1","Il10"))

library(monocle)
data <- as(as.matrix(IM2@assays$RNA@data), 'sparseMatrix')
pd <- IM2@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
pd <- new("AnnotatedDataFrame",
          data=pd)
fd <- new("AnnotatedDataFrame",
          data=fData)
cd <- newCellDataSet(data,
                     phenoData = pd,
                     featureData =fd)#, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower = 0.1))
rm(data)
gc()
sc_cds <- estimateSizeFactors(cd)
sc_cds <- estimateDispersions(sc_cds)
sc_cds$celltype
diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = T)
plot_cell_trajectory(sc_cds, color_by ="celltype",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14),legend.title = element_text(size=16))+
  theme(legend.position = "right")
plot_cell_trajectory(sc_cds,color_by = "Pseudotime")+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))+theme(legend.position = "right")

plot_cell_trajectory(sc_cds,color_by = "State")+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))+theme(legend.position = "right")

saveRDS(sc_cds,"./cell_M_cds.rds")


s.genes <- c("Retnla","Ly6c1","Tnf","Il1b","Cd86","Cd80","Ccl2",
             "Arg1","Il10","Mrc1","Cd163","Mgl2")
p1 <- plot_genes_jitter(sc_cds[s.genes,],grouping = "State",color_by = "State")
p1


disp_table <- dispersionTable(sc_cds)
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
sc_cds <- setOrderingFilter(sc_cds, unsup_clustering_genes$gene_id)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)

plot_cell_trajectory(sc_cds, color_by ="celltype",cell_size = 0.5,splite_by="Diets")+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

plot_cell_trajectory(sc_cds, color_by ="State",cell_size = 0.5,splite_by="Diets")+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))






#devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
BEAM_res=BEAM(sc_cds[ordering_genes,],branch_point = 1,cores =2)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
calculate_NB_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
{
  expr_hint <- expr_selection_func(f_expression)
  if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
    disp_guess_fit <- disp_func(expr_hint)
    
    # For NB: Var(Y)=mu*(1+mu/k)
    f_expression_var <- var(f_expression)
    f_expression_mean <- mean(f_expression)
    
    disp_guess_meth_moments <- f_expression_var - f_expression_mean
    disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
    
    #return (max(disp_guess_fit, disp_guess_meth_moments))
    return (disp_guess_fit)
  }
  return (NULL)
}
library("RColorBrewer", lib.loc="C:/Users/tutu/AppData/Local/R/win-library/4.2")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

df1 <- plot_genes_branched_heatmap(sc_cds[row.names(subset(BEAM_res,qval<1e-4)),],
                                   branch_point = 1,
                                   num_clusters = 4, #这些基因被分成几个group
                                   cores = 1,
                                   branch_labels = c("Cell fate 1", "Cell fate 2"),
                                   #hmcols = NULL, #默认值
                                   hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                   branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                   use_gene_short_name = T,
                                   show_rownames = F,
                                   return_heatmap = T #是否返回一些重要信息
)
df1
saveRDS(df1,"./monocle_df1.rds")
visCluster(object = df1,plot.type = "both")
df1=plot_pseudotime_heatmap2(sc_cds[ordering_genes[1:200],],
                             num_clusters = 4,
                             cores = 1)
visCluster(object = df1,plot.type = "line")
df2 <- plot_multiple_branches_heatmap2(sc_cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                       branches = c(1,3,4,5),
                                       num_clusters = 4,
                                       cores = 1,
                                       use_gene_short_name = T,
                                       show_rownames = T)
visCluster(object = df,plot.type = "heatmap")
# enrich
df <- plot_multiple_branches_heatmap2(sc_cds[row.names(subset(BEAM_res,qval < 0.1)),],
                                      branches = c(1,3,4,5),
                                      num_clusters = 4,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T)

# enrich for clusters
library(org.Mm.eg.db)
enrich <- enrichCluster(object = df1,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.1,
                        topn = 5,
                        seed = 5201314)

write.table(enrich,"enrich.txt",sep="\t",quote=F)
# check
head(enrich[1:3,])

enrich <- read.delim("F:/A_aging_lnc/enrich.txt", row.names=1)
markGenes = sample(unique(df1$wide.res$gene),25,replace = F)

# PLOT
df1
visCluster(object = df1,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           go.col = rep(jjAnno::useMyCol("calm",n = 4),each = 5),
           add.bar = T,
           line.side = "left")
