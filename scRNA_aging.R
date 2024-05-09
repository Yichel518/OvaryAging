library(Seurat)
library(SingleR)
setwd("/home/tuyx/scRNA_aging/output/")
folders=c("M6-B-1","M6-B-2","M6-B-3","M6-W-1","M6-W-2","M6-W-3","M6-P-1","M6-P-2","M6-P-3",
          "M9-B-1","M9-B-2","M9-B-3","M9-W-1","M9-W-2","M9-W-3","M9-P-1","M9-P-3")
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,min.cells = 3, 
                     min.features =200)
})
sce.big <- merge(sceList[[1]], 
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],sceList[[6]],
                       sceList[[7]],sceList[[8]],sceList[[9]],
                       sceList[[10]],sceList[[11]],sceList[[12]],sceList[[13]],sceList[[14]],
                       sceList[[15]],sceList[[16]],sceList[[17]]), 
                 add.cell.ids = folders, 
                 project = "ovary")
rm(sceList)
gc()
sce.big[["percent.mito"]] <- PercentageFeatureSet(object = sce.big, pattern = "^mt-") #线粒体基因以mt开头的选择出来
sce.big.mt=sce.big[["percent.mito"]]
VlnPlot(object = sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
plot1 <- FeatureScatter(object = sce.big, feature1 = "nCount_RNA", feature2 = "percent.mito") 
plot2 <- FeatureScatter(object = sce.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
#CombinePlots(plots = list(plot1,plot2))
sce.big <- subset(x = sce.big, subset = nFeature_RNA > 200 & nFeature_RNA<6000 & percent.mito<15)
sce.big <- SCTransform(sce.big,  verbose = FALSE,conserve.memory = T)
gc()
sce.big <- NormalizeData(object =sce.big, normalization.method = "LogNormalize", scale.factor = 1e4)
sce.big <- FindVariableFeatures(sce.big, selection.method = "vst")

sce.big <- RunPCA(sce.big, features = VariableFeatures(sce.big))
DimPlot(sce.big, reduction = "pca", group.by="orig.ident")
ElbowPlot(sce.big, ndims=50, reduction="pca") 

sce.big <- FindNeighbors(object = sce.big, dims = 1:40)
sce.big <- FindClusters(object = sce.big, resolution = 0.6)
sce.big <- RunUMAP(object = sce.big, dims = 1:40)
DimPlot(object = sce.big, reduction = "umap",label = T,raster=FALSE)
sce.big <- RunTSNE(object = sce.big, dims = 1:40)
DimPlot(object = sce.big, reduction = "tsne",label = T,label.size = 6,raster=FALSE)
DimPlot(object = sce.big, reduction = "tsne",label = F,group.by="orig.ident") 
meta=sce.big@meta.data # pbmc的meta文件，包含了seurat 的聚类结果
ImmGen.se=ImmGenData() #(鼠)
Mouse.se=MouseRNAseqData() #(鼠)


sce.big_for_SingleR <- GetAssayData(sce.big, slot="data") ##获取标准化矩阵
sce.big.hesc <- SingleR(test = sce.big_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
sce.big.hesc2 <- SingleR(test = sce.big_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
sce.big@meta.data$SingleR_Mouse <- sce.big.hesc$labels
sce.big@meta.data$SingleR_ImmGen <- sce.big.hesc2$labels
DimPlot(object = sce.big, reduction = "tsne",label =T,group.by = "SingleR_Mouse",raster=FALSE)
DimPlot(object = sce.big, reduction = "tsne",label =T,group.by = "SingleR_ImmGen",raster=FALSE)
DimPlot(object = sce.big, reduction = "umap",label =T,group.by = "SingleR_Mouse",raster=FALSE)
DimPlot(object = sce.big, reduction = "umap",label =T,group.by = "SingleR_ImmGen",raster=FALSE,label.size = 5)
DimPlot(object = sce.big, reduction = "tsne",label = T,raster=FALSE)
DimPlot(object = sce.big, reduction = "umap",label = T,raster=FALSE,label.size = 5)



DotPlot(object = sce.big, features = c("Ddx4","Dazl"))
FeaturePlot(object = sce.big, features = c("Cd133"),order=T,reduction = "tsne")
FeaturePlot(object = sce.big, features = c("Dazl"),order=T,reduction = "tsne")

DotPlot(object = sce.big, features = c("Nanog","Par6"))

# Theca 1 8 11 15 26
DotPlot(object = sce.big, features = c("Cyp17a1","Cyp11a1","Star","Hsd3b1","Lhcgr","Wt1"))
DotPlot(object = all, features = c("Cyp17a1","Cyp11a1","Star","Hsd3b1","Lhcgr","Wt1"))
FeaturePlot(object = sce.big, features = c("Cyp17a1","Cyp11a1","Star"),order=T,reduction = "umap")
FeaturePlot(object = all, features = c("Cyp17a1","Cyp11a1","Star"),order=T,reduction = "tsne",raster=FALSE)

# SC 1 9 13
DotPlot(object = sce.big, features = c("Star","Col1a1","Col1a2","Mgp"))
DotPlot(object =all, features = c("Col1a1","Col1a2","Mgp","Csrp1","Acta2","Actg2","Cnn1","Myh11","Tagln"))

DotPlot(object = sce.big, features = c("Amhr2","Wt1","Foxl2","Inha","Gja1","Cyp17a1","Cyp11a1","Star","Hsd3b1","Lhcgr"))

FeaturePlot(object = sce.big, features = c("Star"),order = T,reduction = "tsne")
FeaturePlot(object = sce.big, features = c("Cyp11a1"),order = T,reduction = "umap")

DotPlot(object = sce.big, features = c("Tcf21","Star","Aldh1a1","Col1a1","Col1a2","Dcn","Bgn","Mgp","Mylk"))
FeaturePlot(object = sce.big, features = c("Tcf21","Star","Aldh1a1","Col1a1","Col1a2"),order = T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Nr2f2","Tcf21"))
FeaturePlot(object = sce.big, features = c("Tcf21","Nr2f2"),order=T,reduction = "tsne")


# GC  6 7 13 14
DotPlot(object =  sce.big, features = c("Amh","Cyp19a1","Inha","Bmpr2","Kitl","Bmpr1b","Inhba","Kctd14"))
FeaturePlot(object = sce.big, features =  c("Amh","Cyp19a1","Inha","Bmpr2","Kitl"),order=T,reduction = "umap")

FeaturePlot(object = sce.big, features =  c("Camta1"),order=T,reduction = "tsne")


# FBC 0 2 3 4  10 17 18 
DotPlot(object = sce.big, features = c("Col1a1","Pdgfra","Col1a2","Col3a1","Col4a1","Pi16",
                                       "Col15a1","Fsp1","Ddr2","Map2"))
FeaturePlot(object = sce.big, features =  c("Col1a1","Pdgfra","Col1a2","Col3a1","Col4a1"),
            order=T,reduction = "umap")
FeaturePlot(object = sce.big, features =  c("Col1a1"),order=T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Cyp17a1","Cyp11a1","Star","Hsd3b1","Lhcgr","Col1a1","Pdgfra","Col1a2","Col3a1"))



# EDC/BVEDC 5 16

DotPlot(object = sce.big, features = c("Aplnr","Cdh5","Nos3","Cldn11","Wt1","Pecam1","Vwf"))
FeaturePlot(object = sce.big, features =  c("Cdh5"),order=T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Cdh5","Pecam1","Vwf"))

# SMC 3 4 17 18
DotPlot(object = sce.big, features = c("Csrp1","Acta2","Actg2","Cnn1","Myh11","Tagln"))
FeaturePlot(object = sce.big, features =  c("Actg2"),order=T,reduction = "tsne")
DotPlot(object = all, features = c("Csrp1","Acta2","Actg2","Cnn1","Myh11","Tagln","Map2"))
FeaturePlot(object = all, features =  c("Myh11","Tagln","Map2"),order=T,reduction = "umap",raster=FALSE)


# EC 9
DotPlot(object = sce.big, features = c("Cldn11","Gja1","Gjb2","Lgr5","Vim","Krt8","Krt18","Krt19"))
FeaturePlot(object = sce.big, features =  c("Csrp1","Acta2","Actg2","Cnn1"),order=T,reduction = "umap")
FeaturePlot(object = sce.big, features =  c("Lgr5"),order=T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Epcam","Krt19","Prom1","Aldh1a1"))


# DC  21
DotPlot(object = sce.big, features = c("Itgax"))
FeaturePlot(object = sce.big, features =c("Itgax"),order=T,reduction = "tsne")
# M 12
DotPlot(object = sce.big, features = c("Cd68","Nos1","Nos2"))
FeaturePlot(object = sce.big, features =c("Cd68"),order=T,reduction = "tsne")
FeaturePlot(object = sce.big, features =c("Cd68"),order=T,reduction = "tsne")

DotPlot(object = sce.big, features = c("Cd68","Cd163","Cd14"))

# B 24
DotPlot(object = sce.big, features = c("Cd19","Pax5","Cd79a","Ms4a"))
FeaturePlot(object = sce.big, features =c("Pax5"),order=T,reduction = "tsne")

DotPlot(object = sce.big, features = c("Cd19","Tcl1a"))
DotPlot(object = sce.big, features = c("Cd4","Cd25"))

# LC(粒黄体) 11 15
DotPlot(object = sce.big, features = c("Lhcgr","Pgr","Ptgs2","Star","Cyp11a1"))
FeaturePlot(object = sce.big, features =c("Hsd3b1","Pgr","Npy","Pdgfra","Nppa"),order=T,reduction = "tsne")
FeaturePlot(object = sce.big, features =c("Lhcgr"),order=T,reduction = "tsne")


# GLC 粒细胞 25
DotPlot(object = sce.big, features = c("Mpo","Elane","Cebpb","Cxcr2","Csf3r","Itgam","Fpr1","Fpr2"))
FeaturePlot(object = sce.big, features =c("Csf3r"),order=T,reduction = "tsne")

# NK 23
DotPlot(object = sce.big, features = c("Nkg7","Tbx21","Ccl5","Klrb1c"))
FeaturePlot(object = sce.big, features =c("Tbx21"),order=T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Cx3cr1"))
# T 19 22
FeaturePlot(object = sce.big, features =c("Cd4","Trac"),order=T,reduction = "umap")
FeaturePlot(object = sce.big, features =c("Trac"),order=T,reduction = "tsne")
DotPlot(object = sce.big, features = c("Cd3d","Cd3e","Cd8a"))

DotPlot(object = sce.big, features = c("Tbx21"))
# NKT 19 22
DotPlot(object = sce.big, features = c("Il17rb","Cd3e","Cd8b1","Ly6c2","Trbc2",
                                       "Trac","Cd3g","Klrd1","Nkg7"))

DotPlot(object = sce.big, features = c("Cx3cr1","Cd14"))
FeaturePlot(object = sce.big, features =c("Nkg7"),order=T,reduction = "tsne")


# NC 4 5
FeaturePlot(object = sce.big, features =c("Map2"),order=T,reduction = "tsne")

DotPlot(object = sce.big, features = c("Map2"))

DotPlot(object = sce.big, features = c("Cyp17a1","Cyp11a1","Star",
                                       "Col1a1","Col1a2","Mgp",
                                       "Amh","Cyp19a1","Inha",
                                       "Col3a1","Col4a1",
                                       "Cdh5","Pecam1","Vwf",
                                       "Cldn11","Gja1","Gjb2",
                                       "Itgax",
                                       "Cd68","Nos1","Nos2",
                                       "Cd19","Pax5","Cd79a",
                                       "Elane","Cebpb","Cxcr2",
                                       "Nkg7","Tbx21","Ccl5",
                                       "Cd4","Trac"))


#FBC   0  2 10 11 17 18  
#EDC/BVEDC 5 16 20     
#EC  9
#NC  4 
#DC 21 
#T 19  
#NK  23    
#ILC  22 
#M 12
#B  24
#GLC 25
#STC 13
#GC 6 7 14
#TCC 1 8 11 15 26
#SMC 3 

clustertype<-c('0'='Fibroblasts',
               '1'='Theca cells',
               '2'='Fibroblasts',
               '3'="Smooth muscle cells",
               '4'='Neuron',
               '5'='Endothelial cells',
               '6'='Granulosa cells',
               '7'='Granulosa cells',
               '8'='Theca cells',
               '9'="Epithelial cells",
               '10'='Fibroblasts',
               '11'='Theca cells',
               '12'='Macrophages',
               '13'='Stem cells',
               '14'='Granulosa cells',
               '15'='Theca cells',
               '16'='Endothelial cells',
               '17'='Fibroblasts',
               '18'='Fibroblasts',
               '19'='T cells',
               '20'='Endothelial cells',
               "21"="Dendritic cells",
               "22"="Innate lymphoid cells",
               "23"="Natural killer cells",
               "24"="B cells",
               "25"="Granulocytes",
               "26"="Theca cells"
)
sce.big@meta.data$anno_celltype<-unname(unname(clustertype[sce.big@meta.data$seurat_clusters]))

new.cluster.ids <- c('Fibroblasts',
                     'Theca cells',
                     'Fibroblasts',
                     "Smooth muscle cells",
                     'Neuron',
                     'Endothelial cells',
                     'Granulosa cells',
                     'Granulosa cells',
                     'Theca cells',
                     "Epithelial cells",
                     'Fibroblasts',
                     'Theca cells',
                     'Macrophages',
                     'Stem cells',
                     'Granulosa cells',
                     'Theca cells',
                     'Endothelial cells',
                     'Fibroblasts',
                     'Fibroblasts',
                     'T cells',
                     'Endothelial cells',
                     "Dendritic cells",
                     "Innate lymphoid cells",
                     "Natural killer cells",
                     "B cells",
                     "Granulocytes",
                     "Theca cells")
names(x=new.cluster.ids)=levels(x = sce.big)
sce.big <- RenameIdents(object =sce.big, new.cluster.ids)
sce.big@meta.data$orig.ident
unique(c('Fibroblasts',
         'Theca cells',
         'Fibroblasts',
         "Smooth muscle cells",
         'Neuron',
         'Endothelial cells',
         'Granulosa cells',
         'Granulosa cells',
         'Theca cells',
         "Epithelial cells",
         'Fibroblasts',
         'Theca cells',
         'Macrophages',
         'Stem cells',
         'Granulosa cells',
         'Theca cells',
         'Endothelial cells',
         'Fibroblasts',
         'Fibroblasts',
         'T cells',
         'Endothelial cells',
         "Dendritic cells",
         "Innate lymphoid cells",
         "Natural killer cells",
         "B cells",
         "Granulocytes",
         "Theca cells"))
sce.big@active.ident <- factor(sce.big@active.ident,
                               levels = c("Stem cells","Granulosa cells","Theca cells","Fibroblasts","Smooth muscle cells",
                                          "Neuron","Endothelial cells","Epithelial cells",
                                          "Macrophages" ,"Dendritic cells","Granulocytes", "Innate lymphoid cells",
                                          "Natural killer cells" ,"B cells","T cells"))
sce.big@meta.data$Diets <- paste0(substr(colnames(sce.big),1,4),"RD")
sce.big@meta.data$Diets <- factor(sce.big@meta.data$Diets,levels = c("M6-BRD","M6-WRD","M6-PRD",
                                                                     "M9-BRD","M9-WRD","M9-PRD"))
library(ggsci)
DimPlot(object = sce.big, reduction = "tsne",label = T,label.size = 6)
DimPlot(object = sce.big, reduction = "tsne",label = F,label.size = 6,group.by = "Diets",raster=FALSE)+ 
  scale_color_igv()
DimPlot(object = sce.big, reduction = "tsne",label = F,raster=FALSE,
        label.size = 4,cols = c("#6686c6","#a0c1db","#748AA6","#FAE48B","gray","#ce8f5c",	
                                "#7bac80","#75a5c0","#b5181a","#b72d9e",
                                "#e4cccb","#f6f0eb","#e8c755","#d5d456","#cfe74f"))

sce.big@meta.dat$stage <- substr(colnames(sce.big),1,2)
M6 <- subset(sce.big,cells=rownames(sce.big@meta.data[grep("M6",sce.big@meta.data$Diets),]))
length(rownames(sce.big@meta.data[grep("M6",sce.big@meta.data$Diets),]))
length(rownames(sce.big@meta.data[grep("M9",sce.big@meta.data$Diets),]))
colnames(sce.big)%>%length()
M9 <- subset(sce.big,cells=rownames(sce.big@meta.data[grep("M9",sce.big@meta.data$Diets),]))
saveRDS(sce.big,"./all_obj.rds")
saveRDS(M6,"./M6_obj.rds")
saveRDS(M9,"./M9_obj.rds")
rm(sce.big,M6)
gc()
####################################################




############################################## Rename clusters ########################################
clustertype<-c('0'='Stroma cells',
               '1'='Theca cells',
               '2'='Stroma cells',
               '3'="Stroma cells",
               '4'='Stroma cells',
               '5'='Endothelial cells',
               '6'='Granulosa cells',
               '7'='Granulosa cells',
               '8'='Theca cells',
               '9'="Epithelial cells",
               '10'='Stroma cells',
               '11'='Theca cells',
               '12'='Macrophages',
               '13'='Stem cells',
               '14'='Granulosa cells',
               '15'='Theca cells',
               '16'='Endothelial cells',
               '17'='Stroma cells',
               '18'='Stroma cells',
               '19'='T cells',
               '20'='Endothelial cells',
               "21"="Dendritic cells",
               "22"="Innate lymphoid cells",
               "23"="Natural killer cells",
               "24"="B cells",
               "25"="Granulocytes",
               "26"="Theca cells"
)
all@meta.data$anno_celltype<-unname(unname(clustertype[all@meta.data$seurat_clusters]))
all@meta.data$anno_celltype <- factor(all@meta.data$anno_celltype,
                                     levels = c("Stem cells","Granulosa cells","Theca cells","Stroma cells","Endothelial cells",
                                                "Epithelial cells","Macrophages" ,"Dendritic cells","Granulocytes", "Innate lymphoid cells",
                                                "Natural killer cells" ,"B cells","T cells"))
N <- all@meta.data$anno_celltype
names(N) <- rownames(all@meta.data)
N
all@active.ident <- N
all@meta.data$Diets <- paste0(substr(colnames(all),1,4),"RD")
all@meta.data$Diets <- paste0(substr(colnames(all),1,4),"RD")
all@meta.data$Diets <- factor(all@meta.data$Diets,levels = c("M6-BRD","M6-WRD","M6-PRD",
                                                                     "M9-BRD","M9-WRD","M9-PRD"))
library(ggsci)
DimPlot(object = all, reduction = "tsne",label = T,label.size = 6)
DimPlot(object = all, reduction = "tsne",label = F,label.size = 6,group.by = "Diets",raster=FALSE)+scale_color_simpsons()
DimPlot(object = all, reduction = "tsne",label = F,raster=FALSE,
        label.size = 4)+scale_color_simpsons()
DimPlot(object = all, reduction = "umap",label = F,raster=FALSE,
        label.size = 4)+scale_color_simpsons()

DimPlot(object = all, reduction = "umap",label = T,raster=FALSE,group.by="seurat_clusters",
        label.size = 4)
all@meta.dat$stage <- substr(colnames(all),1,2)
######################################################################################################
all <- readRDS("./all_obj.rds")
all@meta.data$group_age_diet <- paste0(all@meta.data$Diets,"_",all@meta.data$anno_celltype)
all@meta.data$group_age <- paste0(substr(all@meta.data$Diets,1,2),"_",all@meta.data$anno_celltype)

abbreviation<-c('0'='SC',
               '1'='TCC',
               '2'='SC',
               '3'="SC",
               '4'='SC',
               '5'='EC',
               '6'='GC',
               '7'='GC',
               '8'='TCC',
               '9'="EC",
               '10'='SC',
               '11'='TCC',
               '12'='M',
               '13'='STC',
               '14'='GC',
               '15'='TCC',
               '16'='EDC',
               '17'='SC',
               '18'='SC',
               '19'='T',
               '20'='EDC',
               "21"="DC",
               "22"="ILC",
               "23"="NK",
               "24"="B",
               "25"="GLC",
               "26"="TCC"
)
all@meta.data$abbreviation<-unname(unname(abbreviation[all@meta.data$seurat_clusters]))

cell <- c("STC","GC","TCC","SC","EDC","EC","M","DC","GLC","ILC","NK","B","T")
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(all,cells=rownames(all@meta.data[all$abbreviation==cell[x],]))
  saveRDS(tmp, paste0(cell[x],".rds"))
})

######################## subset BRD ############################
meta <- all@meta.data[grep("BRD",all$Diets),]
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(all,cells=rownames(meta[meta$abbreviation==cell[x],]))
  saveRDS(tmp, paste0("BRD_",cell[x],".rds"))
})

######################## subset WRD ############################
meta <- all@meta.data[grep("WRD",all$Diets),]
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(all,cells=rownames(meta[meta$abbreviation==cell[x],]))
  saveRDS(tmp, paste0("WRD_",cell[x],".rds"))
})

######################## subset PRD ############################
meta <- all@meta.data[grep("PRD",all$Diets),]
rds <- purrr::map(1:length(cell),function(x){
  tmp <- subset(all,cells=rownames(meta[meta$abbreviation==cell[x],]))
  saveRDS(tmp, paste0("PRD_",cell[x],".rds"))
})


counts <- as.data.frame(table(all@meta.data$abbreviation,all@meta.data$orig.ident))

length(colnames(all))
# 145144

a <- table(all@meta.data$orig.ident)%>%as.data.frame()
a
counts$all <- c(rep(a$Freq[1],13),rep(a$Freq[2],13),rep(a$Freq[3],13),
                rep(a$Freq[4],13),rep(a$Freq[5],13),rep(a$Freq[6],13),
                rep(a$Freq[7],13),rep(a$Freq[8],13),rep(a$Freq[9],13),
                rep(a$Freq[10],13),rep(a$Freq[11],13),rep(a$Freq[12],13),
                rep(a$Freq[13],13),rep(a$Freq[14],13),rep(a$Freq[15],13),
                rep(a$Freq[16],13),rep(a$Freq[17],13))
counts$Frac <- counts$Freq/counts$all
head(counts)
counts$Diets <- paste0(strsplit2(counts$Var2,"-")[,2],"RD")
counts$Diets <- factor(counts$Diets,levels = c("BRD","WRD","PRD"))
counts$Stage <- strsplit2(counts$Var2,"-")[,1]
counts$Celltype <- factor(counts$Var1,levels=c("STC","GC","TCC","SC","EDC","EC","M","DC","GLC","ILC","NK","B","T"))
ggplot(counts,aes(x=Diets,y=Frac*100,fill=Stage))+#stat_boxplot(geom = "errorbar",width=0.6)+
  geom_boxplot(width=0.6)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  xlab("Celltype")+ylab("Cell Component Fractions %")+
  stat_compare_means(aes(group=Stage),label = "p.signif",method = "t.test")+ 
  facet_wrap( ~Celltype, ncol=5,scales = "free_y")+scale_fill_manual(values=c("#E7AE8F","#BF5217"))+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())

ggplot(counts[counts$Diets=="WRD",],aes(x=Stage,y=Frac*100,fill=Stage))+#stat_boxplot(geom = "errorbar",width=0.6)+
  geom_boxplot(width=0.6)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  xlab("Celltype")+ylab("Cell Component Fractions %")+
  geom_signif(comparisons = list(c("M6","M9")),map_signif_level = T,step_increase = 0.1,
              textsize =5,tip_length = 0,vjust = 0.2)+ 
  facet_wrap( ~Celltype, ncol=5,scales = "free_y")+scale_fill_manual(values=c("#E7AE8F","#BF5217"))+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())

all <- readRDS("../output/all_obj.rds")
all$abbreviation2 <- gsub("SC","SC&TC",all$abbreviation)%>%gsub("TCC","SC&TC",.)

all@meta.data$age <- substr(all$group_age,1,2)
DimPlot(object = all, reduction = "tsne",label = F,
        group.by="age",raster=FALSE,pt.size = 0.05) 
all$abbreviation2 <- factor(all$abbreviation2,levels = c("STC",
                                                         "GC","SC&TC","EC","EDC",
                                                         "GLC","M","DC","B","NK","T","ILC"))
all$Group <- paste0(all$age,"-",all$abbreviation2)


DotPlot(object = all,group.by = "Group",
        cols = c("white","#2B589B"), 
        features = c("Amh","Cyp19a1","Inha",
                     "Cyp11a1","Star",
                     "Col1a1","Col1a2","Mgp","Col4a1",
                     "Cdh5","Pecam1","Vwf",
                     "Cebpb","Cxcr2","Csf3r","Itgam",
                     "Cd68","Clec9a","Itgax","Itgae",
                     "Cd19","Pax5","Cd79a",
                     "Nkg7","Tbx21","Ccl5",
                     "Cd4","Trac","Rorc"))+
  theme_bw()+scale_y_discrete(limits=c("M6-STC","M9-STC","M6-GC","M9-GC",
                                       "M6-SC&TC","M9-SC&TC","M6-EC","M9-EC",
                                       "M6-EDC","M9-EDC","M6-GLC","M9-GLC",
                                       "M6-M","M9-M","M6-DC","M9-DC","M6-B","M9-B",
                                       "M6-NK","M9-NK","M6-T","M9-T","M6-ILC","M9-ILC"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.title =element_text(size = 16),
        axis.text =element_text(size = 16, color = 'black'))+
  theme(axis.title =element_text(size = 16),legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text =element_text(size = 16, color = 'black'))+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust = 0.25,hjust = 1))+coord_flip()


DotPlot(object = all,group.by = "Group", 
        features = c("Foxo3"))+
  theme_bw()+scale_y_discrete(limits=c("M6-STC","M9-STC","M6-GC","M9-GC",
                                       "M6-SC&TC","M9-SC&TC","M6-EC","M9-EC",
                                       "M6-EDC","M9-EDC","M6-GLC","M9-GLC",
                                       "M6-M","M9-M","M6-DC","M9-DC","M6-B","M9-B",
                                       "M6-NK","M9-NK","M6-T","M9-T","M6-ILC","M9-ILC"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.title =element_text(size = 16),
        axis.text =element_text(size = 16, color = 'black'))+
  theme(axis.title =element_text(size = 16),legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text =element_text(size = 16, color = 'black'))+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust = 0.25,hjust = 1))+coord_flip()
all$group_age_diet
DotPlot(object = all,group.by = "group_age_diet", 
        features = c("Foxo3"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+coord_flip()



a=all$abbreviation2
names(a)
all <- readRDS("../output/all_obj.rds")
all$abbreviation2 <- gsub("SC","SC&TC",all$abbreviation)%>%gsub("TCC","SC&TC",.)
unique(all$abbreviation2)
all$abbreviation2 <- factor(all$abbreviation2,levels = c("STC",
                                                         "GC","SC&TC","EC","EDC",
                                                         "GLC","M","DC","B","NK","T","ILC"))
all@meta.data$age <- substr(all$group_age,1,2)
all$Group <- paste0(all$age,"-",all$abbreviation2)
all@active.ident <- all$abbreviation2
all$diets <- substr(all$Diets,4,6)%>%factor(.,levels=c("BRD","WRD","PRD"))
DimPlot(object = all, reduction = "tsne",label = F,order =F,seed = 123,
        cols=c("#E33D30","#F9C52F","#307CEE"),
        group.by="diets",raster=FALSE,pt.size = 0.01) 

DimPlot(object = all, reduction = "tsne",label = F,order =F,seed = 12,
        cols=c("#FF9800","gray"),
        group.by="age",raster=FALSE,pt.size = 0.01)

all_markers <- FindAllMarkers(object =all,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers,"../output/all_markers.csv")

all_top50  <- all_markers  %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)
st.data <- prepareDataFromscRNA(object = all,
                                diffData = all_top50,
                                showAverage = TRUE)


#对每个cluster进行富集分析，这里采用GO富集分析，有需要的可以选择kegg富集分析

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 30,
                        seed = 5201314)

#挑选需要展示的marker基因
markGenes = unique(all_top50$gene)[sample(1:length(unique(all_top50$gene)),40,replace = F)]

#绘制cluster基因表达折线图

visCluster(object = st.data,plot.type = "line")
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:9),
           go.col = rep(jjAnno::useMyCol("stallion",n = 9),each = 5),
           add.bar = T)








prop_test <- sc_utils(all)
prop_test <- permutation_test(
  prop_test, cluster_identity = "abbreviation2",
  sample_1 = "M6", sample_2 = "M9",
  sample_identity = "age"
)
?permutation_plot

permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

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
		plot_data[, clusters := fct_reorder(factor(clusters), desc(obs_log2FD))]
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
#########################################################################################
#M6 <- subset(all,cells=rownames(all@meta.data[grep("M6",all@meta.data$Diets),]))
#length(rownames(all@meta.data[grep("M6",all@meta.data$Diets),]))
#length(rownames(all@meta.data[grep("M9",all@meta.data$Diets),]))
#colnames(all)%>%length()
M9 <- subset(all,cells=rownames(all@meta.data[grep("M9",all@meta.data$Diets),]))
saveRDS(all,"./all_obj.rds")
rm(all)
gc()
# saveRDS(M6,"./M6_obj.rds")
saveRDS(M9,"./M9_obj.rds")






############################################## M9 #####################################################
library(ggsci)
M9 <- readRDS("./M9_obj.rds")
M9$abbreviation2 <- gsub("SC","SC&TC",M9$abbreviation)%>%gsub("TCC","SC&TC",.)
############################# SASP  GSVA ###########################################
aging_gene_sets <- read.delim("~/scRNA_aging/output/aging_gene_sets.txt")
head(aging_gene_sets)
gene1 <- as.list(aging_gene_sets[,1])[1:34]
gene2 <- as.list(aging_gene_sets[,2])
gene3 <- as.list(aging_gene_sets[,3])[1:40]
gene4 <- as.list(aging_gene_sets[,4])[1:27][-c(20,8)]






################################## Seurat Scores ###########################################
M9 <- AddModuleScore(object = M9,features = gene1,ctrl = 100,name = "Cell_Cycle")
M9 <- AddModuleScore(object = M9,features = gene2,ctrl = 100,name = "DNA_Repair")
M9 <- AddModuleScore(object = M9,features = gene3,ctrl = 100,name = "SASP")
M9 <- AddModuleScore(object = M9,features = gene4,ctrl = 100,name = "Inflammation")
M9$Group <- paste0(substr(M9$Diets,4,6),"-",M9$abbreviation2)

data1 <- FetchData(M9,vars = c("Group","Cell_Cycle1"))
data2 <- FetchData(M9,vars = c("Group","DNA_Repair1"))
data3 <- FetchData(M9,vars = c("Group","SASP1"))
data4 <- FetchData(M9,vars = c("Group","Inflammation1"))
head(data1)
cells <- rownames(data1)

data <- data.frame(Diets=strsplit2(data1$Group,"-")[,1],
                   Celltype=strsplit2(data1$Group,"-")[,2],Scores=data1$Cell_Cycle1,
                    row.names = cells)
ggplot(data, aes(x = Diets, y =log2(Scores+1),fill=Diets)) +
  geom_violin(aes(fill = Diets)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("BRD","WRD"),c("WRD","PRD"),c("BRD","PRD")),
    map_signif_level = T,step_increase = 0.1,size = 1,textsize = 7,
    tip_length = 0,vjust = 0.2)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1),axis.line.y=element_line(size=1))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))+facet_wrap(~Celltype)




cells <- rownames(data3)

data <- data.frame(Diets=strsplit2(data3$Group,"-")[,1],
                   Celltype=strsplit2(data3$Group,"-")[,2],Scores=data3$SASP1,
                   row.names = cells)
ggplot(data, aes(x = Diets, y =log2(Scores+1),fill=Diets)) +
  geom_violin(aes(fill = Diets)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("BRD","WRD"),c("WRD","PRD"),c("BRD","PRD")),
    map_signif_level = T,step_increase = 0.1,size = 1,textsize = 7,
    tip_length = 0,vjust = 0.2)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1),axis.line.y=element_line(size=1))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))+facet_wrap(~Celltype)


cells <- rownames(data4)

data <- data.frame(Diets=strsplit2(data4$Group,"-")[,1],
                   Celltype=strsplit2(data4$Group,"-")[,2],Scores=data4$Inflammation1,
                   row.names = cells)
ggplot(data, aes(x = Diets, y =log2(Scores+1),fill=Diets)) +
  geom_violin(aes(fill = Diets)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(
    c("BRD","WRD"),c("WRD","PRD"),c("BRD","PRD")),
    map_signif_level = T,step_increase = 0.1,size = 1,textsize = 7,
    tip_length = 0,vjust = 0.2)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24, color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=1),axis.line.y=element_line(size=1))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 24))+facet_wrap(~Celltype)

####################################################################################
M9 <- RunUMAP(object = M9, dims = 1:40)
M9 <- RunTSNE(object = M9, dims = 1:40)
DimPlot(object = M9, reduction = "tsne",label = F,label.size = 6,group.by = "abbreviation")+scale_color_simpsons()
DimPlot(object = M9, reduction = "umap",label = F,label.size = 6,group.by = "abbreviation")+scale_color_simpsons()
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())
DimPlot(object = M9, reduction = "umap",label = F,label.size = 6,group.by = "abbreviation2")+scale_color_simpsons()


## Cell fraction
devtools::install_local("~/scProportionTest.tar.gz")
library("scProportionTest")

library("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = rownames(M9), #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)

head(MtoH)
data <- as.matrix(M9@assays$RNA@data) # Notice: 这里应用的是空间数据，常规转录组数据提取用sp1@assays$RNA@data
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
data[1:5,1:5]
rownames(data) <- NULL
write.table(data, './cellphonedb_M9_count.txt', sep='\t', quote=F)
saveRDS(data,"./M9_count_data.rds")


meta_data <- cbind(rownames(M9@meta.data), M9@meta.data[,'abbreviation2', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M9_meta_data.txt', sep='\t', quote=F)


############################################## M6 #####################################################
M6 <- readRDS("./M6_obj.rds")
unique(M6$abbreviation)
M6$abbreviation2 <- gsub("FBC","SCTC",M6$abbreviation)%>%gsub("TCC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- gsub("SMC","SCTC",M6$abbreviation2)%>%gsub("NRC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- factor(M6$abbreviation2,levels = c("STC",
                                                         "GC","SCTC","EC","EDC",
                                                         "GLC","M","DC","B","NK","T","ILC"))
unique(M6$abbreviation2)
M6@meta.data$age <- substr(M6$group_age,1,2)
M6$Group <- paste0(M6$age,"-",M6$abbreviation2)
M6@active.ident <- M6$abbreviation2

STGC <- subset(M6,idents=c("STC","GC"))
head(STGC@meta.data)

M6_select <- subset(M6,cells = sample(colnames(M6),40000))
M6_markers <- FindAllMarkers(object =M6_select,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(M6_markers,"../output/M6_markers.csv")
saveRDS(M6_select,"../output/M6_select.rds")
devtools::load_all("~/ClusterGVis-main/")

devtools::install_local("/home/tuyx/jjAnno.tar.gz")
devtools::load_all("/home/tuyx/magick_2.0/magick/")
BiocManager::install("igraph")
library(ClusterGVis)
Idents(M6)
M6_top50  <- M6_markers  %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 50, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = M6,
                                diffData = M6_top50,
                                showAverage = TRUE)

saveRDS(st.data,"./st.data.rds")

#对每个cluster进行富集分析，这里采用GO富集分析，有需要的可以选择kegg富集分析

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.5,
                        topn = 30,
                        seed = 5201314)

#挑选需要展示的marker基因
markGenes = unique(M6_top50$gene)[sample(1:length(unique(M6_top50$gene)),40,replace = F)]
markGenes =c("Amh","Cyp19a1","Inha",
             "Cyp11a1","Star",
             "Col1a1","Col1a2","Mgp","Col4a1",
             "Cdh5","Pecam1","Vwf",
             "Cebpb","Cxcr2","Csf3r","Itgam",
             "Cd68","Clec9a","Itgax","Itgae",
             "Cd19","Pax5","Cd79a",
             "Nkg7","Tbx21","Ccl5",
             "Cd4","Trac","Rorc"))
#绘制cluster基因表达折线图

visCluster(object = st.data,plot.type = "line")
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:9),
           go.col = rep(jjAnno::useMyCol("stM6ion",n = 9),each = 5),
           add.bar = T)







DimPlot(object = M6, reduction = "tsne",label = F,label.size = 6,group.by = "abbreviation2")+scale_color_simpsons()+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())

data <- as.matrix(M6@assays$RNA@data) # Notice: 这里应用的是空间数据，常规转录组数据提取用sp1@assays$RNA@data
rm(M6)
gc()
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
data[1:5,1:5]
rownames(data) <- NULL
saveRDS(data,"./M6_count_data.rds")
data <- readRDS("./M6_count_data.rds")
write.table(data, './cellphonedb_M6_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
rm(data)
gc()

meta_data <- cbind(rownames(M6@meta.data), M6@meta.data[,'abbreviation2', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M6_meta_data.txt', sep='\t', quote=F)
write.table(M6_meta, './M6_meta.txt',col.names=T,row.names=T, sep='\t', quote=F)
###################### 细胞数目统计 6M #######################
M6 <- readRDS("./M6_obj.rds")
M6$abbreviation2 <- gsub("SC","SC&TC",M6$abbreviation)%>%gsub("TCC","SC&TC",.)
M6_meta <- M6@meta.data
rm(M6)
gc()
M6_meta$sample <- strsplit2(rownames(M6_meta),"_")[,1]
data1 <- as.data.frame(table(M6_meta$sample))

###################### 细胞数目统计 9M #######################
M9 <- readRDS("./M9_obj.rds")
M9$abbreviation2 <- gsub("SC","SC&TC",M9$abbreviation)%>%gsub("TCC","SC&TC",.)
M9_meta <- M9@meta.data
rm(M9)
gc()
M9_meta$sample <- strsplit2(rownames(M9_meta),"_")[,1]
data2 <- as.data.frame(table(M9_meta$sample))
data2
write.table(M9_meta, './M9_meta.txt',col.names=T,row.names=T, sep='\t', quote=F)

data <- rbind(data1,data2)
dd <- strsplit2(data$Var1,"-")[,1:2]%>%as.data.frame()
colnames(dd) <- c("Stage","Diets")
dd$Diets <- paste0(dd$Diets,"RD")
data <- cbind(data,dd)
data$Diets <- factor(data$Diets,levels = c("BRD","WRD","PRD"))
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
    scale_fill_manual(values = c("#e20612","#ffd401","#00b0eb"),name="")+
  labs(x="Stages",y="Cell number")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+
  annotate(geom = 'text', label="N=8196", x=0.8, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="n=9709", x=1, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="n=8801", x=1.2, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="n=8592", x=1.8, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="n=7315", x=2.0, y=2000, angle=90,size=5)+
  annotate(geom = 'text', label="n=8652", x=2.2, y=2000, angle=90,size=5)
  




meta <- M9@meta.data[grep("BRD",M9$Diets),]
M9_BRD <- subset(M9,cells=rownames(meta))
meta <- M9@meta.data[grep("WRD",M9$Diets),]
M9_WRD <- subset(M9,cells=rownames(meta))
meta <- M9@meta.data[grep("PRD",M9$Diets),]
M9_PRD <- subset(M9,cells=rownames(meta))
saveRDS(M9_BRD,"./M9_BRD.rds")
saveRDS(M9_WRD,"./M9_WRD.rds")
saveRDS(M9_PRD,"./M9_PRD.rds")


M9_BRD <- readRDS("./M9_BRD.rds")
data <- as.matrix(M9_BRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_M9_BRD_count.txt',col.names=T,row.names=F, sep='\t', quote=F)

rm(M9_BRD,data)


meta_data <- cbind(rownames(M9_BRD@meta.data), M9_BRD@meta.data[,'abbreviation2', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M9_BRD_meta_data.txt', sep='\t', quote=F)
rm(M9_BRD,data)
gc()

M9_WRD <- readRDS("./M9_WRD.rds")
data <- as.matrix(M9_WRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
write.table(data, 'cellphonedb_M9_WRD_count.txt',col.names=T,row.names=F, sep='\t', quote=F)


meta_data <- cbind(rownames(M9_WRD@meta.data), M9_WRD@meta.data[,'abbreviation2', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M9_WRD_meta_data.txt', sep='\t', quote=F)
rm(M9_WRD,data)
gc()

M9_PRD <- readRDS("./M9_PRD.rds")
data <- as.matrix(M9_PRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
write.table(data, 'cellphonedb_M9_PRD_count.txt',col.names=T,row.names=F,  sep='\t', quote=F)

meta_data <- cbind(rownames(M9_PRD@meta.data), M9_PRD@meta.data[,'abbreviation2', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M9_PRD_meta_data.txt',sep='\t', quote=F)
rm(M9_PRD,data)
gc()

########################## 差异基因 ##############################
########################## GC ####################################
BRD_GC <- readRDS("./BRD_GC.rds")
df1=as.data.frame(BRD_GC@active.ident)
idents <- gsub("Granulosa cells","GC",BRD_GC$group_age)
names(idents) <- rownames(BRD_GC@meta.data)
BRD_GC@active.ident<-factor(idents)
BRD_GC_markers <- FindMarkers(BRD_GC, ident.1 = "M9_GC", ident.2 = "M6_GC", min.pct = 0.25)
BRD_GC_markers$group <- 0
BRD_GC_markers[BRD_GC_markers$avg_log2FC>0.25,]$group <- "9M"
BRD_GC_markers[BRD_GC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_GC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_GC_markers$group))
data1$diets <- "BRD"

WRD_GC <- readRDS("./WRD_GC.rds")
df1=as.data.frame(WRD_GC@active.ident)
idents <- gsub("Granulosa cells","GC",WRD_GC$group_age)
names(idents) <- rownames(WRD_GC@meta.data)
WRD_GC@active.ident<-factor(idents)
WRD_GC_markers <- FindMarkers(WRD_GC, ident.1 = "M9_GC", ident.2 = "M6_GC", min.pct = 0.25)
WRD_GC_markers$group <- 0
WRD_GC_markers[WRD_GC_markers$avg_log2FC > 0.25,]$group <- "9M"
WRD_GC_markers[WRD_GC_markers$avg_log2FC< -0.25,]$group <- "6M"
WRD_GC_markers$diets <- "WRD"
data2 <- as.data.frame(table(WRD_GC_markers$group))
data2$diets <- "WRD"

PRD_GC <- readRDS("./PRD_GC.rds")
df1=as.data.frame(PRD_GC@active.ident)
idents <- gsub("Granulosa cells","GC",PRD_GC$group_age)
names(idents) <- rownames(PRD_GC@meta.data)
PRD_GC@active.ident<-factor(idents)
PRD_GC_markers <- FindMarkers(PRD_GC, ident.1 = "M9_GC", ident.2 = "M6_GC", min.pct = 0.25)
PRD_GC_markers$group <- 0
PRD_GC_markers[PRD_GC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_GC_markers[PRD_GC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_GC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_GC_markers$group))
data3$diets <- "PRD"
data <- rbind(data1,data2)%>%rbind(.,data3)
data$celltype <- "GC"
d1 <- data
a=setdiff(rownames(BRD_GC_markers),c(rownames(WRD_GC_markers),rownames(PRD_GC_markers)))
BP1<- clusterProfiler::simplify(enrichGO(gene=a,keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

View(BP1@result)

b=intersect(rownames(BRD_GC_markers),rownames(WRD_GC_markers))%>%intersect(.,rownames(PRD_GC_markers))
BP2<- clusterProfiler::simplify(enrichGO(gene=b,keyType = "SYMBOL",
                                        OrgDb= org.Mm.eg.db, ont="BP",
                                        pAdjustMethod="BH",pvalueCutoff=0.05,
                                        qvalueCutoff=0.05))

########################## M ####################################
BRD_M <- readRDS("./BRD_M.rds")
df1=as.data.frame(BRD_M@active.ident)
idents <- gsub("Macrophages","M",BRD_M$group_age)
names(idents) <- rownames(BRD_M@meta.data)
BRD_M@active.ident<-factor(idents)
BRD_M_markers <- FindMarkers(BRD_M, ident.1 = "M9_M", ident.2 = "M6_M", min.pct = 0.25)
BRD_M_markers$group <- 0
BRD_M_markers[BRD_M_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_M_markers[BRD_M_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_M_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_M_markers$group))
data1$diets <- "BRD"

WRD_M <- readRDS("./WRD_M.rds")
df1=as.data.frame(WRD_M@active.ident)
idents <- gsub("Macrophages","M",WRD_M$group_age)
names(idents) <- rownames(WRD_M@meta.data)
WRD_M@active.ident<-factor(idents)
WRD_M_markers <- FindMarkers(WRD_M, ident.1 = "M9_M", ident.2 = "M6_M", min.pct = 0.25)
WRD_M_markers$group <- 0
WRD_M_markers[WRD_M_markers$avg_log2FC > 0.25,]$group <- "9M"
WRD_M_markers[WRD_M_markers$avg_log2FC< -0.25,]$group <- "6M"
WRD_M_markers$diets <- "WRD"
data2 <- as.data.frame(table(WRD_M_markers$group))
data2$diets <- "WRD"

PRD_M <- readRDS("./PRD_M.rds")
df1=as.data.frame(PRD_M@active.ident)
idents <- gsub("Macrophages","M",PRD_M$group_age)
names(idents) <- rownames(PRD_M@meta.data)
PRD_M@active.ident<-factor(idents)
PRD_M_markers <- FindMarkers(PRD_M, ident.1 = "M9_M", ident.2 = "M6_M", min.pct = 0.25)
PRD_M_markers$group <- 0
PRD_M_markers[PRD_M_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_M_markers[PRD_M_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_M_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_M_markers$group))
data3$diets <- "PRD"
data <- rbind(data1,data2)%>%rbind(.,data3)
data$celltype <- "M"
d2 <- data


########################## STC ####################################
BRD_STC <- readRDS("./BRD_STC.rds")
df1=as.data.frame(BRD_STC@active.ident)
idents <- gsub("Stem cells","STC",BRD_STC$group_age)
names(idents) <- rownames(BRD_STC@meta.data)
BRD_STC@active.ident<-factor(idents)
BRD_STC_markers <- FindMarkers(BRD_STC, ident.1 = "M9_STC", ident.2 = "M6_STC", min.pct = 0.25)
BRD_STC_markers$group <- 0
BRD_STC_markers[BRD_STC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_STC_markers[BRD_STC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_STC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_STC_markers$group))
data1$diets <- "BRD"

WRD_STC <- readRDS("./WRD_STC.rds")
df1=as.data.frame(WRD_STC@active.ident)
idents <- gsub("Stem cells","STC",WRD_STC$group_age)
names(idents) <- rownames(WRD_STC@meta.data)
WRD_STC@active.ident<-factor(idents)
WRD_STC_markers <- FindMarkers(WRD_STC, ident.1 = "M9_STC", ident.2 = "M6_STC", min.pct = 0.25)
WRD_STC_markers$group <- 0
WRD_STC_markers[WRD_STC_markers$avg_log2FC > 0.25,]$group <- "9M"
WRD_STC_markers[WRD_STC_markers$avg_log2FC< -0.25,]$group <- "6M"
WRD_STC_markers$diets <- "WRD"
data2 <- as.data.frame(table(WRD_STC_markers$group))
data2$diets <- "WRD"

PRD_STC <- readRDS("./PRD_STC.rds")
df1=as.data.frame(PRD_STC@active.ident)
idents <- gsub("Stem cells","STC",PRD_STC$group_age)
names(idents) <- rownames(PRD_STC@meta.data)
PRD_STC@active.ident<-factor(idents)
PRD_STC_markers <- FindMarkers(PRD_STC, ident.1 = "M9_STC", ident.2 = "M6_STC", min.pct = 0.25)
PRD_STC_markers$group <- 0
PRD_STC_markers[PRD_STC_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_STC_markers[PRD_STC_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_STC_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_STC_markers$group))
data3$diets <- "PRD"
data <- rbind(data1,data2)%>%rbind(.,data3)
data$celltype <- "STC"
d3 <- data

########################## DC ####################################
BRD_DC <- readRDS("./BRD_DC.rds")
df1=as.data.frame(BRD_DC@active.ident)
idents <- gsub("Dendritic cells","DC",BRD_DC$group_age)
names(idents) <- rownames(BRD_DC@meta.data)
BRD_DC@active.ident<-factor(idents)
BRD_DC_markers <- FindMarkers(BRD_DC, ident.1 = "M9_DC", ident.2 = "M6_DC", min.pct = 0.25)
BRD_DC_markers$group <- 0
BRD_DC_markers[BRD_DC_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_DC_markers[BRD_DC_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_DC_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_DC_markers$group))
data1$diets <- "BRD"

WRD_DC <- readRDS("./WRD_DC.rds")
df1=as.data.frame(WRD_DC@active.ident)
idents <- gsub("Dendritic cells","DC",WRD_DC$group_age)
names(idents) <- rownames(WRD_DC@meta.data)
WRD_DC@active.ident<-factor(idents)
WRD_DC_markers <- FindMarkers(WRD_DC, ident.1 = "M9_DC", ident.2 = "M6_DC", min.pct = 0.25)
WRD_DC_markers$group <- 0
WRD_DC_markers[WRD_DC_markers$avg_log2FC > 0.25,]$group <- "9M"
WRD_DC_markers[WRD_DC_markers$avg_log2FC< -0.25,]$group <- "6M"
WRD_DC_markers$diets <- "WRD"
data2 <- as.data.frame(table(WRD_DC_markers$group))
data2$diets <- "WRD"

PRD_DC <- readRDS("./PRD_DC.rds")
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
data3$diets <- "PRD"
data <- rbind(data1,data2)%>%rbind(.,data3)
data$celltype <- "DC"
d4 <- data

########################## B ####################################
BRD_B <- readRDS("./BRD_B.rds")
df1=as.data.frame(BRD_B@active.ident)
idents <- gsub("B cells","B",BRD_B$group_age)
names(idents) <- rownames(BRD_B@meta.data)
BRD_B@active.ident<-factor(idents)
BRD_B_markers <- FindMarkers(BRD_B, ident.1 = "M9_B", ident.2 = "M6_B", min.pct = 0.25)
BRD_B_markers$group <- 0
BRD_B_markers[BRD_B_markers$avg_log2FC > 0.25,]$group <- "9M"
BRD_B_markers[BRD_B_markers$avg_log2FC< -0.25,]$group <- "6M"
BRD_B_markers$diets <- "BRD"
data1 <- as.data.frame(table(BRD_B_markers$group))
data1$diets <- "BRD"

WRD_B <- readRDS("./WRD_B.rds")
df1=as.data.frame(WRD_B@active.ident)
idents <- gsub("B cells","B",WRD_B$group_age)
names(idents) <- rownames(WRD_B@meta.data)
WRD_B@active.ident<-factor(idents)
WRD_B_markers <- FindMarkers(WRD_B, ident.1 = "M9_B", ident.2 = "M6_B", min.pct = 0.25)
WRD_B_markers$group <- 0
WRD_B_markers[WRD_B_markers$avg_log2FC > 0.25,]$group <- "9M"
WRD_B_markers[WRD_B_markers$avg_log2FC< -0.25,]$group <- "6M"
WRD_B_markers$diets <- "WRD"
data2 <- as.data.frame(table(WRD_B_markers$group))
data2$diets <- "WRD"

PRD_B <- readRDS("./PRD_B.rds")
df1=as.data.frame(PRD_B@active.ident)
idents <- gsub("B cells","B",PRD_B$group_age)
names(idents) <- rownames(PRD_B@meta.data)
PRD_B@active.ident<-factor(idents)
PRD_B_markers <- FindMarkers(PRD_B, ident.1 = "M9_B", ident.2 = "M6_B", min.pct = 0.25)
PRD_B_markers$group <- 0
PRD_B_markers[PRD_B_markers$avg_log2FC > 0.25,]$group <- "9M"
PRD_B_markers[PRD_B_markers$avg_log2FC< -0.25,]$group <- "6M"
PRD_B_markers$diets <- "PRD"
data3 <- as.data.frame(table(PRD_B_markers$group))
data3$diets <- "PRD"
data <- rbind(data1,data2)%>%rbind(.,data3)
data$celltype <- "B"
d5 <- data

data <- rbind(d1,d2)%>%rbind(.,d3)%>%rbind(.,d4)%>%rbind(.,d5)
data$diets <- factor(data$diets,levels = c("BRD","WRD","PRD"))
data$celltype <- factor(data$celltype,levels = c("STC","GC","M","DC","B"))

#"#F39EA2","#EC2124"

data
ggplot(data, aes(x =celltype,y = Freq,fill = Var1))+
  geom_bar(stat="identity",position = "stack")+
  facet_wrap(~diets)+
  theme_bw()+scale_fill_manual(values = c("#C4E0ED","#4E7DB8"))+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(color = 'black'),axis.line.y=element_line(color = 'black'))+theme(legend.position = 'none')+
  ylab("Differentially 
       expressed genes (DEGs)")+
  theme(strip.text.x = element_text(size=14))


########################## Cell interaction 1. cellphonedb #########################
library(tidyverse)
library(RColorBrewer)
library(scales)


###################### 9M ###############################
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
#View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
#View(statdf)
statdf$total_number=0
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
# 237
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,150))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+scale_x_discrete(limits=c("STC","GC","M","DC","B"))+
  scale_y_discrete(limits=c("STC","GC","M","DC","B"))
###################### BRD ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_BRD//pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")

statdf$total_number=0
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
View(statdf) # 226
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,160))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+scale_x_discrete(limits=c("STC","GC","DC","B",""))+
  scale_y_discrete(limits=c("STC","GC","M","DC","B"))

###################### WRD ###############################
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_WRD/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
#View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
#View(statdf)
statdf$total_number=0
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
View(statdf)
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,160))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+scale_x_discrete(limits=c("STC","GC","M","DC","B"))+
  scale_y_discrete(limits=c("STC","GC","M","DC","B"))
###################### PRD ###############################
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_PRD//pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
View(statdf)

colnames(statdf)=c("number")
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
View(statdf)
statdf$total_number=0
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
View(statdf) # 241
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,160))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+scale_x_discrete(limits=c("STC","GC","M","DC","B"))+
  scale_y_discrete(limits=c("STC","GC","M","DC","B"))



prop_test <- sc_utils(M9)
prop_test1 <- permutation_test(
  prop_test, cluster_identity = "abbreviation2",
  sample_1 = "M9-BRD", sample_2 = "M9-WRD",
  sample_identity = "Diets"
)

permutation_plot2(prop_test1,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("GLC","STC","B","DC","NK","GC","T","M"))



prop_test3 <- permutation_test(
  prop_test, cluster_identity = "abbreviation2",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)

permutation_plot2(prop_test3,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("GLC","STC","B","DC","NK","GC","T","M"))

rm(M9,counts,meta,rds,a)
gc()
rm(list=ls())
################################ M GC STC DC B ######################################
################################ M ###################################
M <- readRDS("./M.rds")

M@active.ident <- factor(M$group_age)

M_O <- FindMarkers(M, ident.1 = "M9_Macrophages", ident.2 = "M6_Macrophages", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
M_Y <- FindMarkers(M, ident.1 = "M9_Macrophages", ident.2 = "M6_Macrophages", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(M_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_M_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)

View(KEGG_M_O@result)
write.table(KEGG_M_O@result,"KEGG_M_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

KEGG_M_O_select <- read.delim("~/scRNA_aging/output/KEGG_M_O_select.txt", row.names=1)
KEGG_M_O@result <- na.omit(KEGG_M_O_select[KEGG_M_O_select$select==1,])
#KEGG_M_O@result$Description <- gsub(" - Mus musculus (house mouse)","",KEGG_M_O@result$Description)

dotplot(KEGG_M_O,label_format = 45,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))

mean_M <- apply(as.data.frame(M@assays$RNA@data),1,function(x){tapply(x,M$group_age_diet,mean)}) %>% t() %>% as.data.frame()
head(mean_M)
colnames(mean_M)=strsplit2(colnames(mean_M),"_")[,1]

select <-KEGG_M_O@result[KEGG_M_O@result$Description%in%c(
  "NF-kappa B signaling pathway","TNF signaling pathway","Toll-like receptor signaling pathway",
  "Apoptosis","IL-17 signaling pathway","PI3K-Akt signaling pathway","MAPK signaling pathway",
  "Antigen processing and presentation","JAK-STAT signaling pathway"),] 
select <-KEGG_M_O@result[KEGG_M_O@result$Description%in%c(
  
  "Apoptosis"),] 
a=c(strsplit2(select$geneID,"/"))%>%unique()
ncbi[ncbi$ENTREZID%in%a,]$SYMBOL
data <- mean_M[rownames(mean_M)%in%ncbi[ncbi$ENTREZID%in%a,]$SYMBOL,]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1

pheatmap::pheatmap(n,fontsize = 14,border_color = "white")


BP_MO=clusterProfiler::simplify(enrichGO(gene=rownames(M_O),keyType = "SYMBOL",
                                            OrgDb= org.Mm.eg.db, ont="BP",
                                            pAdjustMethod="BH",pvalueCutoff=0.05,
                                            qvalueCutoff=0.05))
View(BP_MO@result)

BP_MY=clusterProfiler::simplify(enrichGO(gene=rownames(M_Y),keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))
View(BP_MY@result)


ElbowPlot(M,ndims = 50)
M <- FindNeighbors(M, dims = 1:22)
M <- FindClusters(M, resolution = 0.4)
M <- RunUMAP(M, dims = 1:22)
M <- RunTSNE(M, dims = 1:22)

#DimPlot(M, reduction = "umap", pt.size = 1, label = T,group.by = "Diets")
DimPlot(M, reduction = "tsne", pt.size = 1, label = T,label.size = 8)
library(SingleR)
ImmGen.se=ImmGenData()
Mouse.se=MouseRNAseqData()
meta <- M@meta.data


M_for_SingleR <- GetAssayData(M, slot="data") ##获取标准化矩阵
M.hesc <- SingleR(test = M_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
M.hesc2 <- SingleR(test = M_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
M@meta.data$SingleR_Mouse <- M.hesc$labels
M@meta.data$SingleR_ImmGen <- M.hesc2$labels
DimPlot(object = M, reduction = "tsne",label =T,group.by = "SingleR_Mouse",raster=FALSE)

M1 <- subset(M,cells=rownames(meta[meta$SingleR_ImmGen%in%c("Macrophages"),]))
M1<- RunTSNE(M1, dims = 1:22)
DimPlot(object = M1, reduction = "tsne",label =T)
# Monocytes 2
# "Il13ra2","Nos1","Nos2" M1 "Cd86","Cd80","Il17b","Il1b","Cxcl9","Cxcl10" 457
# "Il13ra1","Arg1","Retnla" M2 "Cd206","Cd163","Arg1","Ym1","Fizz1","Mgl2"  236

DotPlot(M1,features = c("Gpr18","Fpr2","Fcgr1","Fcgr3","Fcgr2b",
                       "Cd86","Cd80","Il17b","Il1b","Cxcl9","Cxcl10",
                       "Cd206","Cd163","Retnla","Clec13a","Arg1","Ym1",
                       "Egr2 ","Pdcd1lg2","Mgl2","Tgm2","Mrc1","Il10","Chil1","Csf1r","Adgre1","Itgam"))



#### M1
DotPlot(M1,features = c("Gpr18","Fpr2","Fcgr1","Fcgr3","Fcgr2b","Csf1r",
                       "Cd86","Cd80","Il17b","Il1b","Cxcl9","Cxcl10","Il6","Il12"))
### M2
DotPlot(M1,features = c("Cd163","Retnla","Arg1","Il4","Il13",
                       "Pdcd1lg2","Mgl2","Tgm2","Mrc1","Il10","Chil1"))

M1 <- subset(M1,idents=c(0,2,3,4,5,6,7))
################################### M1/M2亚型命名 ###################################
new.cluster.ids <- c('Mo',
                     "M1M2",
                     'M2',
                     "M1",
                     'M1',
                     'M2',
                     'Mo')
names(x=new.cluster.ids)=levels(x = M1)
M1 <- RenameIdents(object =M1, new.cluster.ids)

DotPlot(M1,features = c("Ccr2","Itgam","Cd68","Fcgr1","Fcgr3","Fcgr2b","Mgl2","Mrc1","Il10"),
        cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))+
  scale_y_discrete(limits=rev(c('Mo',"M1","M1M2","M2")))


M1 <- RunTSNE(M1, dims = 1:22)
DimPlot(M1, reduction = "tsne", pt.size = 1, label = T,label.size = 8)+scale_color_futurama()+
  theme(axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),
        axis.text = element_text(size = 18))
DimPlot(M1, reduction = "tsne", pt.size = 1, label = T,label.size = 8,group.by = "seurat_clusters")+scale_color_futurama()
saveRDS(M1,"./M1.rds")

library("scProportionTest")
library(forcats, lib.loc = "/usr/local/lib/R/site-library")
library(plyr)
library(dplyr)
library(IRanges)
M1 <- readRDS("./M1.rds")
df1=as.data.frame(M1@active.ident)
M1@meta.data$celltype=df1$`M1@active.ident`
M1$group_age

prop_test <- sc_utils(M1)
prop_test1 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6_Macrophages", sample_2 = "M9_Macrophages",
  sample_identity = "group_age"
)

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
permutation_plot2(prop_test1,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c('Mo',"M1","M1M2","M2"))

prop_test <- sc_utils(M1)
prop_test2 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-WRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test2,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c('Mo',"M1","M1M2","M2"))


prop_test3 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test3,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c('Mo',"M1","M1M2","M2"))

prop_test4 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-WRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test4,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c('Mo',"M1","M1M2","M2"))

######################## Diff M ##################
M1 <- readRDS("./M1.rds")
df1=as.data.frame(M1@active.ident)
M1@meta.data$celltype=df1$`M1@active.ident`
M1$group <- paste0(substr(colnames(M1),1,2),"_",M1$celltype)

diff1_O <- FindMarkers(M1, ident.1 = "M9_Monocytes", ident.2 = "M6_Monocytes",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
diff1_Y <- FindMarkers(M1, ident.1 = "M9_Monocytes", ident.2 = "M6_Monocytes",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

diff2_O <- FindMarkers(M1, ident.1 = "M9_M1", ident.2 = "M6_M1",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
diff2_Y <- FindMarkers(M1, ident.1 = "M9_M1", ident.2 = "M6_M1",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

diff3_O <- FindMarkers(M1, ident.1 = "M9_M2", ident.2 = "M6_M2",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
diff3_Y <- FindMarkers(M1, ident.1 = "M9_M2", ident.2 = "M6_M2",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]


ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(diff1_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Mo_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)

View(KEGG_Mo_O@result)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(diff1_Y),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_Mo_Y <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)

View(KEGG_Mo_Y@result)










data1 <- as.data.frame(t(M@assays$RNA@data[rownames(M@assays$RNA@data)%in%c("Il13ra2","Nos1","Nos2"),]))
head(data1)
data1 <- apply(data1,1,mean)%>%as.data.frame()
data1$cells <- rownames(data1)
dim(M)
data1$Diets <- M$Diets
ggdata1 <- melt(data1)
head(ggdata1)
a=table(ggdata1[ggdata1$value>(mean(ggdata1$value)+sd(ggdata1$value)),]$Diets)%>%as.data.frame()
b=table(ggdata1[ggdata1$value<(mean(ggdata1$value)+sd(ggdata1$value)),]$Diets)%>%as.data.frame()
a
b

data2 <- as.data.frame(t(M@assays$RNA@data[rownames(M@assays$RNA@data)%in%c("Il13ra1","Arg1","Retnla"),]))
data2 <- apply(data2,1,mean)%>%as.data.frame()
data2$cells <- rownames(data2)
dim(M)
data2$Diets <- M$Diets
ggdata2 <- melt(data2)
head(ggdata2)
c=table(ggdata2[ggdata2$value>(mean(ggdata2$value)+sd(ggdata2$value)),]$Diets)%>%as.data.frame()
d=table(ggdata2[ggdata2$value<(mean(ggdata2$value)+sd(ggdata2$value)),]$Diets)%>%as.data.frame()
c
d

data3 <- as.data.frame(t(M@assays$RNA@data[rownames(M@assays$RNA@data)%in%c("Mafb","Fcgr1","Axl","Ly6c1"),]))
data3 <- apply(data3,1,mean)%>%as.data.frame()

data3$cells <- rownames(data3)
dim(M)
data3$Diets <- M$Diets
ggdata3 <- melt(data3)
head(ggdata3)
e=table(ggdata3[ggdata3$value>(mean(ggdata3$value)+sd(ggdata3$value)),]$Diets)%>%as.data.frame()
f=table(ggdata3[ggdata3$value<(mean(ggdata3$value)+sd(ggdata3$value)),]$Diets)%>%as.data.frame()

e
f
### GSVA #####
library("GSVA")

mmu_kegg <- clusterProfiler::download_KEGG("mmu")
head(mmu_kegg$KEGGPATHID2NAME)
KEGGPATHIDNAME <- as.data.frame(mmu_kegg$KEGGPATHID2NAME)
head(mmu_kegg$KEGGPATHID2EXTID)
PATH2ID <- mmu_kegg$KEGGPATHID2EXTID
PATH2NAME <- mmu_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
path_mmu<- AnnotationDbi::select(org.Mm.eg.db,keys=PATH_ID_NAME$ENTREZID,columns=c("ENTREZID","SYMBOL"),keytype="ENTREZID")
PATH_mmu<- merge(path_mmu,PATH_ID_NAME,all.x=TRUE,by ="ENTREZID")
mmu_hs <- read.delim("/home/tuyx/scRNA-Seq/new_10X/New_0.8/mmu_hs.txt")
head(mmu_hs)
mmu_hs_dedup=mmu_hs[!duplicated(mmu_hs$Gene.name),]

kegggmt <- read.gmt("/home/tuyx/scRNA-Seq/new_10X/New_0.8/c2.cp.kegg.v7.3.symbols.gmt")
head(kegggmt)
colnames(kegggmt)
unique(kegggmt$gene)
kegggmt$Hs.gene<- kegggmt$gene
colnames(mmu_hs_dedup) <- c("mmu.gene","Hs.gene")
tmp<- merge(kegggmt,mmu_hs_dedup,all.x=TRUE,by ="Hs.gene")
tmp <- na.omit(tmp)
kegggmt_mmu<- tmp[,c(2,4)]
kegg_list = split(kegggmt_mmu$mmu.gene, kegggmt_mmu$term)
saveRDS(kegg_list,"./kegg_list.rds")

gsva_M <- GSVA::gsva(as.matrix(M@assays$RNA@data),kegg_list,kcdf="Gaussian",parallel.sz=4,method = "gsva")

dim(gsva_M)
data <- gsva_M
M@meta.data$Diets
meta <- as.data.frame(M@meta.data[,c("Diets")])
rownames(meta) <- rownames(M@meta.data)

data=as.data.frame(t(apply(data,1,function(a){
  tapply(a,meta$`M@meta.data[, c("Diets")]`,mean)})))
rownames(data) <- gsub("_"," ",rownames(data))%>%gsub("KEGG ","",.)
View(data)
n=t(scale(t(data)))
n[n>1]=1
n[n<-1]=-1
set.seed(1234)
p=pheatmap::pheatmap(n,kmeans_k = 6,cluster_cols = F,border_color = "white",fontsize = 18)
d=as.data.frame(p$kmeans$cluster)
d$pathway <- rownames(d)
cluster5 <- d[d$`p$kmeans$cluster`==5,]
pheatmap::pheatmap(n[cluster5$pathway,],cluster_cols = F,border_color = "white",fontsize = 18)
pheatmap::pheatmap(n[c("HEMATOPOIETIC CELL LINEAGE",
                       "NUCLEOTIDE EXCISION REPAIR",
                       "BASE EXCISION REPAIR",
                       "PPAR SIGNALING PATHWAY",
                       "DNA REPLICATION",
                       "BIOSYNTHESIS OF UNSATURATED FATTY ACIDS",
                       "LYSOSOME","NICOTINATE AND NICOTINAMIDE METABOLISM",
                       "RIBOFLAVIN METABOLISM"),],cluster_rows  = T,cluster_cols = T,
                   fontsize = 18,border_color =NULL,main = "Macrophages")


?pheatmap

########################### 识别M衰老细胞 ####################################
SASP_HM <- read.delim("~/scRNA_aging/output/SASP_HM.txt")
colnames(SASP_HM)[2] <- "Gene"
SASP_genes <- read.delim("~/scRNA_aging/output/SASP_genes.txt")
SASP <- merge(SASP_HM,SASP_genes,all.x=T,by="Gene")
head(SASP)
SASP$Gene <- NULL
unique(SASP$Origin)
SASP_select <- subset(SASP,Weights>0.6)

## cell fraction
data <- as.data.frame(M@assays$RNA@data)
tmp <- data[rownames(data)%in%SASP_select$Mouse.gene.name,]%>%t()%>%as.data.frame()
#tmp <- data[rownames(data)%in%SASP_select[SASP_select$Origin%in%c("Senescence Initiators","Senescence Responses","Canonical Senesce Pathway"),]$Mouse.gene.name,]%>%t()%>%as.data.frame()
tmp <- as.data.frame(apply(tmp,1,mean))
colnames(tmp) <- "AVG"
dim(data)
tmp$group <- strsplit2(M@meta.data$group_age_diet,"_")[,1]
tmp$group <- factor(tmp$group,levels = c("M6-BRD","M6-WRD","M6-PRD","M9-BRD","M9-WRD","M9-PRD"))
tmp$cell <- rownames(tmp)

ggdata <- reshape2::melt(tmp)
head(ggdata)
a=table(ggdata[ggdata$value>median(ggdata$value),]$group)%>%as.data.frame()
b=table(ggdata[ggdata$value<median(ggdata$value),]$group)%>%as.data.frame()

#a=table(ggdata[ggdata$value>mean(ggdata$value)+2*sd(ggdata$value),]$group)%>%as.data.frame()
#b=table(ggdata[ggdata$value<mean(ggdata$value)+2*sd(ggdata$value),]$group)%>%as.data.frame()
a
b
frac <- data.frame(cell=a$Var1,ex=a$Freq,Nex=b$Freq)
frac$all <- frac$ex+frac$Nex
frac$f <- frac$ex/frac$all*100
frac$group <- "Aged"
frac2 <- frac
frac2$f <- frac$Nex/frac$all*100
frac2$group <- "Non-Aged"
### 堆叠
data <- rbind(frac2,frac)
data$group <- factor(data$group,levels = c("Non-Aged","Aged"))
data

ggplot(data, aes( x = cell,y=f,fill = group))+
  geom_col(position = 'stack', width = 0.6)+
  theme_classic()+
  scale_fill_manual(values=c("#E7AE8F","#BF5217"))+  
  labs(x="",y="Cell senescence fraction")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=14),legend.position = "right",legend.title =element_text(size=16),
        axis.line.x=element_line(),axis.line.y=element_line())





################################ GC ###################################
GC <- readRDS("./GC.rds")
GC@active.ident <- factor(GC$group_age)

GC_O <- FindMarkers(GC, ident.1 = "M9_Granulosa cells", ident.2 = "M6_Granulosa cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
GC_Y <- FindMarkers(GC, ident.1 = "M9_Granulosa cells", ident.2 = "M6_Granulosa cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC_O <- enrichKEGG(gene =ncbi$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)

View(KEGG_GC_O@result)
#KEGG_GC_O@result$Description <- gsub(" - Mus musculus (house mouse)","",KEGG_GC_O@result$Description)

write.table(KEGG_GC_O@result,"KEGG_GC_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

KEGG_GC_O_select <- read.delim("~/scRNA_aging/output/KEGG_GC_O_select.txt", row.names=1)
KEGG_GC_O@result <- na.omit(KEGG_GC_O_select[KEGG_GC_O_select$select==1,])
KEGG_GC_O <-new("enrichResult",result=KEGG_GC_O@result)
dotplot(KEGG_GC_O,label_format = 60,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))


mean_GC <- apply(as.data.frame(GC@assays$RNA@data),1,function(x){tapply(x,GC$group_age_diet,mean)}) %>% t() %>% as.data.frame()
head(mean_GC)
dim(mean_GC)
colnames(mean_GC)=strsplit2(colnames(mean_GC),"_")[,1]

select <-KEGG_GC_O@result[KEGG_GC_O@result$Description%in%c(
  "PI3K-Akt signaling pathway","Glutathione metabolism ",
  "Chemical carcinogenesis - DNA adducts",
  "mTOR signaling pathway","TNF signaling pathway",
  "IL-17 signaling pathway",
  "Chemical carcinogenesis - reactive oxygen species"),] 

a=c(strsplit2(select$geneID,"/"))%>%unique()
ncbi[ncbi$ENTREZID%in%a,]$SYMBOL
data <- mean_GC[rownames(mean_GC)%in%c(ncbi[ncbi$ENTREZID%in%a,]$SYMBOL,
                                       c("Cdkn2a","Cdkn2d","Cdkn1b")),]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1
rownames(data)
pheatmap::pheatmap(n[-c(10:11),c(1,3,2,4,6,5)],fontsize = 16,border_color = "black",cluster_cols = F)





BP_GCO=clusterProfiler::simplify(enrichGO(gene=rownames(GC_O),keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))
View(BP_GCO@result)
BP_GCY=clusterProfiler::simplify(enrichGO(gene=rownames(GC_Y),keyType = "SYMBOL",
                                         OrgDb= org.Mm.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))
View(BP_GCY@result)


geneList_GC = GC_OY$avg_log2FC
names(geneList_GC) = rownames(GC_OY)
geneList_GC = sort(geneList_GC, decreasing = TRUE)
gsea_GC <- gseGO(geneList     = geneList_GC,keyType = "SYMBOL",
                OrgDb        = org.Mm.eg.db,
                ont          = "BP",
                pvalueCutoff = 1,
                verbose      = 3)
gsea_GC_table <- as.data.frame(gsea_GC)%>%subset(.,pvalue<0.25)
View(gsea_GC_table)


?gsva
GC@meta.data$Diets
meta <- as.data.frame(GC@meta.data[,c("Diets")])
rownames(meta) <- rownames(GC@meta.data)
dim(meta)
mean = apply(as.data.frame(GC@assays$RNA@data),1,function(x){tapply(x,meta$`GC@meta.data[, c("Diets")]`,mean)}) %>% t() %>% as.data.frame()
gsva_GC <- GSVA::gsva(as.matrix(mean),kegg_list,kcdf="Gaussian",parallel.sz=2,method = "gsva")

dim(gsva_GC)
data <- gsva_GC
dim(data)
VlnPlot(GC1, features = c("Myc","Apoe","Nd3"),pt.size = 0,ncol = 3)+NoLegend()

GC1$group_age_diet
VlnPlot(GC1, features = c("Fosl2","Foxo3","Junb","Jund","Stat3","Fos"),pt.size = 0,ncol = 3,group.by = "group_age_diet")+NoLegend()

GC1 <- readRDS("./GC1.rds")
df1=as.data.frame(GC1@active.ident)
GC1@meta.data%>%head()
idents <- paste0(substr(rownames(GC1@meta.data),1,2),"-",df1$`GC1@active.ident`)
names(idents) <- rownames(GC1@meta.data)
GC1@active.ident<-factor(idents)
AN_GC1 <- subset(GC1,idents=c("M6-Antral GC","M9-Antral GC"))
AN_GC1$Group <- AN_GC1@active.ident
AN_GC1$group_age_diet <- gsub("_Granulosa cells","-AnGC",AN_GC1$group_age_diet)
VlnPlot(AN_GC1, features = c("Foxo3"),pt.size = 0,ncol = 2,group.by = "group_age_diet")+NoLegend()
VlnPlot(GC1, features = c("Fosl2","Foxo3","Junb","Jund","Stat3","Fos"),pt.size = 0,ncol = 3,group.by = "group_age_diet")+NoLegend()

AN_GC1


#data=as.data.frame(t(apply(data,1,function(a){
#  tapply(a,meta$`GC@meta.data[, c("Diets")]`,mean)})))

rownames(data) <- gsub("_"," ",rownames(data))%>%gsub("KEGG ","",.)

View(data)
n=t(scale(t(data)))
n[n>1]=1
n[n<-1]=-1
set.seed(1234)
p=pheatmap::pheatmap(n,kmeans_k = 6,cluster_cols = F,border_color = "white",fontsize = 18)
d=as.data.frame(p$kmeans$cluster)
d$pathway <- rownames(d)
cluster23 <- d[d$`p$kmeans$cluster`%in%c(2,3),]
pheatmap::pheatmap(n[c("NITROGEN METABOLISM","P53 SIGNALING PATHWAY",
                       "NOTCH SIGNALING PATHWAY","VEGF SIGNALING PATHWAY","ECM RECEPTOR INTERACTION",
                       "INSULIN SIGNALING PATHWAY"),],cluster_cols = T,border_color = NA,fontsize = 18,
                   main="Granulosa cells")





ElbowPlot(GC,ndims = 50)

GC <- FindNeighbors(GC, dims = 1:22)
GC <- FindClusters(GC, resolution = 0.4)
GC <- RunUMAP(GC, dims = 1:22)
GC <- RunTSNE(GC, dims = 1:22)
GC$group_age
DimPlot(GC, reduction = "tsne", pt.size = 1, label = T,group.by = "group_age")
DimPlot(GC, reduction = "tsne", pt.size = 1, label = T,label.size = 8)
FeaturePlot(GC,features = c("Nanos3","Kit","Nanog","Dazl","Ddx4","Prom1","Dppa3"),reduction = "tsne",pt.size = 1.5,order = TRUE,label = T)& theme(legend.position = "right")
FeaturePlot(GC,features = c("Myc","Apoe"),reduction = "tsne",pt.size = 1.5,order = TRUE,label = T)& theme(legend.position = "right")
FeaturePlot(BP,features = c("Foxo3","Jund"),reduction = "tsne",pt.size = 1.5,order = TRUE,label = T)& theme(legend.position = "right")

# 3 5 Mitotic GC
# 2 6 Preantral GC
# 4 7 Atretic GC
# 1 Antral GC
#  
DotPlot(GC,features = c("Inhbb","Fst","Gja1",
                        "Itih5","Cald1","Pik3ip1",
                        "Igfbp5","Gatm","Col18a1",
                        "Top2a","Cenpa","Pcna",
                        "Cyp11a1","Ptgfr","Onecut2",
                        "Inhba","Bx4","Nr5a2",
                        "Foxl2","Cyp19a1","Fshr"),cols = c("#584EA1","#ED1C26"))+theme(axis.text.x = element_text(angle = 45,hjust =1))

GC_markers <- FindAllMarkers(object = GC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(GC_markers[GC_markers$cluster==0,])
BP=clusterProfiler::simplify(enrichGO(gene=rownames(GC_markers[GC_markers$cluster==0,]),keyType = "SYMBOL",
                                          OrgDb= org.Mm.eg.db, ont="BP",
                                          pAdjustMethod="BH",pvalueCutoff=0.05,
                                          qvalueCutoff=0.05))
View(BP@result)
GC1 <- subset(GC,idents=c(1:7))
new.cluster.ids <- c('Antral GC',
                     'Preantral GC',
                     'Mitotic GC',
                     "Atretic GC",
                     'Mitotic GC',
                     'Preantral GC',
                     'Atretic GC')
names(x=new.cluster.ids)=levels(x = GC1)
GC1 <- RenameIdents(object =GC1, new.cluster.ids)

DotPlot(GC1,features = c("Top2a","Cenpa","Pcna",
                         "Inhbb","Fst","Gja1",
                         "Igfbp5","Gatm","Col18a1",
                         "Inhba","Bx4","Nr5a2",
                        "Foxl2","Cyp19a1","Fshr",
                        "Itih5","Cald1","Pik3ip1"),
        cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))+
  scale_y_discrete(limits=rev(c("Mitotic GC","Preantral GC","Antral GC","Atretic GC")))


GC1 <- RunTSNE(GC1, dims = 1:22)
saveRDS(GC1,"./GC1.rds")
DimPlot(GC1, reduction = "tsne", pt.size = 1, label = T,label.size = 8)+scale_color_d3()


######################### GC1 KEGG 9M ######################
GC1 <- readRDS("./GC1.rds")
GC1$group <- paste0(substr(GC1$Diets,1,6),"-",GC1@active.ident)
GC1@meta.data$group_age

GC1_OY<- FindMarkers(GC1, ident.1 = "M9_Granulosa cells", ident.2 = "M6_Granulosa cells",group.by = "group_age", min.pct = 0.25)
write.csv(GC1_OY,"GC1_OY.csv")


colnames(mean_GC)
median_GC <- apply(as.data.frame(GC1@assays$RNA@data),1,function(x){tapply(x,GC1$group,median)}) %>% t() %>% as.data.frame()
View(t(mean_GC[c("Fosl2","Foxo3","Junb","Jund","Stat3"),]))
pheatmap::pheatmap(t(scale(t(median_GC[c("Fosl2","Foxo3","Junb","Jund","Stat3","Fos"),
                                     c("M6-BRD-Antral GC","M6-WRD-Antral GC","M6-PRD-Antral GC",
                                       "M9-BRD-Antral GC","M9-WRD-Antral GC","M9-PRD-Antral GC")]))),cluster_cols = F)

mean_GC["Fos",]
head(mean_GC)

################## BRD ##########################
dim(Mito)
Mito<- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M6-BRD-Mitotic GC",group.by = "group", min.pct = 0.25)
Pre<- FindMarkers(GC1, ident.1 = "M9-BRD-Preantral GC", ident.2 = "M6-BRD-Preantral GC",group.by = "group", min.pct = 0.25)
Ant<- FindMarkers(GC1, ident.1 = "M9-BRD-Antral GC", ident.2 = "M6-BRD-Antral GC",group.by = "group", min.pct = 0.25)
Atr<- FindMarkers(GC1, ident.1 = "M9-BRD-Atretic GC", ident.2 = "M6-BRD-Atretic GC",group.by = "group", min.pct = 0.25)
data1=data.frame(cell=c("Mit GC","Pre GC","Ant GC","Atr GC"),diff=c(dim(Mito)[1],
                                                                    dim(Pre)[1],
                                                                    dim(Ant)[1],
                                                                    dim(Atr)[1]))
data1$diets <- "BRD"
################## WRD ##########################
Mito2<- FindMarkers(GC1, ident.1 = "M9-WRD-Mitotic GC", ident.2 = "M6-WRD-Mitotic GC",group.by = "group", min.pct = 0.25)
Pre2<- FindMarkers(GC1, ident.1 = "M9-WRD-Preantral GC", ident.2 = "M6-WRD-Preantral GC",group.by = "group", min.pct = 0.25)
Ant2<- FindMarkers(GC1, ident.1 = "M9-WRD-Antral GC", ident.2 = "M6-WRD-Antral GC",group.by = "group", min.pct = 0.25)
Atr2<- FindMarkers(GC1, ident.1 = "M9-WRD-Atretic GC", ident.2 = "M6-WRD-Atretic GC",group.by = "group", min.pct = 0.25)
data2=data.frame(cell=c("Mit GC","Pre GC","Ant GC","Atr GC"),diff=c(dim(Mito2)[1],
                                                                    dim(Pre2)[1],
                                                                    dim(Ant2)[1],
                                                                    dim(Atr2)[1]))
data2$diets <- "WRD"
################## PRD ##########################
Mito3<- FindMarkers(GC1, ident.1 = "M9-PRD-Mitotic GC", ident.2 = "M6-PRD-Mitotic GC",group.by = "group", min.pct = 0.25)
Pre3<- FindMarkers(GC1, ident.1 = "M9-PRD-Preantral GC", ident.2 = "M6-PRD-Preantral GC",group.by = "group", min.pct = 0.25)
Ant3<- FindMarkers(GC1, ident.1 = "M9-PRD-Antral GC", ident.2 = "M6-PRD-Antral GC",group.by = "group", min.pct = 0.25)
Atr3<- FindMarkers(GC1, ident.1 = "M9-PRD-Atretic GC", ident.2 = "M6-PRD-Atretic GC",group.by = "group", min.pct = 0.25)

data3=data.frame(cell=c("Mit GC","Pre GC","Ant GC","Atr GC"),diff=c(dim(Mito3)[1],
                                                                    dim(Pre3)[1],
                                                                    dim(Ant3)[1],
                                                                    dim(Atr3)[1]))
data3$diets <- "PRD"

data <- rbind(data1,data2)%>%rbind(.,data3)
data$diets <- factor(data$diets,levels = c("BRD","WRD","PRD"))
data$cell <- factor(data$cell,levels = c("Mit GC", "Pre GC", "Ant GC", "Atr GC"))
ggplot(data, aes(x =cell,y = diff))+
  geom_bar(stat="identity")+
  facet_wrap(~diets)+
  theme_bw()+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle =45,hjust=1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(color = 'black'),axis.line.y=element_line(color = 'black'))+theme(legend.position = 'none')+
  ylab("Differentially 
       expressed genes (DEGs)")+
  theme(strip.text.x = element_text(size=14))

aa############################### GC KEGG Select ###########################
KEGG_GC1_1_O <- read.delim("~/scRNA_aging/output/KEGG_GC1_1_O.txt", row.names=1)%>%.[.$select==1,]%>%na.omit()
KEGG_GC1_2_O <- read.delim("~/scRNA_aging/output/KEGG_GC1_2_O.txt", row.names=1)%>%.[.$select==1,]%>%na.omit()
KEGG_GC1_3_O <- read.delim("~/scRNA_aging/output/KEGG_GC1_3_O.txt", row.names=1)%>%.[.$select==1,]%>%na.omit()
KEGG_GC1_4_O <- read.delim("~/scRNA_aging/output/KEGG_GC1_4_O.txt", row.names=1)%>%.[.$select==1,]%>%na.omit()
KEGG_GC1_1_O$cell <- "Mit GC"
KEGG_GC1_2_O$cell <- "Pre GC"
KEGG_GC1_3_O$cell <- "Ant GC"
KEGG_GC1_4_O$cell <- "Atr GC"
KEGG_GC1_select <- rbind(KEGG_GC1_1_O,KEGG_GC1_2_O)%>%rbind(.,KEGG_GC1_3_O)%>%rbind(.,KEGG_GC1_4_O)
b <- lapply(str_split(KEGG_GC1_select$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
KEGG_GC1_select$GeneRatio <- b
ggplot(KEGG_GC1_select,aes(cell,Description,size = GeneRatio))+
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
  ) +scale_x_discrete(limits=c("Mit GC","Pre GC","Ant GC","Atr GC"))+
  scale_y_discrete(limits= unique(KEGG_GC1_select$Description)) 

a=KEGG_GC1_select[KEGG_GC1_select$Description%in%c("Reactive oxygen species",
                                                   "Cytokine-cytokine receptor interaction",
                                                   "Insulin signaling pathway",
                                                   "ECM-receptor interaction",
                                                   "Receptor activation",
                                                   "mTOR signaling pathway"),]$geneID%>%str_split("/")%>%unlist()

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=unique(a),columns=c("SYMBOL","ENTREZID"),keytype="ENTREZID")
data <- mean_GC[rownames(mean_GC)%in%ncbi$SYMBOL,]
write.csv(data,"aging_gene_GC.csv")
pheatmap::pheatmap(data,cluster_cols = F,cluster_rows = T,fontsize = 14,border_color = "white")

n <- t(scale(t(data)))
n[n>1]=2
n[n<-1]=-1
pheatmap::pheatmap(n,cluster_cols = T,cluster_rows = T,fontsize = 14,border_color = "black")

################ GSVA ################
gsva_GC <- GSVA::gsva(as.matrix(mean_GC),kegg_list,kcdf="Gaussian",parallel.sz=2,method = "gsva")
gsva_GC




GC1_BW<- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M9-WRD-Mitotic GC",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
GC1_WB<- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M9-WRD-Mitotic GC",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

library(org.Mm.eg.db)
library(clusterProfiler)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_BW),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_BW_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
View(KEGG_GC1_BW_O@result)
KEGG_GC1_BW_O@result <- subset(KEGG_GC1_BW_O@result,pvalue<0.05)
GC1_BP <- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M9-PRD-Mitotic GC",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
GC1_PB <- FindMarkers(GC1, ident.1 = "M9-BRD-Mitotic GC", ident.2 = "M9-PRD-Mitotic GC",group.by = "group", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]


ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_BP),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_BP_O <- enrichKEGG(gene =ncbi$ENTREZID,
                            organism = 'mmu',
                            pvalueCutoff =0.05)
View(KEGG_GC1_BP_O@result)
KEGG_GC1_BP_O@result <- subset(KEGG_GC1_BP_O@result,pvalue<0.05)

intersect(KEGG_GC1_BW_O@result$Description,KEGG_GC1_BP_O@result$Description)

intersect(rownames(GC1_BW),rownames(GC1_BP))
intersect(rownames(GC1_WB),rownames(GC1_PB))
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=intersect(rownames(GC1_WB),rownames(GC1_PB)),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC_O <- enrichKEGG(gene =ncbi$ENTREZID,
                            organism = 'mmu',
                            pvalueCutoff =0.05)
View(KEGG_GC_O@result)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_WB),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_WB_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_PB),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_PB_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)
KEGG_WB_O@result <- subset(KEGG_WB_O@result,pvalue<0.05)
KEGG_PB_O@result <- subset(KEGG_PB_O@result,pvalue<0.05)
intersect(KEGG_WB_O@result$Description,KEGG_PB_O@result$Description)

################################ GC1 count ###################################
table(GC1$orig.ident)
counts <- as.data.frame(table(GC1@active.ident,GC1@meta.data$orig.ident))
counts
length(colnames(counts))
# 145144

a <- table(GC1@meta.data$orig.ident)%>%as.data.frame()
a
counts$all <- c(rep(a$Freq[1],4),rep(a$Freq[2],4),rep(a$Freq[3],4),
                rep(a$Freq[4],4),rep(a$Freq[5],4),rep(a$Freq[6],4),
                rep(a$Freq[7],4),rep(a$Freq[8],4),rep(a$Freq[9],4),
                rep(a$Freq[10],4),rep(a$Freq[11],4),rep(a$Freq[12],4),
                rep(a$Freq[13],4),rep(a$Freq[14],4),rep(a$Freq[15],4),
                rep(a$Freq[16],4),rep(a$Freq[17],4))
counts$Frac <- counts$Freq/counts$all
head(counts)
counts$Diets <- paste0(strsplit2(counts$Var2,"-")[,2],"RD")
counts$Diets <- factor(counts$Diets,levels = c("BRD","WRD","PRD"))
counts$Stage <- strsplit2(counts$Var2,"-")[,1]
counts$Celltype <- factor(counts$Var1,levels=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))

data <- counts[counts$Stage=="M9",]
ggplot(data,aes(x=Diets,y=Frac*100,fill=Diets))+#stat_boxplot(geom = "errorbar",width=0.6)+
  geom_boxplot(width=0.6)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  xlab("Diets")+ylab("Cell Component Fractions %")+
# geom_signif(comparisons = list(c("BRD","WRD"),c("BRD","PRD"),c("WRD","PRD")),
#              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2,na.rm = T)+
  facet_wrap( ~Celltype, ncol=2,scales = "free_y")+scale_fill_manual(values=c("#E33E30","#F9C42C","#307DEE"))+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())



library("scProportionTest")
df1=as.data.frame(GC1@active.ident)
GC1@meta.data$celltype=df1$`GC1@active.ident`
GC1$group_age
prop_test <- sc_utils(GC1)
prop_test1 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6_Granulosa cells", sample_2 = "M9_Granulosa cells",
  sample_identity = "group_age"
)

permutation_plot2(prop_test1,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))

prop_test <- sc_utils(GC1)
prop_test2 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-WRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test2,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))


prop_test3 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test3,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("Mitotic GC","Preantral GC","Antral GC","Atretic GC"))




################################ STC ###################################
STC <- readRDS("./STC.rds")
STC@active.ident <- factor(STC$group_age)

#"Smc4","Mki67","Smc2","Mns1","Top2a","Nppc","Fbxo5"



STC_O <- FindMarkers(STC, ident.1 = "M9_Stem cells", ident.2 = "M6_Stem cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
STC_Y <- FindMarkers(STC, ident.1 = "M9_Stem cells", ident.2 = "M6_Stem cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(STC_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_STC_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)

View(KEGG_STC_O@result)

#KEGG_STC_O@result$Description <- gsub(" - Mus musculus (house mouse)","",KEGG_STC_O@result$Description)

write.table(KEGG_STC_O@result,"KEGG_STC_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

KEGG_STC_O_select <- read.delim("~/scRNA_aging/output/KEGG_STC_O_select.txt", row.names=1)
KEGG_STC_O@result <- na.omit(KEGG_STC_O_select[KEGG_STC_O_select$select==1,])
KEGG_STC_O <-new("enrichResult",result=KEGG_STC_O@result)

KEGG_STC_O@result
dotplot(KEGG_STC_O,label_format = 35,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))

mean_STC <- apply(as.data.frame(STC@assays$RNA@data),1,function(x){tapply(x,STC$group_age_diet,mean)}) %>% t() %>% as.data.frame()
head(mean_STC)
colnames(mean_STC)=strsplit2(colnames(mean_STC),"_")[,1]

select <-KEGG_STC_O@result[KEGG_STC_O@result$Description%in%c(
                                                              "ECM-receptor interaction",
                                                              "PI3K-Akt signaling pathway"),] 

a=c(strsplit2(select$geneID,"/"))%>%unique()
ncbi[ncbi$ENTREZID%in%a,]$SYMBOL
data <- mean_STC[rownames(mean_STC)%in%ncbi[ncbi$ENTREZID%in%a,]$SYMBOL,]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1

pheatmap::pheatmap(n,fontsize = 16,border_color = "white")


BP_STCO=clusterProfiler::simplify(enrichGO(gene=rownames(STC_O),keyType = "SYMBOL",
                                          OrgDb= org.Mm.eg.db, ont="BP",
                                          pAdjustMethod="BH",pvalueCutoff=0.05,
                                          qvalueCutoff=0.05))
View(BP_STCO@result)

BP_STCY=clusterProfiler::simplify(enrichGO(gene=rownames(STC_Y),keyType = "SYMBOL",
                                          OrgDb= org.Mm.eg.db, ont="BP",
                                          pAdjustMethod="BH",pvalueCutoff=0.05,
                                          qvalueCutoff=0.05))
View(BP_STCY@result)
?gsva
kegg_list
gsva_STC <- GSVA::gsva(as.matrix(STC@assays$RNA@data),kegg_list,kcdf="Gaussian",parallel.sz=4,method = "gsva")

dim(gsva_STC)
data <- gsva_STC
STC@meta.data$Diets
meta <- as.data.frame(STC@meta.data[,c("Diets")])
rownames(meta) <- rownames(STC@meta.data)

data=as.data.frame(t(apply(data,1,function(a){
  tapply(a,meta$`STC@meta.data[, c("Diets")]`,mean)})))
rownames(data) <- gsub("_"," ",rownames(data))%>%gsub("KEGG ","",.)
View(data)
n=t(scale(t(data)))
n[n>1]=1
n[n<-1]=-1
set.seed(1234)
p=pheatmap::pheatmap(n,kmeans_k = 6,cluster_cols = F,border_color = "white",fontsize = 18)
d=as.data.frame(p$kmeans$cluster)
d$pathway <- rownames(d)
cluster4 <- d[d$`p$kmeans$cluster`%in%c(4),]
n
pheatmap::pheatmap(n[c("ECM RECEPTOR INTERACTION","CELL ADHESION MOLECULES CAMS",
                       "SULFUR METABOLISM","GLYCAN BIOSYNTHESIS"),],
                   cluster_cols = F,cluster_rows = F,
                   border_color = "white",fontsize = 18,main = "Dendritic cells")


ElbowPlot(STC,ndims = 50)
STC <- FindNeighbors(STC, dims = 1:22)
STC <- FindClusters(STC, resolution = 0.4)
STC <- RunUMAP(STC, dims = 1:22)
STC <- RunTSNE(STC, dims = 1:22)
STC$anno_celltype
DimPlot(STC, reduction = "tsne", pt.size = 1, label = T,group.by = "Diets")
DimPlot(STC, reduction = "tsne", pt.size = 1, label = T,label.size = 8)
DimPlot(STC, reduction = "umap", pt.size = 1, label = T,label.size = 8)

STC_for_SingleR <- GetAssayData(STC, slot="data") ##获取标准化矩阵
STC.hesc <- SingleR(test = STC_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
STC.hesc2 <- SingleR(test = STC_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
STC@meta.data$SingleR_Mouse <- STC.hesc$labels
STC@meta.data$SingleR_ImmGen <- STC.hesc2$labels
DimPlot(object = STC, reduction = "tsne",label =T,group.by = "SingleR_Mouse",raster=FALSE)
DimPlot(object = STC, reduction = "tsne",label =T,group.by = "SingleR_ImmGen",raster=FALSE)


FeaturePlot(STC,features = c("Nanos3","Kit","Nanog","Dazl","Ddx4","Prom1","Dppa3","Ifitm3"),reduction = "tsne",pt.size = 1.5,order = TRUE,label = T)& theme(legend.position = "right")
FeaturePlot(STC,features = c("Ddx4","Ifitm3","Dppa3","Pou5f1","Nanog","Dazl","Prom1"),reduction = "tsne",pt.size = 1.5,order = TRUE,label = T)& theme(legend.position = "right")
DotPlot(STC,features = c("Ddx4","Dppa3","Pou5f1","Nanog","Dazl","Ifitm3"))


# 
DotPlot(STC,features = c("Top2a","Cenpa","Racgap1","Pcna","Itih5","Cald1"))+theme(axis.text.x = element_text(angle = 45,hjust =1))

new.cluster.ids <- c("STC","STC","STC","GC","GC","GC","GC")
names(x=new.cluster.ids)=levels(x = STC)
STC <- RenameIdents(object =STC, new.cluster.ids)
DotPlot(STC,features = c("Top2a","Cenpa","Racgap1","Pcna","Itih5","Cald1"),
        cols = c("#584EA1","#ED1C26"))+
  theme(axis.text.x = element_text(angle = 45,hjust =1),
        axis.text = element_text(size = 16))+
  scale_y_discrete(limits=rev(c("STC","GC")))


STC <- RunTSNE(STC, dims = 1:30)
DimPlot(STC, reduction = "tsne", pt.size = 1, label = T,label.size = 8)+
  scale_color_jco()+
  theme(axis.text.x = element_text(size = 18),axis.text.y = element_text(size = 18),
        axis.text = element_text(size = 18))

df1=as.data.frame(STC@active.ident)
STC@meta.data$celltype=df1$`STC@active.ident`
STC$group_age
prop_test <- sc_utils(STC)
prop_test1 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6_Stem cells", sample_2 = "M9_Stem cells",
  sample_identity = "group_age"
)

permutation_plot2(prop_test1,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("STC","GC"))

prop_test <- sc_utils(STC)
prop_test2 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-WRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test2,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("STC","GC"))


prop_test3 <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M9-BRD", sample_2 = "M9-PRD",
  sample_identity = "Diets"
)
permutation_plot2(prop_test3,log2FD_threshold = log2(1),FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  scale_x_discrete(limits=c("STC","GC"))


########################## GC 衰老轨迹构建 #################################
GC1 <- readRDS("./GC1.rds")
GC1$stage <- substr(GC1$Diets,1,2)
df1=as.data.frame(GC1@active.ident)
GC1@meta.data$celltype=df1$`GC1@active.ident`
GC1@meta.data$stage_cell <- paste0(GC1$stage,"-",GC1@meta.data$celltype)

df1=as.data.frame(GC1@active.ident)
idents <- paste0(GC1$stage,"-",df1$`GC1@active.ident`)
names(idents) <- rownames(GC1@meta.data)
GC1@active.ident<-factor(idents)

GC1_1 <- FindMarkers(GC1, ident.1 = "M9-Mitotic GC", ident.2 = "M6-Mitotic GC", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
GC1_2 <- FindMarkers(GC1, ident.1 = "M9-Preantral GC", ident.2 = "M6-Preantral GC", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
GC1_3 <- FindMarkers(GC1, ident.1 = "M9-Antral GC", ident.2 = "M6-Antral GC", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
GC1_4 <- FindMarkers(GC1, ident.1 = "M9-Atretic GC", ident.2 = "M6-Atretic GC", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
mean_GC1<-AverageExpression(GC1)
head(mean_GC1)
data<- mean_GC1$RNA[c(rownames(GC1_1),rownames(GC1_2),rownames(GC1_3),rownames(GC1_4)),]
View(GC1_3)
data<- mean_GC1$RNA[rownames(GC1_3),]

head(data)
n=t(scale(t(data)))
dim(n)
pheatmap::pheatmap(na.omit(n)[,c(1,5)],show_rownames = F,fontsize = 18,cluster_rows = T,cluster_cols = F)

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_1[GC1_1$avg_log2FC>0,]),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_1_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)


View(KEGG_GC1_1_O@result)
write.table(KEGG_GC1_1_O@result,"KEGG_GC1_1_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_2[GC1_2$avg_log2FC>0,]),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_2_O <- enrichKEGG(gene =ncbi$ENTREZID,
                           organism = 'mmu',
                           pvalueCutoff =0.05)


View(KEGG_GC1_2_O@result)
write.table(KEGG_GC1_2_O@result,"KEGG_GC1_2_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_3[GC1_3$avg_log2FC>0,]),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_3_O <- enrichKEGG(gene =ncbi$ENTREZID,
                           organism = 'mmu',
                           pvalueCutoff =0.05)

BP_GC1_3_O=clusterProfiler::simplify(enrichGO(gene=rownames(GC1_3[GC1_3$avg_log2FC>0,]),keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_GC1_3_O@result)
View(KEGG_GC1_3_O@result)
write.table(KEGG_GC1_3_O@result,"KEGG_GC1_3_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_4[GC1_4$avg_log2FC>0,]),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_GC1_4_O <- enrichKEGG(gene =ncbi$ENTREZID,
                           organism = 'mmu',
                           pvalueCutoff =0.05)


View(KEGG_GC1_4_O@result)
write.table(KEGG_GC1_4_O@result,"KEGG_GC1_4_O.txt",col.names=T,row.names=T,quote=F,sep="\t")


KEGG_GC1_1_O_select <- read.delim("~/scRNA_aging/output/KEGG_GC1_1_O.txt", row.names=1)
KEGG_GC1_1_O@result <- na.omit(KEGG_GC1_1_O_select[KEGG_GC1_1_O_select$select==1,])
KEGG_GC1_1_O <-new("enrichResult",result=KEGG_GC1_1_O@result)
KEGG_GC1_1_O@result
dotplot(KEGG_GC1_1_O,label_format = 50,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))

KEGG_GC1_2_O_select <- read.delim("~/scRNA_aging/output/KEGG_GC1_2_O.txt", row.names=1)
KEGG_GC1_2_O@result <- na.omit(KEGG_GC1_2_O_select[KEGG_GC1_2_O_select$select==1,])
KEGG_GC1_2_O <-new("enrichResult",result=KEGG_GC1_2_O@result)
KEGG_GC1_2_O@result
dotplot(KEGG_GC1_2_O,label_format = 50,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))

KEGG_GC1_3_O_select <- read.delim("~/scRNA_aging/output/KEGG_GC1_3_O.txt", row.names=1)
KEGG_GC1_3_O@result <- na.omit(KEGG_GC1_3_O_select[KEGG_GC1_3_O_select$select==1,])
KEGG_GC1_3_O <-new("enrichResult",result=KEGG_GC1_3_O@result)
KEGG_GC1_3_O@result
dotplot(KEGG_GC1_3_O,label_format = 30,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))


library(limma)
a=c(strsplit2(KEGG_GC1_3_O@result$geneID,"/")[1:3,])%>%unique()
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(GC1_3[GC1_3$avg_log2FC>0,]),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
ncbi[ncbi$ENTREZID%in%a,]$SYMBOL
a
mean <-  apply(as.data.frame(GC1@assays$RNA@data),1,function(x){
  tapply(x,paste0(substr(GC1$group_age_diet,1,7),GC1$celltype),mean)}) %>% t() %>% as.data.frame()
data <- mean[rownames(mean)%in%ncbi[ncbi$ENTREZID%in%a,]$SYMBOL,]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1
colnames(data)
pheatmap::pheatmap(na.omit(n)[,c(1,5,9,13,17,21)],fontsize = 16,border_color = "white",cluster_cols = F)






KEGG_GC1_4_O_select <- read.delim("~/scRNA_aging/output/KEGG_GC1_4_O.txt", row.names=1)
KEGG_GC1_4_O@result <- na.omit(KEGG_GC1_4_O_select[KEGG_GC1_4_O_select$select==1,])
KEGG_GC1_4_O <-new("enrichResult",result=KEGG_GC1_4_O@result)
KEGG_GC1_4_O@result
dotplot(KEGG_GC1_4_O,label_format = 50,color = "pvalue")+theme_classic()+
  theme(axis.title  = element_text(size=18),axis.text = element_text(size=18),
        legend.text = element_text(size=14),legend.title  = element_text(size=14))




########################## GC 衰老轨迹构建 #################################
GC1@active.ident


BRD <- subset(GC1,cells=rownames(GC1@meta.data[GC1@meta.data$Diets%in%c("M6-BRD","M9-BRD"),]))
WRD <- subset(GC1,cells=rownames(GC1@meta.data[GC1@meta.data$Diets%in%c("M6-WRD","M9-WRD"),]))
PRD <- subset(GC1,cells=rownames(GC1@meta.data[GC1@meta.data$Diets%in%c("M6-PRD","M9-PRD"),]))
saveRDS(BRD,"./B_GC.rds")
saveRDS(WRD,"./W_GC.rds")
saveRDS(PRD,"./P_GC.rds")

gc()
GC1$celltype

devtools::load_all("~/monocle")
library(monocle)
GC1 <- readRDS("./GC1.rds")
AN_GC1 <- subset(GC1,idents=c("M6-Antral GC","M9-Antral GC"))
meta1=AN_GC1@meta.data
Mit_GC1 <- subset(GC1,idents=c("M6-Mitotic GC","M9-Mitotic GC"))
meta2 <- Mit_GC1@meta.data
rownames(meta1[grep("BRD",meta1$Diets),])
AnB_GC1 <- subset(AN_GC1,cells=rownames(meta1[grep("BRD",meta1$Diets),]))
AnW_GC1 <- subset(AN_GC1,cells=rownames(meta1[grep("WRD",meta1$Diets),]))
AnP_GC1 <- subset(AN_GC1,cells=rownames(meta1[grep("PRD",meta1$Diets),]))

data <- as(as.matrix(AnB_GC1@assays$RNA@data), 'sparseMatrix')
pd <- AnB_GC1@meta.data
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
sc_cds$group <- paste0(sc_cds$Diets,"-",sc_cds$celltype)
?plot_cell_trajectory
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

saveRDS(sc_cds,"AnB_GC1_sc_cds.rds")
diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))+theme(legend.position = "right")
saveRDS(sc_cds,"AnB_GC2_sc_cds.rds")



data <- as(as.matrix(AnW_GC1@assays$RNA@data), 'sparseMatrix')
pd <- AnW_GC1@meta.data
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
?plot_cell_trajectory
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

saveRDS(sc_cds,"AnW_GC1_sc_cds.rds")
diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~stage_cell")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))+theme(legend.position = "right")
saveRDS(sc_cds,"AnW_GC2_sc_cds.rds")



data <- as(as.matrix(AnP_GC1@assays$RNA@data), 'sparseMatrix')
pd <- AnP_GC1@meta.data
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
?plot_cell_trajectory
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

saveRDS(sc_cds,"AnP_GC1_sc_cds.rds")
diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~stage_cell")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds
plot_cell_trajectory(sc_cds, color_by ="stage_cell",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))+theme(legend.position = "right")
saveRDS(sc_cds,"AnP_GC2_sc_cds.rds")






################################ DC ###################################
DC <- readRDS("./DC.rds")
DC@active.ident <- factor(DC$group_age)

DC_O <- FindMarkers(DC, ident.1 = "M9_Dendritic cells", ident.2 = "M6_Dendritic cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
DC_Y <- FindMarkers(DC, ident.1 = "M9_Dendritic cells", ident.2 = "M6_Dendritic cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(DC_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_DC_O <- enrichKEGG(gene =ncbi$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)

View(KEGG_DC_O@result)
#KEGG_DC_O@result$Description <- gsub(" - Mus musculus (house mouse)","",KEGG_DC_O@result$Description)

write.table(KEGG_DC_O@result,"KEGG_DC_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

KEGG_DC_O_select <- read.delim("~/scRNA_aging/output/KEGG_DC_O_select.txt", row.names=1)
KEGG_DC_O@result <- na.omit(KEGG_DC_O_select[KEGG_DC_O_select$select==1,])
KEGG_DC_O <-new("enrichResult",result=KEGG_DC_O@result)




BP_DCO=clusterProfiler::simplify(enrichGO(gene=rownames(DC_O),keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_DCO@result)

BP_DCY=clusterProfiler::simplify(enrichGO(gene=rownames(DC_Y),keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_DCY@result)


gsva_DC <- GSVA::gsva(as.matrix(DC@assays$RNA@data),kegg_list,kcdf="Gaussian",parallel.sz=4,method = "gsva")

dim(gsva_DC)
data <- gsva_DC
DC@meta.data$Diets
meta <- as.data.frame(DC@meta.data[,c("Diets")])
rownames(meta) <- rownames(DC@meta.data)

data=as.data.frame(t(apply(data,1,function(a){
  tapply(a,meta$`DC@meta.data[, c("Diets")]`,mean)})))
rownames(data) <- gsub("_"," ",rownames(data))%>%gsub("KEGG ","",.)
View(data)
n=t(scale(t(data)))
n[n>1]=1
n[n<-1]=-1
set.seed(1234)
p=pheatmap::pheatmap(n,kmeans_k = 6,cluster_cols = F,border_color = "white",fontsize = 18)
d=as.data.frame(p$kmeans$cluster)
d$pathway <- rownames(d)
cluster25 <- d[d$`p$kmeans$cluster`%in%c(2,5),]
cluster1 <- d[d$`p$kmeans$cluster`%in%c(1),]

pheatmap::pheatmap(n[c("ARACHIDONIC ACID METABOLISM","ASCORBATE AND ALDARATE METABOLISM","STEROID HORMONE BIOSYNTHESIS",
                       "RETINOL METABOLISM","MISMATCH REPAIR","COMPLEMENT AND COAGULATION CASCADES",
                       "ANTIGEN PROCESSING AND PRESENTATION","RNA DEGRADATION","CHEMOKINE SIGNALING PATHWAY",
                       "MAPK SIGNALING PATHWAY","MTOR SIGNALING PATHWAY","NOTCH SIGNALING PATHWAY",
                       "T CELL RECEPTOR SIGNALING PATHWAY","B CELL RECEPTOR SIGNALING PATHWAY"),],
                   cluster_cols = F,cluster_rows = F,
                   border_color = "white",fontsize = 18,main = "Dendritic cells")


################################ B ###################################
B <- readRDS("./B.rds")
B@active.ident <- factor(B$group_age)
B_O <- FindMarkers(B, ident.1 = "M9_B cells", ident.2 = "M6_B cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC>0,]
B_Y <- FindMarkers(B, ident.1 = "M9_B cells", ident.2 = "M6_B cells", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]%>%.[.$avg_log2FC<0,]

ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(B_O),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_B_O <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'mmu',
                        pvalueCutoff =0.05)

View(KEGG_B_O@result)
#KEGG_B_O@result$Description <- gsub(" - Mus musculus (house mouse)","",KEGG_B_O@result$Description)

write.table(KEGG_B_O@result,"KEGG_B_O.txt",col.names=T,row.names=T,quote=F,sep="\t")

BP_BO=clusterProfiler::simplify(enrichGO(gene=rownames(B_O),keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_BO@result)

BP_BY=clusterProfiler::simplify(enrichGO(gene=rownames(B_Y),keyType = "SYMBOL",
                                           OrgDb= org.Mm.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
View(BP_BY@result)
gsva_B <- GSVA::gsva(as.matrix(B@assays$RNA@data),kegg_list,kcdf="Gaussian",parallel.sz=4,method = "gsva")

dim(gsva_B)
data <- gsva_B
B@meta.data$Diets
meta <- as.data.frame(B@meta.data[,c("Diets")])
rownames(meta) <- rownames(B@meta.data)

data=as.data.frame(t(apply(data,1,function(a){
  tapply(a,meta$`B@meta.data[, c("Diets")]`,mean)})))
rownames(data) <- gsub("_"," ",rownames(data))%>%gsub("KEGG ","",.)
View(data)
n=t(scale(t(data)))
n[n>1]=1
n[n<-1]=-1
set.seed(123456)
p=pheatmap::pheatmap(n,kmeans_k = 6,cluster_cols = F,border_color = "white",fontsize = 18)
d=as.data.frame(p$kmeans$cluster)
d$pathway <- rownames(d)
cluster5 <- d[d$`p$kmeans$cluster`==5,]
pheatmap::pheatmap(n[cluster5$pathway,],cluster_cols = F,border_color = "white",fontsize = 18)

########################## Pathway ##############################
KEGG_GC_O@result <- subset(KEGG_GC_O@result,pvalue<0.05)
KEGG_GC_O@result$'-log10(pvalue)' <- -log10(KEGG_GC_O@result$pvalue)
KEGG_GC_O@result$class <- "GC"

KEGG_STC_O@result <- subset(KEGG_STC_O@result,pvalue<0.05)
KEGG_STC_O@result$'-log10(pvalue)' <- -log10(KEGG_STC_O@result$pvalue)
KEGG_STC_O@result$class <- "STC"

KEGG_M_O@result <- subset(KEGG_M_O@result,pvalue<0.05)
KEGG_M_O@result$'-log10(pvalue)' <- -log10(KEGG_M_O@result$pvalue)
KEGG_M_O@result$class <- "M"

KEGG_DC_O@result <- subset(KEGG_DC_O@result,pvalue<0.05)
KEGG_DC_O@result$'-log10(pvalue)' <- -log10(KEGG_DC_O@result$pvalue)
KEGG_DC_O@result$class <- "DC"

KEGG_B_O@result <- subset(KEGG_B_O@result,pvalue<0.05)
KEGG_B_O@result$'-log10(pvalue)' <- -log10(KEGG_B_O@result$pvalue)
KEGG_B_O@result$class <- "B"

write.table(KEGG_GC_O@result,"KEGG_GC_O0.05.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(KEGG_STC_O@result,"KEGG_STC_O0.05.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(KEGG_M_O@result,"KEGG_M_O0.05.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(KEGG_DC_O@result,"KEGG_DC_O0.05.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(KEGG_B_O@result,"KEGG_B_O0.05.txt",col.names=T,row.names=T,quote=F,sep="\t")

KEGG_select <- read.delim("~/scRNA_aging/output/KEGG_select.txt")
View(KEGG_select)
KEGG_select$geneID

ggplot(KEGG_select,aes(x=class,y=Description))+
  geom_tile(aes(fill=log10.pvalue...1),color="black")+
  theme(panel.background = element_blank())+scale_x_discrete(limits=c("GC","STC","M","DC","B"))+xlab("")+
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 14,color = 'black'))+	
  theme(axis.text = element_text(size = 14),axis.text.x = element_text(),	
        legend.text = element_text(size=13),legend.position = "right",legend.title = element_text(size=13),
        axis.line.x=element_line(),axis.line.y=element_line())+
  scale_fill_gradient2(low='white',high = '#0A3B7B')+labs(fill="-log10(pvalue)")
numeric(KEGG_select$GeneRatio)

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
  ) +scale_x_discrete(limits=c("GC","STC","M","DC","B"))+
  scale_y_discrete(limits= unique(KEGG_select$Description)) 


select <-KEGG_select[KEGG_select$Description%in%c(
  "mTOR signaling pathway"),] 
select
a=c(strsplit2(select$geneID,"/"))%>%unique()
a
#M9 <- readRDS("./M9_obj.rds")
ncbi<- AnnotationDbi::select(org.Mm.eg.db,keys=rownames(M9),columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
averageExp <- AverageExpression(M9,group.by = "Diets")
mean=averageExp$RNA

ncbi[ncbi$ENTREZID%in%a,]$SYMBOL

data <- mean[rownames(mean)%in%ncbi[ncbi$ENTREZID%in%a,]$SYMBOL,]
data
n <- t(scale(t(data)))
#n[n>1]=1
#n[n<1]=-1
pheatmap::pheatmap(n,fontsize = 16,border_color = "white")

################################ SC ###################################
library(org.Mm.eg.db)
mmu_kegg <- clusterProfiler::download_KEGG("mmu")
names(mmu_kegg)
head(mmu_kegg$KEGGPATHID2NAME)
KEGGPATHIDNAME <- as.data.frame(mmu_kegg$KEGGPATHID2NAME)
head(mmu_kegg$KEGGPATHID2EXTID)
PATH2ID <- mmu_kegg$KEGGPATHID2EXTID
PATH2NAME <- mmu_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")
head(PATH_ID_NAME)
path_mmu<- AnnotationDbi::select(org.Mm.eg.db,keys=PATH_ID_NAME$ENTREZID,columns=c("ENTREZID","SYMBOL"),keytype="ENTREZID")
PATH_mmu <-  merge(path_mmu,PATH_ID_NAME,all.x=TRUE,by ="ENTREZID")
#View(PATH_mmu)
ECM <- unique(PATH_mmu[grep("ECM",PATH_mmu$DESCRPTION),]$SYMBOL)
ECM
ROS <- unique(PATH_mmu[grep("oxygen species",PATH_mmu$DESCRPTION),]$SYMBOL)
ROS
saveRDS(ROS,"./ROS.rds")
saveRDS(ECM,"./ECM.rds")
SC <- readRDS("./SC.rds")
gc()
ElbowPlot(SC, ndims=50, reduction="pca") 
SC <- FindNeighbors(object = SC, dims = 1:22)
SC <- FindClusters(object = SC, resolution = 0.6)
SC <- RunUMAP(object = SC, dims = 1:22)
DimPlot(object = SC, reduction = "umap",label = T,label.size = 6)+scale_color_simpsons()
DotPlot(object =SC, features = ECM)

DotPlot(object =SC, features = c("Il1a","Il1b","Il6","Tnf","Ccl2","Cxcl2","Cxcl10","Mmp3","Mmp9",
                                 "Pai1","Vegfa","Fgf2","Tgfb1","Spp1","Serpine1","Pdgfb"))
FeaturePlot(object =SC, features = c("Dcn","Col1a1","Col3a1","Adamts1",
                                     "Cd74","Cxcl10","Acta2","Actg2","Ccn2","Cd52"),order=T,reduction = "umap")
FeaturePlot(object =SC, features = c("Fgf2","Tgfb1","Spp1","Serpine1","Pdgfb"),order=T,reduction = "umap")

# Matrix fibroblast   0 1 2 3 5 7 10 11
# Remodeling myofibroblast Adamts1+ 4 9
# SASP fibroblast Spp1+ 8
# Active myofibroblast Acta2+ Actg2+ 6
new.cluster.ids <- c('Matrix fibroblast Dcn+/Cols+',
                     'Matrix fibroblast Dcn+/Cols+',
                     'Matrix fibroblast Dcn+/Cols+',
                     "Matrix fibroblast Dcn+/Cols+",
                     'Remodeling myofibroblast Adamts1+',
                     'Matrix fibroblast Dcn+/Cols+',
                     'Active myofibroblast Acta2+ Actg2+',
                     'Matrix fibroblast Dcn+/Cols+',
                     'SASP fibroblast Spp1+',
                     "Remodeling myofibroblast Adamts1+",
                     'Matrix fibroblast Dcn+/Cols+',
                     'Matrix fibroblast Dcn+/Cols+')
names(x=new.cluster.ids)=levels(x = SC)
SC <- RenameIdents(object =SC, new.cluster.ids)
DimPlot(object = SC, reduction = "umap",label = F,label.size = 6)+scale_color_simpsons()

SC$celltype <-SC@active.ident
SC$age <- substr(SC$group_age,1,2)
prop_test <- sc_utils(SC)
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "M6", sample_2 = "M9",
  sample_identity = "age"
)

permutation_plot2(prop_test,log2FD_threshold =0,FDR_threshold = 0.05)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

#################### 卵丘扩张因子 ######################
library(limma)
M9<- readRDS("./M9_obj.rds")
data <- M9@assays$RNA@data %>% as.data.frame(.)%>%
  .[rownames(.)%in%c("Kitl","Amh","Cx43","Fshr","Inha","Gdf9","Bmp15","Kit","Foxl2",
                     "Notch1","Egfr","Fgfr2","Igf1","Igf1r","Vcam1","Adm","Hbegf",
                     "Cxcr4","S1pr1","Nlrp5"),]

ggdata <- melt(data)
ggdata$diets <- paste0(strsplit2(ggdata$variable,"-")[,2],"RD")

ggdata$diets <- factor(ggdata$diets,levels = c("BRD","WRD","PRD"))
ggplot(ggdata,aes(x=diets,y=log2(value+1),fill=diets,))+#stat_boxplot(geom = "errorbar",width=0.6)+
  geom_boxplot(width=0.6)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  xlab("")+ylab("Relative Expression")+
  geom_signif(comparisons = list(c("BRD","WRD"),c("WRD","PRD")),map_signif_level = T,step_increase = 0.1,
              textsize =5,tip_length = 0,vjust = 0.2)+scale_color_igv()+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())





#Kitl（Stem Cell Factor）
#Amh（Anti-Müllerian Hormone）
#Cx43（Connexin 43）
#Fshr（Follicle-Stimulating Hormone Receptor）
#Inha（Inhibin Alpha Subunit）
#Gdf9（Growth Differentiation Factor 9）
#Bmp15（Bone Morphogenetic Protein 15）
#Kit（Kit proto-oncogene receptor tyrosine kinase）
#Foxl2（Forkhead Box L2）
#Notch1（Notch Receptor 1）
#Egfr（Epidermal Growth Factor Receptor）
#Fgfr2（Fibroblast Growth Factor Receptor 2）
#Igf1（Insulin-like Growth Factor 1）
#Igf1r（Insulin-like Growth Factor 1 Receptor）
#Vcam1（Vascular Cell Adhesion Molecule 1）
#Adm（Adrenomedullin）
#Hbegf（Heparin-binding EGF-like growth factor）
#Cxcr4（C-X-C chemokine receptor type 4）
#S1pr1（Sphingosine-1-phosphate receptor 1）
#Nlrp5（Nucleotide-binding oligomerization domain-like receptor protein 5）



#################### ECM ###############################
data <- M9@assays$RNA@data %>% as.data.frame(.)%>%
  .[rownames(.)%in%ECM,]

ggdata <- melt(data)
ggdata$diets <- paste0(strsplit2(ggdata$variable,"-")[,2],"RD")

ggdata$diets <- factor(ggdata$diets,levels = c("BRD","WRD","PRD"))
ggplot(ggdata,aes(x=diets,y=log(value+1),fill=diets,))+#stat_boxplot(geom = "errorbar",width=0.6)+
  geom_violin(width=0.6)+theme_classic()+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = element_text(size=14),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+
  xlab("")+ylab("Relative Expression")+
  geom_signif(comparisons = list(c("BRD","WRD"),c("WRD","PRD"),c("BRD","PRD")),map_signif_level = T,step_increase = 0.1,
              textsize =5,tip_length = 0,vjust = 0.2)+scale_color_igv()+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())

library(reshape2)
library(ggplot2)
library(ggridges)


ggplot(data_melt, aes(x = value , y = variable , fill = variable)) +
  geom_density_ridges(alpha = 0.5) +
  theme_ridges() +
  theme(legend.position = "none")

########################## Cell interaction 2. iTALK #########################
GC1 <- readRDS("./GC1.rds")
M1 <- readRDS("./M1.rds")
df1=as.data.frame(GC1@active.ident)
GC1@meta.data$celltype=df1$`GC1@active.ident`
df1=as.data.frame(M1@active.ident)
M1@meta.data$celltype=df1$`M1@active.ident`

GC_M <- merge(GC1,M1)
unique(GC_M$celltype)
rm(GC1,M1)
gc()

library("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = rownames(GC_M), #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)

saveRDS(MtoH,"./MtoH.rds")

MtoH <- readRDS("./MtoH.rds")
############################### GC_M 6M #############################
cell <- GC_M@meta.data[grep("M6",GC_M$Diets),]%>%rownames(.)
GC_M_6 <- subset(GC_M,cells=cell)

data <- as.matrix(GC_M_6@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_M6_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_6$celltype

meta_data <- cbind(rownames(GC_M_6@meta.data), GC_M_6@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M6_GCM_meta_data.txt', sep='\t', quote=F)
#View(meta_data)



############################### GC_M 9M #############################
cell <- GC_M@meta.data[grep("M9",GC_M$Diets),]%>%rownames(.)
GC_M_9 <- subset(GC_M,cells=cell)
GC_M_9$Diets <- factor(GC_M_9$Diets,levels = c("M9-BRD","M9-WRD","M9-PRD"))
GC_M_9$celltype <- factor(GC_M_9$celltype,levels = c("Mitotic GC","Preantral GC","Antral GC","Atretic GC",
                                                     "Mo","M1","M1M2","M2"))
VlnPlot(GC_M_9, features = c("Vegfa","Itgav"),cols = c("#F8766D","#078992","#CD9600"),
        group.by = "celltype",split.by = "Diets",pt.size = 0,ncol = 1)


AverageExp<-AverageExpression(GC_M_9)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap::pheatmap(coorda$r,fontsize=16,border_color="black",
                   color = colorRampPalette(c( "white", "firebrick3"))(100))


data <- as.matrix(GC_M_9@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_M9_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)

meta_data <- cbind(rownames(GC_M_9@meta.data), GC_M_9@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_M9_GCM_meta_data.txt', sep='\t', quote=F)
rm(data,GC_M_6,GC_M_9)
gc()

############################### GC_M BRD #############################
cell <- GC_M@meta.data[grep("M9-BRD",GC_M$Diets),]%>%rownames(.)
GC_M_BRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_BRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_BRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_BRD$celltype

meta_data <- cbind(rownames(GC_M_BRD@meta.data), GC_M_BRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_BRD_GCM_meta_data.txt', sep='\t', quote=F)



############################### GC_M WRD #############################
cell <- GC_M@meta.data[grep("M9-WRD",GC_M$Diets),]%>%rownames(.)
GC_M_WRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_WRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_WRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_WRD$celltype

meta_data <- cbind(rownames(GC_M_WRD@meta.data), GC_M_WRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_WRD_GCM_meta_data.txt', sep='\t', quote=F)


############################### GC_M PRD #############################
cell <- GC_M@meta.data[grep("M6-PRD",GC_M$Diets),]%>%rownames(.)
GC_M_PRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_PRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_PRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_PRD$celltype

meta_data <- cbind(rownames(GC_M_PRD@meta.data), GC_M_PRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_PRD_GCM_meta_data.txt', sep='\t', quote=F)

############################### GC_M BRD #############################
cell <- GC_M@meta.data[grep("M6-BRD",GC_M$Diets),]%>%rownames(.)
GC_M_BRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_BRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_6M_BRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_BRD$celltype

meta_data <- cbind(rownames(GC_M_BRD@meta.data), GC_M_BRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_6M_BRD_GCM_meta_data.txt', sep='\t', quote=F)



############################### GC_M WRD #############################
cell <- GC_M@meta.data[grep("M6-WRD",GC_M$Diets),]%>%rownames(.)
GC_M_WRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_WRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_6M_WRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_WRD$celltype

meta_data <- cbind(rownames(GC_M_WRD@meta.data), GC_M_WRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_6M_WRD_GCM_meta_data.txt', sep='\t', quote=F)


############################### GC_M PRD #############################
cell <- GC_M@meta.data[grep("M6-PRD",GC_M$Diets),]%>%rownames(.)
GC_M_PRD <- subset(GC_M,cells=cell)
data <- as.matrix(GC_M_PRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_6M_PRD_GCM_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
GC_M_PRD$celltype

meta_data <- cbind(rownames(GC_M_PRD@meta.data), GC_M_PRD@meta.data[,'celltype', drop=F])  
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
write.table(meta_data, 'cellphonedb_6M_PRD_GCM_meta_data.txt', sep='\t', quote=F)

###################### BRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/BRD_6M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))


###################### WRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/WRD_6M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,120))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))


###################### PRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/PRD_6M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,120))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))




###################### BRD ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/BRD_9M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))


###################### WRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/WRD_9M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))

###################### PRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/PRD_9M/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))
###################### 9M  ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/M9//pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))

###################### 6M  ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/GC_M/M6//pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)

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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,180))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)+
  scale_x_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))+
  scale_y_discrete(limits=c("Mitotic_GC","Preantral_GC","Antral_GC","Atretic_GC","Mo","M1","M1M2","M2"))


####################### cellphone ############################

###################### BRD ###############################
library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_BRD/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))

head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,250))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)


###################### WRD ###############################

library(stringr)
pvalues=read.table("/home/tuyx/scRNA_aging/output/cellphonedb/M9_WRD/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
colnames(pvalues)
colnames(pvalues) <- gsub("SC.TC","SCTC",colnames(pvalues))

head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
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
statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,230))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)
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
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,250))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)
p2

###################### 6M ###############################
M6 <- readRDS("./M6_BP.rds")
unique(M6$abbreviation)
M6$abbreviation2 <- gsub("FBC","SCTC",M6$abbreviation)%>%gsub("TCC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- gsub("SMC","SCTC",M6$abbreviation2)%>%gsub("NRC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- factor(M6$abbreviation2,levels = c("STC",
                                                       "GC","SCTC","EC","EDC",
                                                       "GLC","M","DC","B","NK","T","ILC"))
unique(M6$abbreviation2)
M6@meta.data$age <- substr(M6$group_age,1,2)
M6$Group <- paste0(M6$age,"-",M6$abbreviation2)
M6@active.ident <- M6$abbreviation2 

meta <- M6@meta.data[grep("BRD",M6$Diets),]
M6_BRD <- subset(M6,cells=rownames(meta))
meta <- M6@meta.data[grep("WRD",M6$Diets),]
M6_WRD <- subset(M6,cells=rownames(meta))
meta <- M6@meta.data[grep("PRD",M6$Diets),]
M6_PRD <- subset(M6,cells=rownames(meta))
rm(M6)
gc()
saveRDS(M6_BRD,"./M6_BRD.rds")
saveRDS(M6_WRD,"./M6_WRD.rds")
saveRDS(M6_PRD,"./M6_PRD.rds")
MtoH <- readRDS("./MtoH.rds")
cell <- M6_BRD@meta.data[grep("M6-BRD",M6_BRD$Diets),]%>%rownames(.)
M6_BRD<- subset(M6_BRD,cells=cell)
data <- as.matrix(M6_BRD@assays$RNA@data)
data <- data.frame(gene=rownames(data), data, check.names = F)
data$Gene <- MtoH[match(data$gene, MtoH[,1]),2]
data <- subset(data, Gene!='NA')
data <- dplyr::select(data, Gene, everything())
data <- data[, !(colnames(data) %in% 'gene')]
rownames(data) <- NULL
data[1:5,1:5]
write.table(data, 'cellphonedb_M6_BRD_count.txt',col.names=T,row.names=F, sep='\t', quote=F)
M6_BRD$abbreviation2

meta_data <- cbind(rownames(M6_BRD@meta.data), M6_BRD@meta.data[,'abbreviation2', drop=F]) 
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" 
head(meta_data)
write.table(meta_data, 'cellphonedb_M6_BRD_meta_data.txt', sep='\t', quote=F)





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
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,190))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number)


saveRDS(GC_M,"./GC_M.rds")
devtools::load_all("~/monocle")

library(monocle)
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
# iTALK 要求的矩阵: 行为细胞，列为基因
GC_M <- readRDS("./GC_M.rds")
MtoH <- readRDS("./MtoH.rds")
meta <- GC_M@meta.data
count_O <- as.data.frame(GC_M@assays$RNA@data)[,grep("M9",rownames(meta))]
count_O <- data.frame(gene=rownames(count_O), count_O, check.names = F)
count_O$Gene <- MtoH[match(count_O$gene, MtoH[,1]),2]
count_O <- subset(count_O, Gene!='NA')
count_O <- dplyr::select(count_O, Gene, everything())
count_O <- count_O[, !(colnames(count_O) %in% 'gene')]
count_O[1:5,1:5]
count_O <- count_O[!duplicated(count_O$Gene),]
rownames(count_O) <- count_O$Gene
count_O$Gene <- NULL
iTalk_O <- as.data.frame(t(count_O))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_O$cell_type <- meta[grep("M9",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)

iTalk_O$cell_type <- factor(iTalk_O$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                        "Atr GC","Mo","M1", "M1M2", "M2"))


#iTalk_O <- iTalk_O[iTalk_O$cell_type%in%c(c("Mitotic GC","Preantral GC","Antral GC",
#                                            "Atretic GC","Monocytes")),]
unique(iTalk_O$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_O, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')

cell_types <- unique(iTalk_O$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
#View(iTalk_res)
res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]

res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)
dim(res)

iTalk_O <- res
res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
iTalk_O_30 <- res1
LRPlot(iTalk_O_30,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)


count_Y <- as.data.frame(GC_M@assays$RNA@data)[,grep("M6",rownames(meta))]
count_Y <- data.frame(gene=rownames(count_Y), count_Y, check.names = F)
count_Y$Gene <- MtoH[match(count_Y$gene, MtoH[,1]),2]
count_Y <- subset(count_Y, Gene!='NA')
count_Y <- dplyr::select(count_Y, Gene, everything())
count_Y <- count_Y[, !(colnames(count_Y) %in% 'gene')]
count_Y[1:5,1:5]
count_Y <- count_Y[!duplicated(count_Y$Gene),]
rownames(count_Y) <- count_Y$Gene
count_Y$Gene <- NULL
iTalk_Y <- as.data.frame(t(count_Y))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_Y$cell_type <- meta[grep("M6",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)

iTalk_Y$cell_type <- factor(iTalk_Y$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                        "Atr GC","Mo","M1", "M1M2", "M2"))


#iTalk_O <- iTalk_O[iTalk_O$cell_type%in%c(c("Mitotic GC","Preantral GC","Antral GC",
#                                            "Atretic GC","Monocytes")),]
unique(iTalk_Y$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_Y, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')

cell_types <- unique(iTalk_Y$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
#View(iTalk_res)
res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]

res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]

LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)
iTalk_Y <- res
res2 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:200,]
iTalk_Y_200 <- res2
dim(iTalk_Y_200)
b=unique(iTalk_Y_200$ligand,iTalk_Y_200$receptor)
ncbi2<- AnnotationDbi::select(org.Hs.eg.db,keys=b,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG2 <- enrichKEGG(gene =ncbi2$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff =0.05)

View(KEGG2@result)




count_BRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M6-B",rownames(meta))]
count_BRD <- data.frame(gene=rownames(count_BRD), count_BRD, check.names = F)
count_BRD$Gene <- MtoH[match(count_BRD$gene, MtoH[,1]),2]
count_BRD <- subset(count_BRD, Gene!='NA')
count_BRD <- dplyr::select(count_BRD, Gene, everything())
count_BRD <- count_BRD[, !(colnames(count_BRD) %in% 'gene')]
count_BRD[1:5,1:5]
count_BRD <- count_BRD[!duplicated(count_BRD$Gene),]
rownames(count_BRD) <- count_BRD$Gene
count_BRD$Gene <- NULL
iTalk_BRD <- as.data.frame(t(count_BRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_BRD$cell_type <- meta[grep("M6-B",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_BRD$cell_type <- factor(iTalk_BRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_BRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_BRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_BRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)

count_WRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M6-W",rownames(meta))]
count_WRD <- data.frame(gene=rownames(count_WRD), count_WRD, check.names = F)
count_WRD$Gene <- MtoH[match(count_WRD$gene, MtoH[,1]),2]
count_WRD <- subset(count_WRD, Gene!='NA')
count_WRD <- dplyr::select(count_WRD, Gene, everything())
count_WRD <- count_WRD[, !(colnames(count_WRD) %in% 'gene')]
count_WRD[1:5,1:5]
count_WRD <- count_WRD[!duplicated(count_WRD$Gene),]
rownames(count_WRD) <- count_WRD$Gene
count_WRD$Gene <- NULL
iTalk_WRD <- as.data.frame(t(count_WRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_WRD$cell_type <- meta[grep("M6-W",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_WRD$cell_type <- factor(iTalk_WRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_WRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_WRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_WRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)


count_PRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M6-P",rownames(meta))]
count_PRD <- data.frame(gene=rownames(count_PRD), count_PRD, check.names = F)
count_PRD$Gene <- MtoH[match(count_PRD$gene, MtoH[,1]),2]
count_PRD <- subset(count_PRD, Gene!='NA')
count_PRD <- dplyr::select(count_PRD, Gene, everything())
count_PRD <- count_PRD[, !(colnames(count_PRD) %in% 'gene')]
count_PRD[1:5,1:5]
count_PRD <- count_PRD[!duplicated(count_PRD$Gene),]
rownames(count_PRD) <- count_PRD$Gene
count_PRD$Gene <- NULL
iTalk_PRD <- as.data.frame(t(count_PRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_PRD$cell_type <- meta[grep("M6-P",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_PRD$cell_type <- factor(iTalk_PRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_PRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_PRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_PRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)



count_BRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M9-B",rownames(meta))]
count_BRD <- data.frame(gene=rownames(count_BRD), count_BRD, check.names = F)
count_BRD$Gene <- MtoH[match(count_BRD$gene, MtoH[,1]),2]
count_BRD <- subset(count_BRD, Gene!='NA')
count_BRD <- dplyr::select(count_BRD, Gene, everything())
count_BRD <- count_BRD[, !(colnames(count_BRD) %in% 'gene')]
count_BRD[1:5,1:5]
count_BRD <- count_BRD[!duplicated(count_BRD$Gene),]
rownames(count_BRD) <- count_BRD$Gene
count_BRD$Gene <- NULL
iTalk_BRD <- as.data.frame(t(count_BRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_BRD$cell_type <- meta[grep("M9-B",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_BRD$cell_type <- factor(iTalk_BRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_BRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_BRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_BRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)

count_WRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M9-W",rownames(meta))]
count_WRD <- data.frame(gene=rownames(count_WRD), count_WRD, check.names = F)
count_WRD$Gene <- MtoH[match(count_WRD$gene, MtoH[,1]),2]
count_WRD <- subset(count_WRD, Gene!='NA')
count_WRD <- dplyr::select(count_WRD, Gene, everything())
count_WRD <- count_WRD[, !(colnames(count_WRD) %in% 'gene')]
count_WRD[1:5,1:5]
count_WRD <- count_WRD[!duplicated(count_WRD$Gene),]
rownames(count_WRD) <- count_WRD$Gene
count_WRD$Gene <- NULL
iTalk_WRD <- as.data.frame(t(count_WRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_WRD$cell_type <- meta[grep("M9-W",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_WRD$cell_type <- factor(iTalk_WRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_WRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_WRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_WRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)


count_PRD <- as.data.frame(GC_M@assays$RNA@data)[,grep("M9-P",rownames(meta))]
count_PRD <- data.frame(gene=rownames(count_PRD), count_PRD, check.names = F)
count_PRD$Gene <- MtoH[match(count_PRD$gene, MtoH[,1]),2]
count_PRD <- subset(count_PRD, Gene!='NA')
count_PRD <- dplyr::select(count_PRD, Gene, everything())
count_PRD <- count_PRD[, !(colnames(count_PRD) %in% 'gene')]
count_PRD[1:5,1:5]
count_PRD <- count_PRD[!duplicated(count_PRD$Gene),]
rownames(count_PRD) <- count_PRD$Gene
count_PRD$Gene <- NULL
iTalk_PRD <- as.data.frame(t(count_PRD))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_PRD$cell_type <- meta[grep("M9-P",rownames(meta)),]$celltype%>%gsub("Mitotic GC","Mit GC",.)%>%
  gsub("Preantral GC","Pre GC",.)%>%gsub("Antral GC","Ant GC",.)%>%gsub("Atretic GC","Atr GC",.)
iTalk_PRD$cell_type <- factor(iTalk_PRD$cell_type,levels= c("Mit GC","Pre GC","Ant GC",
                                                            "Atr GC","Mo","M1", "M1M2", "M2"))
unique(iTalk_PRD$cell_type)
my10colors <-c('#2CA02C', '#FF7F0E', '#1F77B4', '#D62728', '#FF6F00', '#008EA0', '#C71000','#8A4198')
highly_exprs_genes <- rawParse(iTalk_PRD, top_genes=50, stats="mean")

# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_PRD$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- levels(cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

res <- iTalk_res[iTalk_res$cell_from%in%c("Mo","M1", "M1M2", "M2"),]
res <- res[res$cell_to%in%c("Mit GC","Pre GC","Ant GC",
                            "Atr GC"),]
res <- res[res$comm_type%in%c('growth factor'),]


res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,
       link.arr.width=res1$cell_to_mean_exprs)


################ STC ###################
M6 <- readRDS("./M6_obj.rds")
unique(M6$abbreviation)
M6$abbreviation2 <- gsub("FBC","SCTC",M6$abbreviation)%>%gsub("TCC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- gsub("SMC","SCTC",M6$abbreviation2)%>%gsub("NRC","SCTC",.)
unique(M6$abbreviation2)
M6$abbreviation2 <- factor(M6$abbreviation2,levels = c("STC",
                                                       "GC","SCTC","EC","EDC",
                                                       "GLC","M","DC","B","NK","T","ILC"))
unique(M6$abbreviation2)
M6@meta.data$age <- substr(M6$group_age,1,2)
M6$Group <- paste0(M6$age,"-",M6$abbreviation2)
M6@active.ident <- M6$abbreviation2


STGC <- subset(M6,idents=c("STC","GC"))
rm(M6)
gc()
head(STGC@meta.data)
STGC <- RunPCA(STGC, features = VariableFeatures(STGC))
DimPlot(STGC, reduction = "pca", group.by="orig.ident")
ElbowPlot(STGC, ndims=50, reduction="pca") 
STGC <- FindNeighbors(object = STGC, dims = 1:8)
STGC <- FindClusters(object = STGC, resolution = 0.2)
STGC <- RunUMAP(object = STGC, dims = 1:8)
DimPlot(object = STGC, reduction = "umap",label = T,raster=FALSE)
STGC <- RunTSNE(object = STGC, dims = 1:8)
DimPlot(object = STGC, reduction = "tsne",label = T,label.size = 6,raster=FALSE)
DimPlot(object = STGC, reduction = "tsne",label = F,group.by="orig.ident") 
meta=STGC@meta.data 
library(SingleR)
ImmGen.se=ImmGenData() #(鼠)
Mouse.se=MouseRNAseqData() #(鼠)

STGC_for_SingleR <- GetAssayData(STGC, slot="data") ##获取标准化矩阵
STGC.hesc <- SingleR(test = STGC_for_SingleR, ref = Mouse.se, labels = Mouse.se$label.main) 
STGC.hesc2 <- SingleR(test = STGC_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main) 
STGC@meta.data$SingleR_Mouse <- STGC.hesc$labels
STGC@meta.data$SingleR_ImmGen <- STGC.hesc2$labels
DimPlot(object = STGC, reduction = "tsne",label =T,group.by = "SingleR_Mouse",raster=FALSE)
DimPlot(object = STGC, reduction = "tsne",label =T,group.by = "SingleR_ImmGen",raster=FALSE)
DimPlot(object = STGC, reduction = "umap",label =T,group.by = "SingleR_Mouse",raster=FALSE)
DimPlot(object = STGC, reduction = "umap",label =T,group.by = "SingleR_ImmGen",raster=FALSE,label.size = 5)
DimPlot(object = STGC, reduction = "tsne",label = T,raster=FALSE)
DimPlot(object = STGC, reduction = "umap",label = T,raster=FALSE,label.size = 5)

STGC@meta.data$SingleR_ImmGen[grep("Stem cells",STGC@meta.data$SingleR_ImmGen,invert = T)] <- "Granulosa cells"
DimPlot(object = STGC, reduction = "umap",label =T,group.by = "SingleR_ImmGen",
        raster=FALSE,label.size = 5,cols=c("#F99502","#569F86"))

ident <- STGC@meta.data$SingleR_ImmGen
names(ident) <- rownames(STGC@meta.data)
STGC@active.ident <- factor(ident,levels = c("Stem cells","Granulosa cells"))


data <- as(as.matrix(STGC@assays$RNA@data), 'sparseMatrix')
pd <- STGC@meta.data
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

sc_cds$group <- sc_cds$SingleR_ImmGen
?plot_cell_trajectory
plot_cell_trajectory(sc_cds, color_by ="group",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))
saveRDS(sc_cds,"STGC1_sc_cds.rds")

diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)
sc_cds

plot_cell_trajectory(sc_cds, color_by ="group",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")

saveRDS(sc_cds,"STGC2_sc_cds.rds")
plot_cell_trajectory(sc_cds, color_by ="State",cell_size = 1)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),
        legend.text = element_text(size=14))+theme(legend.position = "right")


sc_cds <- readRDS("./STGC1_sc_cds.rds")
sc_cds <- orderCells(sc_cds,reverse = T)
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

rm(STGC)
gc()
colnames(sc_cds)
sc_cds <- sc_cds[,sc_cds$State==1]
diff_celltype<- differentialGeneTest(sc_cds,fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_celltype, qval < 0.01))
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = T)
sc_cds
plot_cell_trajectory(sc_cds, color_by ="group",cell_size = 1)+
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
saveRDS(sc_cds,"./STGC_select1.rds")
?BEAM
BEAM_res=BEAM(sc_cds[ordering_genes,],branch_point = 2,cores = 1)
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_STGC_select.rds")

plot_genes_branched_heatmap(sc_cds[row.names(subset(BEAM_res,qval<1e-4)),],
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


p <- plot_pseudotime_heatmap(sc_cds[row.names(subset(diff_celltype, qval <1e-4)),],
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
devtools::install_github("junjunlab/ClusterGVis")
devtools::install_local("~/ClusterGVis.tar.gz")
library(ClusterGVis)
devtools::install_github("junjunlab/jjAnno")
library(jjAnno)
saveRDS(diff_celltype,"./diff_celltype.rds")

df <- plot_pseudotime_heatmap2(sc_cds[row.names(subset(diff_celltype, qval <1e-4)),],
                               num_clusters = 4,
                               cores = 1)

visCluster(object = p,plot.type = "both")
library(dplyr)

df_cluster_monocle <- read.csv("~/scRNA_aging/output/df_cluster_monocle.csv", row.names=1)
M9 <- readRDS("./M9_obj.rds")
unique(M9$abbreviation)
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
ggdata$group <- factor(ggdata$group,levels = c("M9-BRD","M9-WRD","M9-PRD"))
head(ggdata)
df_cluster_monocle$SYMBOL <- df_cluster_monocle$gene
data <- merge(df_cluster_monocle,ggdata,all.x=T,by="SYMBOL")%>%na.omit(.)
data1 <- data[data$value>0.1,]
data1$group <- substr(data1$group,4,6)
data1$group <- factor(data1$group,levels = c("BRD","WRD","PRD"))
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

ggplot(data2,aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("BRD","WRD"),c("WRD","PRD"),c("BRD","PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0.1,test="t.test")+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 22),axis.text =element_text(size = 22, color = 'black'))+	
  theme(axis.text = element_text(size = 22),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 22),
        strip.text.x = element_text(size=16),strip.background = element_blank())+
  facet_wrap(~cluster,ncol=4,scales = "free_y")

all <- readRDS("./all_obj.rds")
sc_cds <- readRDS("./STGC_select1.rds")
unique(all$abbreviation)
all$abbreviation2 <- gsub("SC","SCTC",all$abbreviation)%>%gsub("TCC","SCTC",.)
unique(all$abbreviation2)

all$abbreviation2 <- factor(all$abbreviation2,levels = c("STC",
                                                       "GC","SCTC","EC","EDC",
                                                       "GLC","M","DC","B","NK","T","ILC"))
unique(all$abbreviation2)
all@meta.data$age <- substr(all$group_age,1,2)
all$Group <- paste0(all$age,"-",all$abbreviation2)
all@active.ident <- all$abbreviation2
STGC <- subset(all,idents=c("STC","GC"))
rm(all)
gc()
STGC
STGC <- subset(STGC,cells=colnames(sc_cds))
data <- as.data.frame(STGC@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
rm(data,STGC)
gc()
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-BRD","M6-WRD","M6-PRD",
                                               "M9-BRD","M9-WRD","M9-PRD"))
head(ggdata)
df_cluster_monocle$SYMBOL <- df_cluster_monocle$gene
data <- merge(df_cluster_monocle,ggdata,all.x=T,by="SYMBOL")%>%na.omit(.)
rm(ggdata)
gc()
data1 <- data[data$value>1,]
#data1$group <- substr(data1$group,4,6)
#data1$group <- factor(data1$group,levels = c("BRD","WRD","PRD"))
data1$cluster <- paste0("cluster",data1$cluster)
library(dplyr)



cluster <- data1[data1$cluster=="cluster1",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data1[data1$cluster=="cluster2",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data1[data1$cluster=="cluster3",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))
cluster <- data1[data1$cluster=="cluster4",]
cluster %>% group_by(group) %>%dplyr::summarise(n=mean(value))


ggplot(data1,aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-BRD","M9-BRD"),c("M6-WRD","M9-WRD"),c("M6-PRD","M9-PRD"),
                                 c("M9-BRD","M9-WRD"),c("M9-WRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.6,textsize = 4,
              tip_length = 0,vjust = 0)+ylab("Relative expression")+xlab("")+	
  theme(axis.title =element_text(size = 18),axis.text =element_text(size = 18, color = 'black'))+	
  theme(axis.text = element_text(size = 18),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = 'none')+ggtitle("")+ 
  theme(plot.title = element_text(size = 18),
        strip.text.x = element_text(size=16),strip.background = element_blank())+
  facet_wrap(~cluster,ncol=4,scales = "free_y")


plotdf=pData(sc_cds)
plotdf$group
ggplot(plotdf, aes(x=Pseudotime,y=group,fill=group))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(axis.text = element_text(colour = "black",size = 16),
    panel.grid = element_blank()
  )


ggplot(plotdf, aes(x=Pseudotime,y=group,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(axis.text = element_text(colour = "black",size = 16),
    panel.grid = element_blank()
  )



disp_table <- dispersionTable(sc_cds)
# first method for genes selected
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
sc_cds <- setOrderingFilter(sc_cds, unsup_clustering_genes$gene_id)
sc_cds<- reduceDimension(sc_cds, max_components = 2,
                         method = 'DDRTree')
sc_cds <- orderCells(sc_cds,reverse = F)

plot_cell_trajectory(sc_cds, color_by ="group",cell_size = 0.5)+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 20, color = 'black'))+
  theme(axis.text.x = element_text(angle =0, hjust = 0.5),legend.text = element_text(size=14))

sc_cds <- reduceDimension(sc_cds, max_components = 2, num_dim = 2,
                        reduction_method = 'tSNE',
                        verbose = T)
sc_cds <- clusterCells(sc_cds, num_clusters = 2)
?plot_cell_clusters
plot_cell_clusters(sc_cds, 1, 2, color = "group",cell_size = 0.5)


######################### Gene sets ###########################
# install packages from CRAN
cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")
if (!requireNamespace(cran.packages, quietly = TRUE)) {
  install.packages(cran.packages, ask = F, update = F)
}

# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")
if (!requireNamespace(bioconductor.packages, quietly = TRUE)) {
  BiocManager::install(bioconductor.packages, ask = F, update = F)
}
# install packages from Github
if (!requireNamespace("UCell", quietly = TRUE)) {
  devtools::install_github("carmonalab/UCell")
}
if (!requireNamespace("irGSEA", quietly = TRUE)) {
  devtools::install_github("chuiqin/irGSEA")
}
