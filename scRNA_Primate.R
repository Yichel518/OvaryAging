setwd("/home/tuyx/scRNA_aging/")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
# Bulk
Homo_mouse_human_macaque <- read.delim("~/scRNA_aging/Homo_mouse_human_macaque.txt")
human_mouse_exp_cor <- read.csv("~/scRNA_aging/human_mouse_exp_cor.txt", sep="")
human_mouse_exp_cor <- human_mouse_exp_cor[abs(human_mouse_exp_cor$log2FoldChange)>0.5,]
dim(human_mouse_exp_cor)
#View(human_mouse_exp_cor)
# Aging genes
head(human_mouse_exp_cor)
human_gene <- human_mouse_exp_cor[human_mouse_exp_cor$direction%in%c("Aging"),]$gene_id
age_gene <- Homo_mouse_human_macaque[Homo_mouse_human_macaque$Human.gene.stable.ID%in%human_gene,]$Macaque.gene.name
age_gene <- intersect(human_mouse_exp_cor$Gene.name,age_gene)
human_gene <- human_mouse_exp_cor[human_mouse_exp_cor$direction%in%c("Early"),]$gene_id
early_gene <- Homo_mouse_human_macaque[Homo_mouse_human_macaque$Human.gene.stable.ID%in%human_gene,]$Macaque.gene.name
length(age_gene)
length(early_gene)
# scRNA
GSE130664_barcode_information <- read.delim("~/scRNA_aging/GSE130664_barcode_information.txt")
GSE130664_merge_UMI_count <- read.delim("~/scRNA_aging/GSE130664_merge_UMI_count.txt", row.names=1)
View(GSE130664_barcode_information)
colnames(GSE130664_merge_UMI_count)
ovary <- CreateSeuratObject(counts = GSE130664_merge_UMI_count, project = "ovary", min.cells = 3, 
                                      min.features =200)
ovary[["percent.mito"]] <- PercentageFeatureSet(object = ovary, pattern = "^MT")
ovary@assays$RNA@data[grep("MT",rownames(ovary@assays$RNA@data)),]%>%rownames()

head(ovary[["percent.mito"]])
ovary.mt=ovary[["percent.mito"]]
head(x = ovary@meta.data, 5)
VlnPlot(object = ovary, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
plot1 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# Filter cells
ovary <- subset(ovary, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 500000)
# Normalize the data
ovary <- NormalizeData(object = ovary, normalization.method = "LogNormalize", scale.factor = 1e4)
# Find top 2000 variable genes
ovary <- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2000)
# ID top 10 variable genes
top10 <- head(VariableFeatures(ovary), 10)
# Plot variable genes
plot1 <- VariableFeaturePlot(ovary)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
# Scale data based on all genes
all.genes <- rownames(ovary)
ovary <- ScaleData(ovary, features = all.genes)
# Use ID'd 2000 Variable genes to run PCA analysis
ovary <- RunPCA(ovary, features = VariableFeatures(object = ovary))
# Plot heatmaps based on PCs
DimHeatmap(ovary, dims = 1:20, cells = 500, balanced = TRUE)
# Elbow plot and jackstraw plots are used to determine the number of PC to use 
ElbowPlot(ovary,ndims = 50)
# Elbow plot and jackstraw plots are used to determine the number of PC to use 
plan("multiprocess", workers = 4)
ovary <- JackStraw(ovary, num.replicate = 100)
ovary <- ScoreJackStraw(ovary, dims = 1:20)
JackStrawPlot(ovary, dims = 1:20)

ovary <- FindNeighbors(ovary, dims = 1:20)
ovary <- FindClusters(ovary, resolution = 0.4)
ovary <- RunUMAP(ovary, dims = 1:20)
ovary <- RunTSNE(ovary, dims = 1:20)
DimPlot(ovary, reduction = "umap", pt.size = 1, label = T)
DimPlot(ovary, reduction = "tsne", pt.size = 1, label = T)

# OO 5
DotPlot(object = ovary, features = c("SYCP3","DDX4","GDF9","ZP3"))
# GC 2 4
DotPlot(object = ovary, features = c("AMH","NR5A2","CYP19A1","INHA","WT1"))
# SC 0 11 13
DotPlot(object = ovary, features = c("TCF21","STAR","ALDH1A1","COL1A1","COL1A2"))
# SMC 1  14
DotPlot(object = ovary, features = c("DES","RARB","CSRP1","ACTA2","ACTG2"))
# EC 6
DotPlot(object = ovary, features = c("CDH5","ERG","VWF","RNASE1"))
# M 8 
DotPlot(object = ovary, features = c("CD68","CD14","CD163"))
# NKT 3
DotPlot(object = ovary, features = c("CD3D","KLRB","REL","TRAC","CD3G","KLRD1"))

ovary1 <- ovary
df1=as.data.frame(ovary1@active.ident)
ovary1@meta.data$celltype=paste0(df1$`ovary1@active.ident`,"_",ovary1$group)
idents <- paste0(df1$`ovary1@active.ident`,"_",ovary1$group)
names(idents) <- rownames(ovary1@meta.data)
ovary1@active.ident<-factor(idents)
DotPlot(object = ovary1, features = c("SYCP3","DDX4","GDF9","ZP3",
                                     "AMH","NR5A2","CYP19A1","INHA",
                                     "TCF21","STAR","COL1A1","COL1A2",
                                     "DES","CSRP1","ACTG2",
                                     "CDH5","ERG","VWF","RNASE1",
                                     "CD68","CD14","CD163",
                                     "CD3D","KLRB","TRAC","CD3G","KLRD1"))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.title =element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.title =element_text(size = 18),legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        axis.text =element_text(size = 18, color = 'black'))+
  theme(axis.text.x = element_text(size=18,angle = 90, vjust = 0.25,hjust = 1))+coord_flip()+
  scale_y_discrete(limits=c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
### Aging cell
OO <- FindMarkers(ovary1, ident.1 = "OO_O", ident.2 = "OO_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
GC <- FindMarkers(ovary1, ident.1 = "GC_O", ident.2 = "GC_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
SC <- FindMarkers(ovary1, ident.1 = "SC_O", ident.2 = "SC_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
SMC <- FindMarkers(ovary1, ident.1 = "SMC_O", ident.2 = "SMC_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
EC <- FindMarkers(ovary1, ident.1 = "EC_O", ident.2 = "EC_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
M <- FindMarkers(ovary1, ident.1 = "M_O", ident.2 = "M_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]
NKT <- FindMarkers(ovary1, ident.1 = "NKT_O", ident.2 = "NKT_Y", min.pct = 0.25)%>%.[order(.$avg_log2FC,decreasing = T),]


intersect(rownames(OO),rownames(GC))%>%intersect(.,rownames(SC))%>%intersect(.,rownames(SMC))%>%
  intersect(.,rownames(EC))%>%intersect(.,rownames(M))%>%intersect(.,rownames(NKT))


a=c("OO","GC","SC","SMC","EC","M","NKT")
get(a[1])
BP<- purrr::map(1:length(a),function(i){clusterProfiler::simplify(enrichGO(gene=rownames(get(a[i])),keyType = "SYMBOL",
                                         OrgDb= org.Hs.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))})

View(BP[[3]]@result)









ovary <- subset(x = ovary,idents=c(0,1,2,3,4,5,6,8,11,13,14))




new.cluster.ids <- c("SC","SMC","GC","NKT","GC",
                     "OO","EC","M","SC","SC","SMC")
names(x=new.cluster.ids)=levels(x = ovary)
ovary <- RenameIdents(object = ovary, new.cluster.ids)
DimPlot(ovary, reduction = "tsne", pt.size = 1, label = T,label.size = 6)
df1=as.data.frame(ovary@active.ident)
ovary@meta.data$celltype=df1$`ovary@active.ident`
ovary_new.markers <- FindAllMarkers(object = ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(ovary_new.markers)
write.table(ovary_new.markers,"ovary_new.markers.txt",col.names=T,row.names=T,quote=F,sep="\t")
ovary@meta.data$group <- substr(ovary@meta.data$orig.ident,1,1)
DimPlot(ovary,reduction='tsne',group.by = "group")

ovary <- ScaleData(object = ovary, features = rownames(ovary))
ovary@assays$RNA@scale.data
DoHeatmap(ovary1, features =intersect(c(age_gene,"CDKN2A","CDKN1A","CDKN2D","GPNMB"),rownames(ovary)))+NoLegend()+theme(axis.text = element_text(size = 14))
?DoHeatmap
data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c(c(age_gene,"CDKN2A","CDKN1A","CDKN2D","GPNMB","ARF")),]%>%t()%>%as.data.frame()
cdkn2a$group <- ovary@meta.data$group
cdkn2a$group <- factor(cdkn2a$group,levels = c("Y","O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- melt(cdkn2a)

ggplot(ggdata,aes(x = group,y = log2(value+1),fill=group))+
  stat_boxplot(geom = "errorbar", width=0.1)+
  geom_boxplot(width=0.4)+theme_classic()+
  geom_signif(comparisons = list(c("Y","O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2,na.rm = T)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Aging biomarkers relative expression")+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())






?AverageExpression
averageExp<-apply(ovary@assays$RNA@data,1,function(x){tapply(x,paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group),mean)}) %>% t() %>% as.data.frame()
# Expage <- averageExp$RNA[rownames(averageExp$RNA)%in%monkey_gene,]
Expage <- averageExp[rownames(averageExp)%in%age_gene,]
n=t(scale(t(Expage)))
n[n>2]=2
n[n<-2]=-2
p <- pheatmap::pheatmap(na.omit(n)[,c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O")],
                   show_rownames = F,fontsize = 18,cluster_rows = T,cluster_cols = F,
                   colorRampPalette(c("white","#FEF5B1","#E24C36"))(100))

pheatmap::pheatmap(na.omit(n)[p$tree_row$order,c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O")],
                   show_rownames = F,fontsize = 18,cluster_rows = F,cluster_cols = F,
                   colorRampPalette(c("white","#FEF5B1","#E24C36"))(100),legend = "bottum")

col_fun = colorRamp2(c(-1, 0, 1), c("white","#FEF5B1","#E24C36"))

ha = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:3), 
                                              height = unit(4, "cm")))
v = rnorm(50)
data=na.omit(n)[p$tree_row$order,c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O")]
ha <- HeatmapAnnotation(
  density =anno_density(
    data,
    height = unit(2, "cm"),
    gp = gpar(
      fill =c("#6686c6","#a0c1db","gray","#E65F92","#c77364","#ce8f5c",						
              "#7bac80","#75a5c0","#b5181a","#b72d9e",						
              "#e4cccb","#f6f0eb","#e8c755","#d5d456")
    )
  )
)
ha


Heatmap(data,
        cluster_rows = F,cluster_columns = F,show_row_names = F,col = col_fun,heatmap_legend_param = list(
          title = "scale",
          title_position = "leftcenter-rot"),top_annotation = ha)




ovary@meta.data$group <- factor(ovary@meta.data$group,levels = c("Y","O"))
FeaturePlot(ovary,features = c("CDKN2A"),reduction = "tsne",pt.size = 1.5,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")
FeaturePlot(ovary,features = c("CDKN1A"),reduction = "tsne",pt.size = 1.5,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")
FeaturePlot(ovary,features = c("CDKN2D"),reduction = "tsne",pt.size = 1.5,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")
FeaturePlot(ovary,features = c("GPNMB"),reduction = "tsne",pt.size = 1.5,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")

c("CDKN2A","CDKN1A","CDKN2D","GPNMB")
FeaturePlot(ovary,features = c("MYC","APOE"),reduction = "tsne",pt.size = 1.5,ncol = 4,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")

FeaturePlot(ovary,features = c("ND3","ENSMFAG00000011609"),reduction = "tsne",pt.size = 1.5,ncol = 4,order = TRUE,split.by = "group",label = T)& theme(legend.position = "right")

meta_Y <- ovary@meta.data[ovary@meta.data$group=="Y",]
meta_Y <- data.frame(Cell=rownames(meta_Y),cell_type=meta_Y$celltype)
count_Y <- as.data.frame(ovary@assays$RNA@data)[,meta_Y$Cell]
count_Y <- cbind(Gene=rownames(count_Y),count_Y)
rownames(count_Y) <- NULL

meta_O <- ovary@meta.data[ovary@meta.data$group=="O",]
meta_O <- data.frame(Cell=rownames(meta_O),cell_type=meta_O$celltype)
count_O <- as.data.frame(ovary@assays$RNA@data)[,meta_O$Cell]
count_O <- cbind(Gene=rownames(count_O),count_O)
rownames(count_O) <- NULL
write.table(meta_Y,"meta_Y.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(count_Y,"count_Y.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(meta_O,"meta_O.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(count_O,"count_O.txt",col.names=T,row.names=F,quote=F,sep="\t")
count_O[1:5,1:5]

library(tidyverse)
library(RColorBrewer)
library(scales)

pvalues=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/O_out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.01))
#View(statdf)

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

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,150))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
    axis.text.x.bottom = element_text(),
    panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number,size =5)

pvalues=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/Y_out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
head(pvalues)
pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
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

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=total_number))+geom_tile(color="white")+
  scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,150))+
  theme_minimal()+
  theme(axis.text = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x.bottom = element_text(),
        panel.grid = element_blank()
  )+xlab("")+ylab("")+geom_text(label=statdf$total_number,size =5)

sigmean_Y=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/Y_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
#View(sigmean_Y)
a <- na.omit(sigmean_Y[,c('gene_a',"gene_b","GC.OO")])
colnames(a) <- c('gene_a',"gene_b","level")
a$class <- "GC.OO"
b <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.GC")])
colnames(b) <- c('gene_a',"gene_b","level")
b$class <- "OO.GC"
GCOO_Y <- rbind(a,b)
dim(GCOO_Y)
GCOO_Y$pairs <- paste0(GCOO_Y$gene_a,"_",GCOO_Y$gene_b)

sigmean_O=read.table("/home/tuyx/scRNA_aging/cellphonedb-data-4.1.0/O_out/significant_means.txt",header = T,sep = "\t",stringsAsFactors = F)
c <- na.omit(sigmean_O[,c('gene_a',"gene_b","GC.OO")])
colnames(c) <- c('gene_a',"gene_b","level")
c$class <- "GC.OO"
d <- na.omit(sigmean_Y[,c('gene_a',"gene_b","OO.GC")])
colnames(d) <- c('gene_a',"gene_b","level")
d$class <- "OO.GC"
GCOO_O <- rbind(c,d)
GCOO_O$pairs <- paste0(GCOO_O$gene_a,"_",GCOO_Y$gene_b)

intersect(GCOO_Y$pairs,GCOO_O$pairs)
strsplit2(setdiff(GCOO_O$pairs,GCOO_Y$pairs),"_")[,1]

list_o <- strsplit2(setdiff(GCOO_O$pairs,GCOO_Y$pairs),"_")[,1]
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list_o,
                             columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
BP <- clusterProfiler::simplify(enrichGO(gene=list_o,keyType = "SYMBOL",
                                   OrgDb= org.Hs.eg.db, ont="BP",
                                   pAdjustMethod="BH",pvalueCutoff=0.05,
                                   qvalueCutoff=0.05))

View(BP@result)
write.table(BP@result, 'BP_GC_OO_pairs.txt', sep='\t', quote=F)

KEGG <- enrichKEGG(gene =sym$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff =0.05)
View(KEGG@result)
write.table(KEGG@result, 'KEGG_GC_OO_pairs.txt', sep='\t', quote=F)


list_y <- strsplit2(setdiff(GCOO_Y$pairs,GCOO_O$pairs),"_")[,1]
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list_y,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
BP <- clusterProfiler::simplify(enrichGO(gene=list_y,keyType = "SYMBOL",
                                         OrgDb= org.Hs.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))

View(BP@result)
KEGG <- enrichKEGG(gene =sym$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff =0.05)
View(KEGG@result)
write.table(KEGG@result, 'KEGG_GC_YY_pairs.txt', sep='\t', quote=F)

write.table(BP@result, 'BP_GC_YY_pairs.txt', sep='\t', quote=F)





KEGG_GC_OO_pairs <- read.delim("~/scRNA_aging/KEGG_GC_OO_pairs.txt", row.names=1)


s <- KEGG_GC_OO_pairs[KEGG_GC_OO_pairs$select==1,]%>%na.omit()

s
b <- lapply(str_split(s$GeneRatio,"/"),as.numeric)
b <- sapply(b,function(x) temp=x[[1]]/x[[2]] )
s$GeneRatio <- b

ggplot(s,aes(GeneRatio,Description,size = Count))+
  geom_point(shape=21,aes(fill= pvalue),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,8))+
  theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=18,angle=0),
        axis.text.y  = element_text(color="black",size=18),
        legend.text = element_text(size =14,color="black"),
        legend.title = element_text(size =14,color="black"),
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))+scale_y_discrete(limits= unique(s$Description)) 



MtoH <- readRDS("MtoH.rds")
sym<- AnnotationDbi::select(org.Hs.eg.db,keys=list,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
strsplit2(s$geneID,"/")%>%c()
sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL
list_M <- MtoH[MtoH$HGNC.symbol%in%strsplit2(setdiff(GCOO_O$pairs,GCOO_Y$pairs),"_")[,1],]$MGI.symbol
list_M <- MtoH[MtoH$HGNC.symbol%in%sym[sym$ENTREZID%in%strsplit2(s$geneID,"/")%>%c(),]$SYMBOL,]$MGI.symbol

View(s)

GC1$Diets

data <- as.data.frame(GC1@assays$SCT@data)
data$SYMBOL <- rownames(data)
ggdata <- reshape2::melt(data)
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M6-BRD","M9-PRD","M9-BRD"))

head(ggdata)
data <- ggdata[ggdata$SYMBOL%in%c("Pdgfb","Col4a2","Ncam1","Fgfr2","Kitl","Efna4", "Tyro3","Axl","Nectin1","Fgfr1","Spp1","Vegfa","Efna1","Notch2"),]


DoHeatmap(GC1,features =list_M,group.by = "Diets",slot = "scale.data")+NoLegend()

DotPlot(GC1,features = c("Pdgfb","Col4a2","Ncam1","Fgfr2","Kitl","Efna4", "Tyro3","Axl","Nectin1","Fgfr1","Spp1","Vegfa","Efna1","Notch2"),group.by = "Diets")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,hjust = 1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(size=0.5),axis.line.y=element_line(size=0.5))

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












## iTALK
#devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)

# iTALK 要求的矩阵: 行为细胞，列为基因
#old <- c(rownames(OO[OO$avg_log2FC>0,]),rownames(GC[GC$avg_log2FC>0,]),
#  rownames(SC[SC$avg_log2FC>0,]),rownames(SMC[SMC$avg_log2FC>0,]),
#  rownames(M[M$avg_log2FC>0,]),rownames(NKT[NKT$avg_log2FC>0,]),
#  rownames(EC[EC$avg_log2FC>0,]))
count_O <- as.data.frame(ovary@assays$RNA@data)[old,meta_O$Cell]
iTalk_O <- as.data.frame(t(count_O))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_O$cell_type <- meta_O$cell_type
unique(iTalk_O$cell_type)
my10colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87')
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

NetView(iTalk_res,col=cell_col,vertex.label.cex=0.5,
        arrow.width=0.5,edge.max.width=2)

res <- iTalk_res[iTalk_res$cell_from%in%c("SC","GC","M","NKT","EC","SMC"),]
res <- res[res$cell_to%in%c("OO"),]

#res <- res[res$comm_type=="growth factor",]

italk_oo<- res

NetView(res,col=cell_col,vertex.label.cex=1,
        arrow.width=1,edge.max.width=2)

#View(res)

#res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]

res1 <- italk_oo[italk_oo$ligand%in%setdiff(italk_oo$ligand,italk_yy$ligand),]

LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,link.arr.width=res1$cell_to_mean_exprs)

res3 <- res2[order(res2$cell_from_mean_exprs*res2$cell_to_mean_exprs,decreasing=T),][1:20,]

LRPlot(res3,datatype='mean count',cell_col=cell_col,link.arr.lwd=res3$cell_from_mean_exprs,link.arr.width=res3$cell_to_mean_exprs)


NetView(res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),],col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:50,]
LRPlot(iTalk_res,datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs,link.arr.width=iTalk_res$cell_to_mean_exprs)
?NetView
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)


# iTALK 要求的矩阵: 行为细胞，列为基因
#young <- c(rownames(OO[OO$avg_log2FC<0,]),rownames(GC[GC$avg_log2FC<0,]),
#        rownames(SC[SC$avg_log2FC<0,]),rownames(SMC[SMC$avg_log2FC<0,]),
#        rownames(M[M$avg_log2FC<0,]),rownames(NKT[NKT$avg_log2FC<0,]),
#        rownames(EC[EC$avg_log2FC<0,]))
count_Y <- as.data.frame(ovary@assays$RNA@data)[young,meta_Y$Cell]
iTalk_Y <- as.data.frame(t(count_Y))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_Y$cell_type <- meta_Y$cell_type
unique(iTalk_Y$cell_type)
my10colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87')
highly_exprs_genes <- rawParse(iTalk_Y, top_genes=50, stats="mean")
# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- levels(iTalk_Y$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
names(cell_col) <- cell_types

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}

NetView(iTalk_res,col=cell_col,vertex.label.cex=0.5,
        arrow.width=0.5,edge.max.width=0.5)


res <- iTalk_res[iTalk_res$cell_from%in%c("SC","GC","M","NKT","EC","SMC"),]
res <- res[res$cell_to%in%c("OO"),]
italk_yy <- res
#res <- res[res$comm_type=="growth factor",]
write.csv(italk_yy,"italk_yy.csv")
write.csv(italk_oo,"italk_oo.csv")

NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=2)


#res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:30,]
res1 <- italk_yy[italk_yy$ligand%in%setdiff(italk_yy$ligand,italk_oo$ligand),]

#res1 <- res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
LRPlot(res1,datatype='mean count',cell_col=cell_col,link.arr.lwd=res1$cell_from_mean_exprs,link.arr.width=res1$cell_to_mean_exprs)
NetView(res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),],col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

unique(res$cell_to)


iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:20,]
iTalk_res
LRPlot(iTalk_res,datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs,link.arr.width=iTalk_res$cell_to_mean_exprs)


old_l <- italk_oo[italk_oo$ligand%in%setdiff(italk_oo$ligand,italk_yy$ligand),]$ligand
young_l <- italk_yy[italk_yy$ligand%in%setdiff(italk_yy$ligand,italk_oo$ligand),]$ligand
saveRDS(old_l,"./old_l.rds")
saveRDS(young_l,"./young_l.rds")
MtoH <- readRDS("./MtoH.rds")
old_l <- readRDS("./old_l.rds")
young_l <- readRDS("./young_l.rds")

a=MtoH[MtoH$HGNC.symbol%in%young_l,]$MGI.symbol
b=MtoH[MtoH$HGNC.symbol%in%old_l,]$MGI.symbol
BP <- readRDS("./BP.rds")


BP@active.ident

DotPlot(BP,features =a,cols = c("#584EA1","#ED1C26"),group.by = "Diets")+
  theme(axis.text.x = element_text(angle = 45,hjust =1))+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))+coord_flip()

 DotPlot(BP,features =c("Cyp19a1","Cyp17a1","Cyp11a1","Cyp11c1","Cyp21a1",
                        "Hsd3b1", "Hsd3b2","Hsd11b2","Hsd17b1", "Hsd17b3",
                        "Srd5a1","Srd5a2","Esr1","Esr2"),cols = c("#584EA1","#ED1C26"),group.by = "Diets")+
  theme(axis.text.x = element_text(angle = 45,hjust =1))+
  scale_y_discrete(limits=c("M6-PRD","M9-PRD","M9-BRD"))+coord_flip()

data <- as.data.frame(BP@assays$SCT@data)[
  rownames(BP)%in%a,]
ggdata <- reshape2::melt(data)
rm(BP,data)
gc()
head(ggdata)
ggdata$group <- paste0(substr(ggdata$variable,1,4),"RD")
ggdata$group <- factor(ggdata$group,levels = c("M6-PRD","M9-PRD","M9-BRD"))

head(ggdata)
data <- subset(ggdata,value>0.1)
ggplot(na.omit(data), aes(x = group, y =log(value+1),fill=group)) +
  stat_boxplot(geom = "errorbar", width=0.3,lwd=0.5)+
  geom_boxplot(lwd=0.5,fatten = 0.5,width=0.6,notch = T,outlier.colour = NA)+
  theme_classic()+
  scale_fill_manual(values=c("#FED439","#709AE1","#8A9197","#d4AF81","#FD7446","#D5E4A2"))+
  geom_signif(comparisons = list(c("M6-PRD","M9-PRD"),c("M9-BRD","M9-PRD")),
              map_signif_level = T,step_increase = 0.1,size = 0.5,textsize = 3,
              tip_length = 0,vjust = 0.1)+ylab("Ligands Expression Level")+xlab("")+	
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 16, color = 'black'))+	
  theme(axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,
                                                                       hjust = 1),	
        legend.text = element_text(size=12),
        axis.line.x=element_line(size=0.6),axis.line.y=element_line(size=0.6))+
  theme(legend.position = "top")+ggtitle("")+ NoLegend()
  theme(plot.title = element_text(size = 16),strip.text.x =element_text(size = 16))





library(org.Mm.eg.db)
MtoH <- readRDS("./MtoH.rds")
italk_oo$LRs <- paste0(italk_oo$ligand,"_",italk_oo$receptor)
italk_yy$LRs <- paste0(italk_yy$ligand,"_",italk_oo$receptor)

italk_oo_gc <- italk_oo[italk_oo$cell_from%in%c("GC"),]
italk_yy_gc <- italk_oo[italk_yy$cell_from%in%c("GC"),]

italk_oo_sc <- italk_oo[italk_oo$cell_from%in%c("SC"),]
italk_yy_sc <- italk_oo[italk_yy$cell_from%in%c("SC"),]
library(limma)
library(clusterProfiler)
a=setdiff(italk_oo_gc$ligand,italk_yy_gc$ligand)

b=setdiff(italk_yy_gc$ligand,italk_oo_gc$ligand)
a
b
c=setdiff(italk_oo$ligand,italk_yy$ligand)

d=setdiff(italk_yy$ligand,italk_oo$ligand)

e=setdiff(italk_oo_sc$ligand,italk_yy_sc$ligand)

f=setdiff(italk_yy_sc$ligand,italk_oo_sc$ligand)



sym<- AnnotationDbi::select(org.Mm.eg.db,keys=MtoH[MtoH$HGNC.symbol%in%a,]$MGI.symbol,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_oo_gc <- enrichKEGG(gene =sym$ENTREZID,
                          organism = 'mmu',
                          pvalueCutoff =0.05)
View(KEGG_oo_gc@result)
write.table(KEGG_oo_gc@result, 'KEGG_oo_gc.txt', sep='\t', quote=F)

sym<- AnnotationDbi::select(org.Mm.eg.db,keys=MtoH[MtoH$HGNC.symbol%in%b,]$MGI.symbol,
                            columns=c("ENSEMBL","SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_yy_gc <- enrichKEGG(gene =sym$ENTREZID,
                         organism = 'mmu',
                         pvalueCutoff =0.05)
View(KEGG_yy_gc@result)
write.table(KEGG_yy_gc@result, 'KEGG_yy_gc.txt', sep='\t', quote=F)

BP_oo_gc<- clusterProfiler::simplify(enrichGO(gene=a,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="BP",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))
BP_yy_gc<- clusterProfiler::simplify(enrichGO(gene=b,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="BP",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))

MF_oo_gc<- clusterProfiler::simplify(enrichGO(gene=a,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="MF",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))
MF_yy_gc<- clusterProfiler::simplify(enrichGO(gene=b,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="MF",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))
BP_oo<- clusterProfiler::simplify(enrichGO(gene=c,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="BP",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))
BP_yy<- clusterProfiler::simplify(enrichGO(gene=d,keyType = "SYMBOL",
                                              OrgDb= org.Hs.eg.db, ont="BP",
                                              pAdjustMethod="BH",pvalueCutoff=0.05,
                                              qvalueCutoff=0.05))
write.table(BP_yy_gc@result, 'BP_yy_gc.txt', sep='\t', quote=F)
write.table(BP_oo_gc@result, 'BP_oo_gc.txt', sep='\t', quote=F)
write.table(MF_yy_gc@result, 'MF_yy_gc.txt', sep='\t', quote=F)
write.table(MF_oo_gc@result, 'MF_oo_gc.txt', sep='\t', quote=F)
write.table(BP_yy@result, 'BP_yy.txt', sep='\t', quote=F)
write.table(BP_oo@result, 'BP_oo.txt', sep='\t', quote=F)

BP_oo_sc<- clusterProfiler::simplify(enrichGO(gene=e,keyType = "SYMBOL",
                                           OrgDb= org.Hs.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
BP_yy_sc<- clusterProfiler::simplify(enrichGO(gene=f,keyType = "SYMBOL",
                                           OrgDb= org.Hs.eg.db, ont="BP",
                                           pAdjustMethod="BH",pvalueCutoff=0.05,
                                           qvalueCutoff=0.05))
write.table(BP_yy_sc@result, 'BP_yy_sc.txt', sep='\t', quote=F)
write.table(BP_oo_sc@result, 'BP_oo_sc.txt', sep='\t', quote=F)



View(BP_yy_gc@result)
View(BP_oo_gc@result)

KEGG_oo_gc <- read.delim("./KEGG_oo_gc.txt", row.names=1)

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










## cell cycle
cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ovary1 <- RunPCA(ovary1, features = VariableFeatures(ovary1), 
                 ndims.print = 1:20, nfeatures.print = 10)
ovary1@
DimPlot(ovary1,reduction='tsne',group.by = )


## cdkn2a expression
data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("CDKN2A","CDKN1A","CDKN2D","GPNMB","ARF"),]%>%t()%>%as.data.frame()
cdkn2a$group <- paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group)
cdkn2a$group <- factor(cdkn2a$group,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- melt(cdkn2a)

ggplot(ggdata,aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+scale_fill_manual(values=c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                               "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"))+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),c("GC_Y","GC_O"),
                                 c('SC_Y',"SC_O"),c("SMC_Y","SMC_O"),
                                 c("EC_Y","EC_O"),c("M_Y","M_O"),c("NKT_Y","NKT_O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Relative expression")+ facet_wrap( ~ variable, ncol=2,scales = "free_y")+
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())

ggplot(ggdata[ggdata$variable=="GPNMB",],aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+scale_fill_manual(values=c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                                     "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"))+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),c("GC_Y","GC_O"),
                                 c('SC_Y',"SC_O"),c("SMC_Y","SMC_O"),
                                 c("EC_Y","EC_O"),c("M_Y","M_O"),c("NKT_Y","NKT_O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Relative expression")+ facet_wrap( ~ variable, ncol=3,scales = "free_y")+
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())





data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),]%>%t()%>%as.data.frame()
cdkn2a$group <- ovary@meta.data$group
cdkn2a$group <- factor(cdkn2a$group,levels = c("Y","O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- melt(cdkn2a)

ggplot(ggdata,aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+
  geom_signif(comparisons = list(c("Y","O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Relative expression")+scale_fill_manual(values=c("#E7AE8F","#BF5217"))+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())+ 
  facet_wrap( ~ variable, ncol=2,scales = "free_y")

data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),]%>%t()%>%as.data.frame()
cdkn2a$group <- paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group)
cdkn2a$group <- factor(cdkn2a$group,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- melt(cdkn2a)

ggplot(ggdata,aes(x = group,y = value,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3)+
  geom_boxplot(width=0.6)+theme_classic()+scale_fill_manual(values=c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                                     "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"))+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),c("GC_Y","GC_O"),
                                 c('SC_Y',"SC_O"),c("SMC_Y","SMC_O"),
                                 c("EC_Y","EC_O"),c("M_Y","M_O"),c("NKT_Y","NKT_O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Aging biomarkers relative expression")+
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())


ggplot(cdkn2a,aes(x = group,y = CDKN2A,fill=group))+
  stat_boxplot(geom = "errorbar", width=0.3)+
  geom_boxplot(width=0.6)+theme_classic()+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),
                                 c('SC_Y',"SC_O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24,color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("CDKN2A Relative expression")



data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("MYC","APOE","ND3","ENSMFAG00000011609"),]%>%t()%>%as.data.frame()
cdkn2a$group <- paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group)
cdkn2a$group <- factor(cdkn2a$group,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- reshape2::melt(cdkn2a)

ggplot(ggdata,aes(x = group,y = log2(value+1),fill=group))+
  geom_violin()+theme_classic()+
  scale_fill_manual(values=c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                                     "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"))+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),c("GC_Y","GC_O"),
                                 c('SC_Y',"SC_O"),c("SMC_Y","SMC_O"),
                                 c("EC_Y","EC_O"),c("M_Y","M_O"),c("NKT_Y","NKT_O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2,test = "t.test")+
  stat_summary(fun=mean, geom="point", size=2) +
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Relative expression")+
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())+ 
  facet_wrap( ~ variable, ncol=2,scales = "free_y")


data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("MYC","APOE","ND3","ENSMFAG00000011609"),]%>%t()%>%as.data.frame()
cdkn2a$group <- ovary@meta.data$group
cdkn2a$group <- factor(cdkn2a$group,levels = c("Y","O"))
cdkn2a$cell <- rownames(cdkn2a)
ggdata <- reshape2::melt(cdkn2a)
ggdata$variable <- factor(ggdata$variable,levels = c("MYC","APOE","ND3","ENSMFAG00000011609"))
ggplot(ggdata,aes(x = group,y = value,fill=group))+
  geom_violin()+theme_classic()+
  geom_signif(comparisons = list(c("Y","O")),
              map_signif_level = T,step_increase = 0,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+
  theme(axis.title =element_text(size = 20),axis.text =element_text(size = 20,color = 'black'))+	
  theme(axis.text = element_text(size = 20),axis.text.x = element_text(),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Expression Level")+scale_fill_manual(values=c("#E7AE8F","#BF5217"))+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())+ 
  facet_wrap( ~ variable, ncol=2,scales = "free_y")

## cdkn2a cell fraction
data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("MYC","APOE","ND3","ENSMFAG00000011609"),]%>%t()%>%as.data.frame()
cdkn2a$group <- paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group)
cdkn2a$group <- factor(cdkn2a$group,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
cdkn2a$cell <- rownames(cdkn2a)

ggdata <- reshape2::melt(cdkn2a)
ggdata

mean(ggdata$value)+sd(ggdata$value)

a=table(ggdata[ggdata$value>mean(ggdata$value)+sd(ggdata$value),]$group)%>%as.data.frame()
b=table(ggdata[ggdata$value<mean(ggdata$value)+sd(ggdata$value),]$group)%>%as.data.frame()
a
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


## cdkn2a cell fraction
data <- as.data.frame(ovary@assays$RNA@data)
cdkn2a <- data[rownames(data)%in%c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),]%>%t()%>%as.data.frame()
cdkn2a$group <- paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group)
cdkn2a$group <- factor(cdkn2a$group,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
cdkn2a$cell <- rownames(cdkn2a)

ggdata <- reshape2::melt(cdkn2a)
ggdata

mean(ggdata$value)+sd(ggdata$value)

a=table(ggdata[ggdata$value>mean(ggdata$value)+sd(ggdata$value),]$group)%>%as.data.frame()
b=table(ggdata[ggdata$value<mean(ggdata$value)+sd(ggdata$value),]$group)%>%as.data.frame()
a
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



base::load("~/scRNA_aging/20230405.RData")
ovary@meta.data$age <- 0
ovary@meta.data[ovary@meta.data$orig.ident=="YF1",]$age <- 5
ovary@meta.data[ovary@meta.data$orig.ident=="YF2",]$age <- 5
ovary@meta.data[ovary@meta.data$orig.ident=="YF3",]$age <- 5
ovary@meta.data[ovary@meta.data$orig.ident=="YF4",]$age <- 4
ovary@meta.data[ovary@meta.data$orig.ident=="OF1",]$age <- 18
ovary@meta.data[ovary@meta.data$orig.ident=="OF2",]$age <- 19
ovary@meta.data[ovary@meta.data$orig.ident=="OF3",]$age <- 19
ovary@meta.data[ovary@meta.data$orig.ident=="OF4",]$age <- 20

data <- as.data.frame(ovary@assays$RNA@data)
dd <- data[rownames(data)%in%c("MYC","APOE","ND3","ENSMFAG00000011609"),]
dd[1:4,1:4]
dd$gene <- rownames(dd)
ggdata <- reshape2::melt(dd)
ggdata$age <- ovary@meta.data$age
ggdata$celltype <- ovary1@meta.data$celltype
ggdata$g <- substr(ggdata$variable,1,3)
head(ggdata)
d1 <- ggdata[ggdata$g%in%c("YF1","YF2","YF3","YF4"),]
d1 <- d1[grep("_Y",d1$celltype),]

d2 <- ggdata[ggdata$g%in%c("OF1","OF2","OF3","OF4"),]
d2 <- d2[grep("_O",d2$celltype),]
d2
ggplot(data=rbind(d1,d2), aes(x=age, y=value))+
  geom_point(color="black",size=1)+						
  stat_smooth(method="lm",se=T,color="red3")+	#stat_smooth(method="loess",se=T)+						
  stat_cor(data=ggdata, method = "spearman")+theme_classic()+	
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24,color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(),
        legend.text = element_text(size=16),legend.position = "top",legend.title = element_blank(),
        axis.line.x=element_line(size=1),axis.line.y=element_line(size=1))+ 
  facet_wrap( ~ gene, ncol=2,scales = "free_y")+ylab("Relative expression")+ 
  theme(strip.text.x = element_text(size=14),strip.background = element_blank())










IM <- subset(x = ovary,idents=c("NKT","M"))
ElbowPlot(IM,ndims = 50)
IM <- RunUMAP(IM, dims = 1:10)
IM <- RunTSNE(IM, dims = 1:10)
DimPlot(IM, reduction = "tsne", pt.size = 1, label = T)
DimPlot(IM, reduction = "tsne",group.by = "group")
IM@meta.data$group
?FindAllMarkers
df1=as.data.frame(IM@active.ident)
IM@meta.data$celltype=paste0(df1$`IM@active.ident`,"_",IM$group)
idents <- paste0(df1$`IM@active.ident`,"_",IM$group)
names(idents) <- rownames(IM@meta.data)
IM@active.ident<-factor(idents)
IM.markers <- FindAllMarkers(object = IM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
M_OY <- FindMarkers(IM, ident.1 = "M_O", ident.2 = "M_Y", min.pct = 0.25)
NKT_OY <- FindMarkers(IM, ident.1 = "NKT_O", ident.2 = "NKT_Y", min.pct = 0.25)
M_OY <- M_OY[order(M_OY$avg_log2FC),]
NKT_OY <- NKT_OY[order(NKT_OY$avg_log2FC),]

dim(M_OY[M_OY$avg_log2FC>0,])
dim(NKT_OY[NKT_OY$avg_log2FC>0,])
MO=M_OY[M_OY$avg_log2FC>0,]
NKTO=NKT_OY[NKT_OY$avg_log2FC>0,]
top100_MO=M_OY[M_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
top100_NKTO=NKT_OY[NKT_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 

mean_IM<-AverageExpression(IM)
data<- mean_IM$RNA[c(rownames(M_OY),rownames(NKT_OY)),]
n=t(scale(t(data)))
pheatmap(na.omit(n),show_rownames = F,fontsize = 18,cluster_rows = F,cluster_cols = F,
         border_color = NA,legend = NULL)
col_fun = colorRamp2(c(-1, 0, 1), c("#4575B4","#FEEBA1","#D73027"))
Heatmap(na.omit(n)[,c(2,1,4,3)],col = col_fun,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = FALSE,# 在热图上边增加注释
        column_title = NULL,heatmap_legend_param = list(
          title = "scale",
          title_position = "leftcenter-rot") ) # 不需要列标题

data<- mean_IM$RNA[rownames(mean_IM$RNA)%in%c(rownames(M_OY),rownames(NKT_OY)),]
n=t(scale(t(data)))
n[n>2]=2
n[n<-2]=-2
dim(data)
set.seed(200)
p=pheatmap(na.omit(n),show_rownames = T,fontsize = 18,cluster_rows = T,cluster_cols = F,border_color = "white",kmeans_k = 6)
p
cluster5 <- p$kmeans$cluster[p$kmeans$cluster==5]%>%as.data.frame()
cluster5
BP5<- clusterProfiler::simplify(enrichGO(gene=rownames(cluster5),keyType = "SYMBOL",
                                        OrgDb= org.Hs.eg.db, ont="BP",
                                        pAdjustMethod="BH",pvalueCutoff=0.05,
                                        qvalueCutoff=0.05))
View(BP5@result)
BP <- BP5@result[BP5@result$Description%in%c("RNA splicing","aromatic compound catabolic process","cellular component disassembly",
                                "response to oxidative stress","circadian regulation of gene expression",
                                "rhythmic process","histone acetylation","histone modification",
                                "response to topologically incorrect protein","macroautophagy",
                                "RNA destabilization","response to endoplasmic reticulum stress",
                                "response to reactive oxygen species","intrinsic apoptotic signaling pathway",
                                "response to lipopolysaccharide","positive regulation of miRNA-mediated gene silencing"),]

BP5@result <- rbind(BP,BP4@result[BP4@result$Description=="response to estradiol",])
#x <- enrichplot::pairwise_termsim(BP5)
#?emapplot
#enrichplot::emapplot(x)
dotplot(BP5,label_format = 100,showCategory = 50)+theme_classic()+theme(axis.title  = element_text(size=18),axis.text = element_text(size=18))

cluster1 <- p$kmeans$cluster[p$kmeans$cluster==1]%>%as.data.frame()
cluster1
BP1<- simplify(enrichGO(gene=rownames(cluster1),keyType = "SYMBOL",
                        OrgDb= org.Hs.eg.db, ont="BP",
                        pAdjustMethod="BH",pvalueCutoff=0.05,
                        qvalueCutoff=0.05))
View(BP1@result)
cluster4 <- p$kmeans$cluster[p$kmeans$cluster==4]%>%as.data.frame()
cluster4
BP4<- simplify(enrichGO(gene=rownames(cluster4),keyType = "SYMBOL",
                        OrgDb= org.Hs.eg.db, ont="BP",
                        pAdjustMethod="BH",pvalueCutoff=0.05,
                        qvalueCutoff=0.05))
View(BP4@result)
BP4@result[BP4@result$Description=="response to estradiol",]
## aging

## 

BP4





View(M_OY)
dim(M_OY[M_OY$avg_log2FC>0.25,])
BP_MO<- simplify(enrichGO(gene=rownames(M_OY[M_OY$avg_log2FC>0,]),keyType = "SYMBOL",
                                         OrgDb= org.Hs.eg.db, ont="BP",
                                         pAdjustMethod="BH",pvalueCutoff=0.05,
                                         qvalueCutoff=0.05))
View(BP_MO@result)
BP_NKTO<- simplify(enrichGO(gene=rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),keyType = "SYMBOL",
                                        OrgDb= org.Hs.eg.db, ont="BP",
                                        pAdjustMethod="BH",pvalueCutoff=0.05,
                                        qvalueCutoff=0.05))
View(BP_NKTO@result)

write.table(IM.markers,"IM.markers.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(M_OY,"M_OY.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(NKT_OY,"NKT_OY.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(BP_NKTO,"BP_NKTO.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(BP_MO,"BP_MO.txt",col.names=T,row.names=T,quote=F,sep="\t")

BP_NKTO <- read.delim("~/scRNA_aging/BP_NKTO.txt", row.names=1)
BP_NKTO <- BP_NKTO[BP_NKTO$select==1,]%>%na.omit()
View(BP_NKTO)
BP_NKTO$group <- "NKT_O"

BP_MO <- read.delim("~/scRNA_aging/BP_MO.txt", row.names=1)
BP_MO <- BP_MO[BP_MO$select==1,]%>%na.omit()
View(BP_MO)
BP_MO$group <- "M_O"
intersect(BP_MO$Description,BP_NKTO$Description)
int <- intersect(rownames(M_OY),rownames(NKT_OY))
all <- c(rownames(NKT_OY[NKT_OY$avg_log2FC>0.25,]),rownames(M_OY[M_OY$avg_log2FC>0.25,]))
intersect(all,age_gene)
ggdata = list(age_gene,rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,]))
names(ggdata) <- c("ovary_aging","NKT_O","M_O")
ggdata
#BiocManager::install("ggvenn")
library(ggvenn)
p=ggvenn(ggdata,      
         show_percentage = F,show_elements=F,text_size = 10,
         stroke_color = "black",
         fill_color = c("#85AEC9","#A0B3A1","#B2662A"),
         set_name_color = c("black","black","black")) 
p

OO_295 <- setdiff(intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,])),age_gene)
BP_OO_295<- simplify(enrichGO(gene=OO_295,keyType = "SYMBOL",
                                          OrgDb= org.Hs.eg.db, ont="BP",
                                          pAdjustMethod="BH",pvalueCutoff=0.05,
                                          qvalueCutoff=0.05))

#View(BP_OO_295@result)
OO_322 <- intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,]))

BP_OO_322<- simplify(enrichGO(gene=OO_322,keyType = "SYMBOL",
                              OrgDb= org.Hs.eg.db, ont="BP",
                              pAdjustMethod="BH",pvalueCutoff=0.05,
                              qvalueCutoff=0.05))

View(BP_OO_322@result)



BP <- as.data.frame(BP_OO_295)
BP <- BP_OO_295[BP_OO_295$Description%in%c("regulation of translation","regulation of mitotic cell cycle","response to oxidative stress",
                               "response to reactive oxygen species","intrinsic apoptotic signaling pathway",
                               "circadian regulation of gene expression","RNA splicing","NIK/NF-kappaB signaling",
                               "negative regulation of transferase activity","cellular response to chemical stress","macroautophagy",
                               "interleukin-10 production","negative regulation of B cell activation",
                               "negative regulation of DNA-binding transcription factor activity","rhythmic process"),]

BP_OO_295@result <- BP
dotplot(BP_OO_295,label_format = 45)+theme_classic()+theme(axis.title  = element_text(size=18),axis.text = element_text(size=18))
#?BiocManager::install
#devtools::install_github("YuLab-SMU/GOSemSim")
#library(GOSemSim)
#devtools::install_github("YuLab-SMU/DOSE")
#devtools::install_github("YuLab-SMU/HDO.db")
#devtools::install_github('YuLab-SMU/clusterProfiler')
library(clusterProfiler)
ncbi<- AnnotationDbi::select(org.Hs.eg.db,keys=OO_322,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL")
KEGG_OO_322 <- enrichKEGG(gene =ncbi$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff =0.05)



library("GSVA")

kegggmt <- read.gmt("/home/tuyx/scRNA-Seq/new_10X/New_0.8/c2.cp.kegg.v7.3.symbols.gmt")
kegg_list = split(kegggmt$gene,gsub("KEGG_","",kegggmt$term))
gsva_IM <- GSVA::gsva(as.matrix(IM@assays$RNA@data),kegg_list,kcdf="Gaussian",parallel.sz=4)
data <- gsva_IM

library(dplyr)
meta <- as.data.frame(IM@meta.data[,c('orig.ident',"celltype")])
meta <- meta %>%
  arrange(meta$celltype)
dim(data)
meta$celltype%>%table(.)
data=as.data.frame(t(apply(data,1,function(a){
  tapply(a,c(rep("M_O",14),rep("M_Y",86),rep("NKT_O",108),rep("NKT_Y",101)),mean)})))
data
rownames(data) <- gsub("_"," ",rownames(data))%>%tolower(.)
write.table(data,"GSVA_IM.txt",col.names=T,row.names=T,quote=F,sep="\t")


M_OY <- FindMarkers(ovary1, ident.1 = "M_O", ident.2 = "M_Y", min.pct = 0.25)
NKT_OY <- FindMarkers(ovary1, ident.1 = "NKT_O", ident.2 = "NKT_Y", min.pct = 0.25)
OO_OY <- FindMarkers(ovary1, ident.1 = "OO_O", ident.2 = "OO_Y", min.pct = 0.25)
dim(OO_OY[OO_OY$avg_log2FC>0,])
GC_OY <- FindMarkers(ovary1, ident.1 = "GC_O", ident.2 = "GC_Y", min.pct = 0.25)
dim(GC_OY[GC_OY$avg_log2FC>0,])
SC_OY <- FindMarkers(ovary1, ident.1 = "SC_O", ident.2 = "SC_Y", min.pct = 0.25)
EC_OY <- FindMarkers(ovary1, ident.1 = "EC_O", ident.2 = "EC_Y", min.pct = 0.25)
SMC_OY <- FindMarkers(ovary1, ident.1 = "SMC_O", ident.2 = "SMC_Y", min.pct = 0.25)


library(dplyr)
OOO=OO_OY[OO_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC)
GCO=GC_OY[GC_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
SCO=SC_OY[SC_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
SMCO=SMC_OY[SMC_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
ECO=EC_OY[EC_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
NKTO=NKT_OY[NKT_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
MO=M_OY[M_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
a=intersect(rownames(OOO),rownames(GCO))%>%intersect(.,rownames(SCO))%>%
  intersect(.,rownames(SMCO))%>%intersect(.,rownames(ECO))%>%
  intersect(.,rownames(NKTO))%>%intersect(.,rownames(MO))
top_gene_OOO=OO_OY[OO_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC)%>%rownames()
top_gene_GCO=GC_OY[GC_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
top_gene_SCO=SC_OY[SC_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
top_gene_SMCO=SMC_OY[SMC_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
top_gene_ECO=EC_OY[EC_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
top_gene_NKTO=NKT_OY[NKT_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
top_gene_MO=M_OY[M_OY$avg_log2FC>0,]%>% top_n(n = 300, wt = avg_log2FC) %>%rownames()
c=intersect(top_gene_OOO,age_gene)
d=intersect(top_gene_GCO,age_gene)
e=intersect(top_gene_SCO,age_gene)
f=intersect(top_gene_SMCO,age_gene)
g=intersect(top_gene_ECO,age_gene)
h=intersect(top_gene_NKTO,age_gene)
i=intersect(top_gene_MO,age_gene)
a



VlnPlot(ovary1, features = c("MYC","APOE","ND3","ENSMFAG00000011609"),pt.size = 0,ncol = 2,cols =c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                          "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217") )+NoLegend()


VlnPlot(ovary1, features = c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),pt.size = 0,ncol = 2,cols =c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                                                                                   "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217") )+NoLegend()

plots_violins <- VlnPlot(ovary, features =c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),group.by = "group",
                         ncol = 2,cols =c("#E7AE8F","#BF5217"), 
                         log = FALSE,
                         combine = FALSE)+NoLegend()


VlnPlot(ovary1, features = c("MYC","APOE","ND3","ENSMFAG00000011609"),pt.size = 0,
        ncol = 2,cols =c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                         "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"),log = FALSE)+NoLegend()
plots_violins <- VlnPlot(ovary, 
                         cols = c("#E7AE8F","#BF5217"),
                         pt.size = 0,
                         group.by = "group",
                         features = c("CDKN2A","CDKN1A","CDKN2D","GPNMB"), 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
A <- singlecell_gene_test(ovary, 
                          genes.use = c("CDKN2A","CDKN1A","CDKN2D","GPNMB"),
                          group.by = 'group', 
                          comp = c("Y", "O"))
A
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
    scale_x_discrete(labels = c("Y","O"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

plots_violins[[1]]+plots_violins[[2]]+plots_violins[[3]]+plots_violins[[4]]


plots_violins <- VlnPlot(ovary1, 
                         cols = c("#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217",
                                  "#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217","#E7AE8F","#BF5217"),
                         pt.size = 0,
                         group.by = "celltype",
                         features = c("MYC","APOE","ND3","ENSMFAG00000011609"), 
                         ncol = 2, 
                         log = FALSE,
                         combine = FALSE)
A1 <- singlecell_gene_test(ovary1, 
                          genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                          group.by = 'celltype', 
                          comp = c("OO_Y","OO_O"))
A2 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("GC_Y","GC_O"))

A3 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("SC_Y","SC_O"))
A4 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("SMC_Y","SMC_O"))
A5 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("EC_Y","EC_O"))
A6 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("NKT_Y","NKT_O"))
A7 <- singlecell_gene_test(ovary1, 
                           genes.use = c("MYC","APOE","ND3","ENSMFAG00000011609"),
                           group.by = 'celltype', 
                           comp = c("M_Y","M_O"))

A <- rbind(A1,A2)%>%rbind(.,A3)%>%rbind(.,A4)%>%rbind(.,A5)%>%rbind(.,A6)%>%rbind(.,A7)
A
anno_pvalue <- format(A$p_val, scientific = T,digits = 3) 
anno_sig <- A$sig

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
    scale_x_discrete(labels = c("Y","O"))+
    geom_signif(annotations = anno_sig[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}

plots_violins[[1]]+plots_violins[[2]]+plots_violins[[3]]+plots_violins[[4]]








OOY=OO_OY[OO_OY$avg_log2FC<0,]
GCY=GC_OY[GC_OY$avg_log2FC<0,]
SCY=SC_OY[SC_OY$avg_log2FC<0,]
SMCY=SMC_OY[SMC_OY$avg_log2FC<0,]
ECY=EC_OY[EC_OY$avg_log2FC<0,]
NKTY=NKT_OY[NKT_OY$avg_log2FC<0,]
MY=M_OY[M_OY$avg_log2FC<0,]
b=intersect(rownames(OOY),rownames(GCY))%>%intersect(.,rownames(SCY))%>%
  intersect(.,rownames(SMCY))%>%intersect(.,rownames(ECY))%>%
  intersect(.,rownames(NKTY))%>%intersect(.,rownames(MY))
b
ovary1@active.ident <- factor(ovary1@active.ident,levels = c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O"))
DoHeatmap(ovary1, features =c(a,"CDKN2A","CDKN1A","CDKN2D","GPNMB",c,d,e,f,g,h,i))+NoLegend()+theme(axis.text = element_text(size = 14))

averageExp<-apply(ovary1@assays$RNA@data,1,function(x){tapply(x,paste0(ovary@meta.data$celltype,"_",ovary@meta.data$group),mean)}) %>% t() %>% as.data.frame()
# Expage <- averageExp$RNA[rownames(averageExp$RNA)%in%monkey_gene,]
Expage <- averageExp[rownames(averageExp)%in%c(a,"CDKN2A","CDKN1A","CDKN2D","GPNMB",c,d,e,f,g,h,i),]
n=t(scale(t(Expage)))
n[n>1]=1
n[n<-1]=-1
p <- pheatmap::pheatmap(na.omit(n)[c(a,"CDKN2A","CDKN1A","CDKN2D","GPNMB",c,d,e,f,g,i,h),c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O")],
                        show_rownames =T,fontsize =12,cluster_rows =F,cluster_cols = F)

col_fun = colorRamp2(c(-1, 0, 1), c("#4575B4","#FEEBA1","#D73027"))
Heatmap(na.omit(n)[unique(c(a,"CDKN2A","CDKN1A","CDKN2D","GPNMB",c,d,e,f,g,i,h)),c("OO_Y","OO_O","GC_Y","GC_O",'SC_Y',"SC_O","SMC_Y","SMC_O","EC_Y","EC_O","M_Y","M_O","NKT_Y","NKT_O")],
        col = col_fun,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names =T,# 在热图上边增加注释
        column_title = NULL,heatmap_legend_param = list(
          title = "scale",
          title_position = "leftcenter-rot") ) # 不需要列标题

########################################### select cell type ##########################################################
library(Seurat)
select <- subset(x = ovary,idents=c("OO","GC","NKT","M"))
#select <- ScaleData(select, features = intersect(intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,])),age_gene))
#select <- RunPCA(select, features = intersect(intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,])),age_gene))

#ElbowPlot(select,ndims = 50)
# Elbow plot and jackstraw plots are used to determine the number of PC to use 
#plan("multiprocess", workers = 4)
#select <- JackStraw(select, num.replicate = 100)
#select <- ScoreJackStraw(select, dims = 1:10)
#JackStrawPlot(select, dims = 1:10)
#select <- FindNeighbors(select, dims = 1:10)
#select <- FindClusters(select, resolution = 0.4)
select <- RunUMAP(select,features = intersect(intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,])),age_gene))
select <- RunTSNE(select,features = intersect(intersect(rownames(NKT_OY[NKT_OY$avg_log2FC>0,]),rownames(M_OY[M_OY$avg_log2FC>0,])),age_gene))
DimPlot(select, reduction = "tsne", pt.size = 1, label = T)
DimPlot(select, reduction = "tsne",group.by = "group")
# OO 5
FeaturePlot(object = select, features = c("SYCP3","DDX4","FIGLA"))
# GC 2 4
FeaturePlot(object = select, features = c("AMH","NR5A2"))

select@meta.data$group
?FindAllMarkers
df1=as.data.frame(select@active.ident)
select@meta.data$celltype=paste0(df1$`select@active.ident`,"_",select$group)
idents <- paste0(df1$`select@active.ident`,"_",select$group)
names(idents) <- rownames(select@meta.data)
select@active.ident<-factor(idents)
select.markers <- FindAllMarkers(object = select, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DoHeatmap(select, features =top10_select$gene) + NoLegend()


#M_OY <- FindMarkers(select, ident.1 = "M_O", ident.2 = "M_Y", min.pct = 0.25)
#NKT_OY <- FindMarkers(select, ident.1 = "NKT_O", ident.2 = "NKT_Y", min.pct = 0.25)
OO_OY <- FindMarkers(select, ident.1 = "OO_O", ident.2 = "OO_Y", min.pct = 0.25)
dim(OO_OY[OO_OY$avg_log2FC>0,])
GC_OY <- FindMarkers(select, ident.1 = "GC_O", ident.2 = "GC_Y", min.pct = 0.25)
dim(GC_OY[GC_OY$avg_log2FC>0,])
library(dplyr)
top100_OOO=OO_OY[OO_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
top100_GCO=GC_OY[GC_OY$avg_log2FC>0,]%>% top_n(n = 100, wt = avg_log2FC) 
OO_O <- OO_OY[OO_OY$avg_log2FC>0,]
GC_O <- GC_OY[GC_OY$avg_log2FC>0,]

ggdata = list(age_gene,rownames(OO_OY[OO_OY$avg_log2FC<0,]),rownames(GC_OY[GC_OY$avg_log2FC<0,]))
names(ggdata) <- c("ovary_aging","OO_O","GC_O")
ggdata
library(ggvenn)
p=ggvenn(ggdata,      
         show_percentage = F,show_elements=F,text_size = 10,
         stroke_color = "black",
         fill_color = c("#85AEC9","#E53F31","#FAC72E"),
         set_name_color = c("black","black","black")) 
p


data <- select@assays$RNA@data[rownames(select@assays$RNA@data)%in%age_gene,]%>%t()%>%as.data.frame()
dim(data)
data$celltype <- paste0(df1$`select@active.ident`,"_",select$group)
ggdata <- melt(data)

levels(data$celltype)
ggdata$celltype <- factor(ggdata$celltype,levels=c("OO_Y","OO_O","GC_Y","GC_O","M_Y","M_O","NKT_Y","NKT_O"))

ggplot(na.omit(ggdata),aes(x = celltype,y = value,fill=celltype))+
  stat_boxplot(geom = "errorbar", width=0.3)+
  geom_boxplot(width=0.6)+theme_classic()+
  scale_fill_manual(values=c("#0ECBC5","#0293A1","#FFDD93","#FFDD93","#E56A54","#F73E00","#6AA6D8","#8DB5BE"))+
  geom_signif(comparisons = list(c("OO_Y","OO_O"),c("GC_Y","GC_O"),
                                 c("M_Y", "M_O"),c("MKT_Y","NKT_O")),
              map_signif_level = T,step_increase = 0.1,tip_length = 0,vjust = 0.2)+
  stat_summary(fun=mean, geom="point", size=2) +xlab("")+guides(fill=FALSE)+
  theme(axis.title =element_text(size = 24),axis.text =element_text(size = 24,color = 'black'))+	
  theme(axis.text = element_text(size = 24),axis.text.x = element_text(angle = 45, hjust =1),	
        legend.text = element_text(size=12),legend.position = "right",legend.title = NULL,
        axis.line.x=element_line(),axis.line.y=element_line())+theme(legend.position = 'none')+
  ylab("Relative expression")


averageExp<-AverageExpression(select)
Expage <- averageExp$RNA[rownames(averageExp$RNA)%in%age_gene,]
n=t(scale(t(Expage)))
n[n>2]=2
n[n<-2]=-2
p <- pheatmap::pheatmap(na.omit(n),show_rownames = F,fontsize = 18,cluster_rows = T,cluster_cols = F)
?pheatmap

pheatmap::pheatmap(n[p$tree_row$order,c("OO_Y","OO_O","GC_Y","GC_O","M_Y","M_O","NKT_Y","NKT_O")],show_rownames = F,fontsize = 18,cluster_rows = F,cluster_cols = T)
select@active.ident
?DoHeatmap
DoHeatmap(select, features = age_gene) + NoLegend()
top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#0ECBC5","#0293A1","#FFDD93","#F5C31E","#E56A54","#F73E00","#6AA6D8","#8DB5BE")), # 设置填充色
                       labels = levels(factor(select$celltype,levels=c("OO_Y","OO_O","GC_Y","GC_O","M_Y","M_O","NKT_Y","NKT_O"))), 
                       labels_gp = gpar(cex = 0.5, col = "white",fontsize=24)))
col_fun = colorRamp2(c(-1, 0, 1), c("cornflowerblue","white","red"))

Heatmap(na.omit(as.data.frame(select@assays$RNA@data)),col = col_fun,
        cluster_rows = T,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = select$celltype,
        top_annotation = top_anno, # 在热图上边增加注释
        column_title = NULL,heatmap_legend_param = list(
          title = "scale",
          title_position = "leftcenter-rot") ) # 不需要列标题



####################### 添加显著性 ################################
singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "***")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "***")))
  }
  
  return(dfOUT)
}
