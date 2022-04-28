library(Seurat)
library(dplyr)
library(patchwork)
load("../resources/TF.CD.receptors.human.mouse.RData")


EightWk.data <-Read10X(data.dir = "count_tommy10x_SM22a8W1/")
EightWk<-CreateSeuratObject(counts = EightWk.data, project = "WK8",min.cells=2,min.features=0)
EightWk<-RenameCells(EightWk, add.cell.id = "WK8")

EightWk<-NormalizeData(EightWk)
EightWk<-ScaleData(EightWk)
hist(EightWk$nFeature_RNA,breaks=1e2,xlab="Num of genes per cell",
		main="distributions of num-of-genes per cell in 8Wk")
abline(v=1350,col=2,lwd=3)

#############################################
EightWk<-FindVariableFeatures(EightWk)
topV10 <- VariableFeatures(EightWk)
topV10 <-topV10[!grepl("Gm|Rik|^Hb",topV10)]
plot1 <- VariableFeaturePlot(EightWk)
anyG<-setdiff(topV10,c(TF.mouse,surface.mouse))[1:50]
TFG<-intersect(topV10,TF.mouse)[1:20]
CDG<- intersect(topV10,surface.mouse)[1:20]
pdf("8wk.top.variable.genes.pdf")
	LabelPoints(plot = plot1, points = c(anyG,TFG,CDG), repel = TRUE)
	LabelPoints(plot = plot1, points = anyG, repel = TRUE)
	LabelPoints(plot = plot1, points = TFG , repel = TRUE)
	LabelPoints(plot = plot1, points = CDG, repel = TRUE)
dev.off()

EightWk<-RunPCA(EightWk)
EightWk<-FindNeighbors(EightWk, dims = 1:10)
EightWk<-FindClusters(EightWk, resolution = 0.5)
EightWk<-RunUMAP(EightWk, reduction= "pca", dims= 1:40)
EightWk<-RunTSNE(EightWk, reduction= "pca", dims.use = 1:40, do.fast = T)


TSNEPlot(EightWk, pt.size = 2)

#############################################
load("old-embed/LEMBED.RData")

EightWk2<-subset(EightWk, subset = nFeature_RNA>1350)
EightWk2<-RunPCA(EightWk2)
EightWk2<-FindNeighbors(EightWk2, dims = 1:10)
EightWk2<-FindClusters(EightWk2, resolution = 0.15)
EightWk2<-RunUMAP(EightWk2, reduction= "pca", dims= 1:40)
EightWk2<-RunTSNE(EightWk2, reduction= "pca", dims.use = 1:40, do.fast = T)


TSNE2<-EightWk2@reductions$tsne@cell.embeddings
ind12<-match(rownames(TSNE2),LEMBED[[2]]$Barcode)
MAT2<-as.matrix(LEMBED[[2]][ind12,2:3])
colnames(MAT2)<-colnames(TSNE2)
rownames(MAT2)<-rownames(TSNE2)
EightWk2@reductions$tsne@cell.embeddings<-MAT2


TSNEPlot(EightWk2,label=T,label.size=12, pt.size = 2)
FeaturePlot(EightWk2,label=F, features=c("Clec3a","Clec3b"),
		cols=c("lightgrey","brown"),reduction="tsne",pt.size=2)

UMAPPlot(EightWk2,label=T,label.size=12, pt.size = 2)
#############################################
EightWk3<-subset(EightWk2,idents=c(0:4,6))
cells.located <- CellSelector(plot = DimPlot(EightWk3,reduction="tsne"))
EightWk4<-subset(EightWk3,cells=setdiff(colnames(EightWk2),"WK8_GCATCGGCAGTCGGTC-1"))
EightWk4<-FindNeighbors(EightWk4, dims = 1:10)
EightWk4<-FindClusters(EightWk4, resolution = 0.2)

TSNEPlot(EightWk3,label=T,label.size=12, pt.size = 2)

EiWk.markers <- FindAllMarkers(EightWk3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EiWk.markers4 <- FindAllMarkers(EightWk4, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

save(EightWk4,file="NP.only.EightWk4.RData")
