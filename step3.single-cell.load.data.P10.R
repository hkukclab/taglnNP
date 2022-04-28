library(Seurat)
library(dplyr)
library(patchwork)
load("../resources/TF.CD.receptors.human.mouse.RData")

PTEN.data <-Read10X(data.dir = "count_tommy10x_SM22a8W1/")
PTEN<-CreateSeuratObject(counts = PTEN.data, project = "WK8",min.cells=2,min.features=0)
PTEN<-RenameCells(PTEN, add.cell.id = "WK8")

PTEN<-NormalizeData(PTEN)
PTEN<-ScaleData(PTEN)
hist(PTEN$nFeature_RNA,breaks=1e2,xlab="Num of genes per cell",
		main="distributions of num-of-genes per cell in 8Wk")
abline(v=1350,col=2,lwd=3)

#############################################
PTEN<-FindVariableFeatures(PTEN)
topV10 <- VariableFeatures(PTEN)
topV10 <-topV10[!grepl("Gm|Rik|^Hb",topV10)]
plot1 <- VariableFeaturePlot(PTEN)
anyG<-setdiff(topV10,c(TF.mouse,surface.mouse))[1:50]
TFG<-intersect(topV10,TF.mouse)[1:20]
CDG<- intersect(topV10,surface.mouse)[1:20]
pdf("P10.top.variable.genes.pdf")
	LabelPoints(plot = plot1, points = c(anyG,TFG,CDG), repel = TRUE)
	LabelPoints(plot = plot1, points = anyG, repel = TRUE)
	LabelPoints(plot = plot1, points = TFG , repel = TRUE)
	LabelPoints(plot = plot1, points = CDG, repel = TRUE)
dev.off()

PTEN<-RunPCA(PTEN)
PTEN<-FindNeighbors(PTEN, dims = 1:10)
PTEN<-FindClusters(PTEN, resolution = 0.5)
PTEN<-RunUMAP(PTEN, reduction= "pca", dims= 1:40)
PTEN<-RunTSNE(PTEN, reduction= "pca", dims.use = 1:40, do.fast = T)


TSNEPlot(PTEN, pt.size = 2)

		"Lox","Alpl","Smpd3","Gdf5","Col10a1","Ihh","Epyc","Col9a1",
		"Dcn","Fmod","Col3a1","Prg4","Mkx","Scx","Omd","Thbs1","Cilp",
		"Aspn"), reduction="umap",ncol=6,
		cols=c("lightgrey","brown"))
#############################################
load("old-embed/LEMBED.RData")

PTEN2<-subset(PTEN, subset = nFeature_RNA>1350)
PTEN2<-RunPCA(PTEN2)
PTEN2<-FindNeighbors(PTEN2, dims = 1:10)
PTEN2<-FindClusters(PTEN2, resolution = 0.15)
PTEN2<-RunUMAP(PTEN2, reduction= "pca", dims= 1:40)
PTEN2<-RunTSNE(PTEN2, reduction= "pca", dims.use = 1:40, do.fast = T)

FeaturePlot(PTEN2,label=F, features=c("T","Hopx","Cd24a","Krt8","Krt18"),
		cols=c("lightgrey","brown"))
FeaturePlot(PTEN2,label=F, features=c("Clec3a","Clec3b"),
		cols=c("lightgrey","brown"),reduction="tsne",pt.size=2)

FeaturePlot(PTEN2,label=F, features=c("Col2a1","Col1a1","Bglap","Runx2","Spp1",
		"Tnmd","Omd","Cilp","Aspn","Thbs1","Prg4",
		"Sod3","Gdf5","Col9a1"),cols=c("lightgrey","brown"))


TSNE2<-PTEN2@reductions$tsne@cell.embeddings
ind12<-match(rownames(TSNE2),LEMBED[[5]]$Barcode)
MAT2<-as.matrix(LEMBED[[5]][ind12,2:3])
colnames(MAT2)<-colnames(TSNE2)
rownames(MAT2)<-rownames(TSNE2)
PTEN2@reductions$tsne@cell.embeddings<-MAT2
PL1<-TSNEPlot(PTEN2,label=T,label.size=12, pt.size = 2)+
	UMAPPlot(PTEN2,label=T,label.size=12, pt.size = 2)
#############################################
PL2<-DimPlot(object =PTEN2,reduction="tsne")
cells.located <- CellSelector(plot = PL2)
#############################################
PTEN3<-subset(PTEN2,idents=c(0:6,9))
PTEN3<-RunUMAP(PTEN3, reduction= "pca", dims= 1:40)
PTEN3<-FindNeighbors(PTEN3, dims = 1:10)
PTEN3<-FindClusters(PTEN3, resolution = 0.2)
TSNEPlot(PTEN3,label=T,label.size=12, pt.size = 2)+
	UMAPPlot(PTEN3,label=T,label.size=12, pt.size = 2)
pdf("P10.TSNE.wo.HSC.LYMPH.pdf")
	DimPlot(object =PTEN3,reduction="tsne",label=T,label.size=12, pt.size = 2)
	DimPlot(object =PTEN3,reduction="umap",label=T,label.size=12, pt.size = 2)
dev.off()

PTEN.markers <- FindAllMarkers(PTEN2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PTEN.markers3 <- FindAllMarkers(PTEN3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

LSIG<-split(as.character(PTEN.markers3$gene),PTEN.markers3$cluster)

LSIG2<-sapply(levels(PTEN.markers3$cluster),function(i){
	thisDEG<-PTEN.markers3[PTEN.markers3$cluster==i,]
	thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25& thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(thisDEG$gene)
})

TOP10<-unlist(sapply(LSIG2,function(x)head(x[!grepl("Gm|Rik|^a$",x)],n=10)))
pdf("heatmap.p10.pdf",height=10)
	DoHeatmap(PTEN3,  features = TOP10)
dev.off()

save(PTEN3,file="NP.only.PTEN3.RData")
