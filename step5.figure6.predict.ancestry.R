ibrary(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gplots)
library(ggalluvial)

load("../DATA/Tommy.NPP10.RData")
load("../DATA/Tommy.NPWK8.RData")

NPP10[["BioName"]]<-c("2"="1NLC",  "1"="2progC1", "0"="3progC2", "3"="4CLC")[as.character(NPP10@meta.data$seurat_clusters)]
NPWK8[["BioName"]]<-c("0"="1NLC1", "3"="2NLC2","1"="4progC2", "2"="3progC1",
	"6"="5progC3", "4"="6cycling", "5"="7CLC")[as.character(NPWK8@meta.data$seurat_clusters)]

VlnPlot(NPP10, features="Tagln", pt.size=0, group.by="BioName")
VlnPlot(NPWK8, features="Tagln", pt.size=0, group.by="BioName")

TWONP@assays$RNA@var.features

set.seed(0)
indsub<-sample(ncol(NPWK8),ncol(NPP10))
subNPWK8<-subset(NPWK8, cells=colnames(NPWK8)[indsub])

set.seed(0)
TWONP<- RunCCA(NPP10, NPWK8,num.cc = 20)
TWONP<- NormalizeData(TWONP)
TWONP<-ScaleData(TWONP)
TWONP<-FindVariableFeatures(TWONP)
TWONP<- RunPCA(TWONP)
TWONP<- RunUMAP(TWONP, reduction= "cca", dims=1:20)
TWONP<- RunTSNE(TWONP, reduction= "pca", dims.use = 1:40, do.fast = T)
TWONP<- FindNeighbors(TWONP, dims = 1:10)
TWONP<- FindClusters(TWONP, resolution = 0.5)

UMAPPlot(TWONP, group.by="BioName", split.by="orig.ident", label=T, label.size=6, pt.size=6)

featuresTwoNPs <- intersect(TWONP@assays$RNA@var.features,intersect(rownames(NPP10),rownames(NPWK8)))

####################################



DIST1<-as.matrix(dist(TWONP@reductions$umap@cell.embeddings))
DIST2<-DIST1[seq(dim(NPP10)[2]), -seq(dim(NPP10)[2])]

predNeighbor<-apply(DIST2,2,function(x){
	BION10<-NPP10[["BioName"]][,1]
	ORD5<-order(x)[1:5]
	N5<-BION10[ORD5]
	tabIT<-sort(table(N5), decreasing=T)[1]
	names(tabIT)
})

table3<-table(NPWK8[["BioName"]][,1],predNeighbor)
pdf("alluvial.chart.with.pie.pdf")
	ggplot(as.data.frame(table3),
      	 aes(y = Freq, axis1 = predNeighbor, axis2 = Var1)) +
	  geom_alluvium(aes(fill = predNeighbor), width = 1/12) +
	  geom_stratum(width = 1/12, fill = "grey", color = "white") +
	  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
	  scale_fill_brewer(type = "qual", palette = "Set1") 

	for(i in 1:7)
		pie(table3[i,], main=i, col=hue_pal()(4))
dev.off()
