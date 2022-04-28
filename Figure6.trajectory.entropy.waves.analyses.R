library(entropy)
library(monocle3)

cds.p10  <- load_cellranger_data("../DATA/count_tommy10x_SM22ap10a")
cds.WK8<- load_cellranger_data("../DATA/count_tommy10x_SM22a8W1")


load("../resources/TF.CD.receptors.human.mouse.RData")
###################################################
load("../DATA/Tommy.NPP10.RData")
load("../DATA/Tommy.NPWK8.RData")

NPP10[["BioName"]]<-c("2"="1mNPC",  "1"="2progC1", "0"="3progC2", "3"="4CLC")[as.character(NPP10@meta.data$seurat_clusters)]
NPWK8[["BioName"]]<-c("0"="1mNPC1", "3"="2mNPC2","1"="4progC2", "2"="3progC1", "6"="5progC3", "4"="6cycling", "5"="7CLC")[as.character(NPWK8@meta.data$seurat_clusters)]

DimPlot(object = NPP10, label=T, 
		label.size=12,reduction= "tsne",group.by="NewName", pt.size=6)
DimPlot(object = NPP10, label=T, 
		label.size=12,reduction= "tsne",group.by="BioName", pt.size=6)

DimPlot(object = NPWK8, label=T, 
		label.size=12,reduction= "tsne",group.by="NewName", pt.size=6)
DimPlot(object = NPWK8, label=T, 
		label.size=12,reduction= "tsne",group.by="BioName", pt.size=6)

metaP10<-NPP10@meta.data
metaWK8<-NPWK8@meta.data

metaTwo<-rbind(metaP10[,c(1:3,5,8:9)], metaWK8[,c(1:3,5,9:10)])

###################################################

TWONP<-RunCCA(NPP10, NPWK8)

TWONP[["OldClust"]]<-paste0(TWONP@meta.data$orig.ident,"_",TWONP@meta.data$seurat_clusters)

TWONP<-NormalizeData(TWONP)
TWONP<-ScaleData(TWONP,fearures=rownames(TWONP))
TWONP<-FindVariableFeatures(TWONP)
TWONP<-RunPCA(TWONP)
TWONP<-RunUMAP(TWONP, reduction= "cca", dims= 1:10)
TWONP<-RunTSNE(TWONP, reduction= "cca", dims.use = 1:40, do.fast = T)
TWONP<-FindNeighbors(TWONP, reduction= "cca", dims = 1:10)
TWONP<-FindClusters(TWONP, resolution = 0.25)
###############################################
#       TRAJECTORY ANALYSES
###############################################
big_cds <- combine_cds(list(cdsNP.P10, cdsNP.WK8))

table(gsub(".*_","",colnames(TWONP))==gsub("_[12]","",colnames(big_cds)))

big_cds@colData[["orig.ident"]]<-metaTwo[["orig.ident"]]
big_cds@colData[["nCount_RNA"]]<-metaTwo[["nCount_RNA"]]
big_cds@colData[["nFeature_RNA"]]<-metaTwo[["nFeature_RNA"]]
big_cds@colData[["seurat_clusters"]]<-metaTwo[["seurat_clusters"]]
big_cds@colData[["NewName"]]<-metaTwo[["NewName"]]
big_cds@colData[["BioName"]]<-metaTwo[["BioName"]]
big_cds@colData[["BioName2"]]<-paste0(metaTwo[["orig.ident"]],"_",metaTwo[["BioName"]])


big_cds <- preprocess_cds(big_cds , num_dim = 20)
big_cds<- align_cds(big_cds, alignment_group = "orig.ident")


plot_pc_variance_explained(big_cds)
big_cds<- reduce_dimension(big_cds)


big_cds<- cluster_cells(big_cds)
big_cds<- learn_graph(big_cds, use_partition = F, close_loop = TRUE)

###############################################
#       ENTROPY ANALYSES
###############################################
plot_cells(big_cds, color_cells_by="BioName2", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)

indCells<-which(big_cds@int_colData$reducedDims$UMAP[,1]<6 & !big_cds@colData[["BioName2"]]%in%c("P10_4CLC", "WK8_7CLC"))
big_cds2<-big_cds[, indCells]

TWONP2<-subset(TWONP, cells =colnames(TWONP)[indCells])

ENTRO<-c()
for(i in 1:dim(TWONP2)[2]){
	if(i%%200==0)cat(i,"\t",date(),"\n")
	y<-as.matrix(GetAssayData(TWONP2)[,i])[,1]
	z<-log2(y[y>0])
	ENTRO[i]<-entropy(discretize(z,100))
}
table(gsub(".*_","",colnames(TWONP2))==gsub("_.*","",colnames(big_cds2)))
bymedian <- reorder(big_cds2@colData$BioName2, ENTRO, median)
boxplot(ENTRO~bymedian )

pdf("pseudotime.2022-apr.pdf")
	plot_cells(big_cds2, color_cells_by="BioName2", label_groups_by_cluster=FALSE,  cell_size =2, group_label_size = 6,  graph_label_size=5)
dev.off()

big_cds2<- order_cells(big_cds2)

bymedian <- reorder(big_cds2@colData$BioName2, pseudotime(big_cds2), median)
boxplot(pseudotime(big_cds2)~bymedian )

plot(pseudotime(big_cds2),ENTRO)

#########################################################
#
#           WAVES ANALYSES
#
#########################################################
library(ggplot2)
library(ggridges)
    
DAT.pseudo<-data.frame(bymedian, BioName2=big_cds2@colData$BioName2, pseudotime=pseudotime(big_cds2), entropy=ENTRO)
pdf("ridge.pseudo.pop.pdf", width=10)
	ggplot(DAT.pseudo, aes(x = pseudotime, y = bymedian)) +
		geom_density_ridges(scale = 4) + 
		#stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
		scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
		scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
		coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
		theme_ridges()

	bymedian3<-bymedian
	levels(bymedian3)<-rev(levels(bymedian))

	ggplot(DAT.pseudo, aes(x = entropy, y = bymedian)) +
		geom_density_ridges(scale = 4) + 
		#stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7) +
		scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
		scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
		coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
		theme_ridges()

dev.off()

#########################################################
NumWave<-7
NumPoint=1000

DNORM<-dnorm(seq(-6,6,len=(NumWave-1)*NumPoint*2+1))
plot(DNORM)


DROMmat<-matrix(NA, (NumWave-1)*NumPoint, NumWave)
for(i in 1:NumWave){
	Ni<- NumWave*NumPoint - i*NumPoint 
	DNORMi<-DNORM[-seq(Ni)][1:nrow(DROMmat)]
	DROMmat[,i]<-DNORMi
}
matplot(DROMmat, type="l", lty=1)



DATnorm<-data.frame(x=rep(seq(nrow(DROMmat)),NumWave),
	y=as.vector(DROMmat), 
	z=rep(seq(NumWave),rep(nrow(DROMmat),NumWave)))

ggplot(DATnorm, aes(x, y, group=z))+
		geom_point(aes(x,y)) +
		facet_wrap(vars(z), ncol=1)

pdf("wave.variables.pdf")

	matplot(DROMmat, type="l", lty=1)

dev.off()
####################################################

####################################################
PSEUDO<-round(pseudotime(big_cds2)/ max(pseudotime(big_cds2))*nrow(DROMmat))
range(PSEUDO)
hist(PSEUDO)

big_cds_500<-big_cds2[which(rowSums(exprs(big_cds2))>500),]
big_cds_TF<-big_cds2[which(rowData(big_cds2)$gene_short_name%in%TF.mouse & rowSums(exprs(big_cds2))>500),]

LCOEF_TF<-list()
LCOEF_ALL<-list()

for(i in 1:NumWave){
	cat(i, ":\t", date(), "\n")
	#if(i%%10==0)cat("\n")
	big_cds_TF@colData$DNORM<-DROMmat[,i][PSEUDO+1]
	big_cds_500@colData$DNORM<-DROMmat[,i][PSEUDO+1]

	TF_fits <- fit_models(big_cds_TF, model_formula_str = "~DNORM")
	TF_coefs <- coefficient_table(TF_fits)
	LCOEF_TF[[i]]<-TF_coefs

	ALL_fits <- fit_models(big_cds_500, model_formula_str = "~DNORM")
	ALL_coefs <- coefficient_table(ALL_fits)
	LCOEF_ALL[[i]]<-ALL_coefs
}

LHEAD_TF<-list()
LHEAD_ALL<-list()
for(i in 1:NumWave){
	TF_coefs<-LCOEF_TF[[i]]
	TF_coefs_top<-TF_coefs%>%filter(status=="OK" & q_value<0.05 & term=="DNORM"& (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_TF[[i]]<-TF_coefs_top$gene_short_name

	ALL_coefs<-LCOEF_ALL[[i]]
	ALL_coefs_top<-ALL_coefs%>%filter(status=="OK"& q_value<0.05 & term=="DNORM" & (estimate)>1) %>% arrange(desc(estimate))
	LHEAD_ALL[[i]]<-ALL_coefs_top$gene_short_name
}


####################################################
pdf("Five.waves.trendline.loess.pdf")
	for(i in 1:NumWave){
		#i<-6

		genei<-intersect(LHEAD_ALL[[5]], rownames(TWONP2))
		#genei<-c("Epas1")
		yy<-apply(as.matrix(GetAssayData(TWONP2)[genei,,drop=F]),2,mean, na.rm=T)
		pseu<-pseudotime(big_cds2)
		fitL<-loess(yy~pseu)

		if(0){
		for(j in 1:length(genei)){
			geneij<-genei[j]
			yyy<-as.vector(GetAssayData(TWONP)[geneij,])

			fitLj<-loess(yyy~pseu)
			predj<-predict(fitLj, newdata=sort(pseu), se=T)
			if(i==1){
				plot(sort(pseu), fitLj$fitted[order(pseu)])
			}else{
				lines(sort(pseu), fitLj$fitted[order(pseu)])
			}
		}}

		predi<-predict(fitL, newdata=sort(pseu), se=T)
		plot(pseu, yy)
		lines(sort(pseu), predi$fit, col="green")
		#lines(sort(pseu), predi$fit - qt(0.975,predi$df)* predi$se.fit)
		#lines(sort(pseu), predi$fit + qt(0.975,predi$df)* predi$se.fit)

	}
dev.off()

