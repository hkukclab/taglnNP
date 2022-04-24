library(DESeq)

DESEQ<-function(ind.control,ind.case,pref1,pref2){
	lab<-as.factor(rep(c(pref1,pref2),c(length(ind.control),length(ind.case))))
	print(lab)
	datin<-dat.counts[,c(ind.control,ind.case)]


	ind.NZ<-which(rowSums(datin>0)>1)
	print(length(ind.NZ))
	countTable<-datin[ind.NZ,]

	GTF.noz<-annot.gtf[ind.NZ,]
	print(dim(countTable))

	de.experiment <- newCountDataSet( countTable, lab)

	sf.est <-estimateSizeFactors( de.experiment)
	sizes<-sizeFactors( sf.est)
	print(sizes)
	sf.est.disp<-  estimateDispersions( sf.est, fitType ="local" )

	str( fitInfo(sf.est.disp) )

	plotDispEsts(sf.est.disp)

	print(head( fData(sf.est.disp) ))
	print(c(2,4)[lab])

	res.test <- nbinomTest( sf.est.disp, pref1, pref2)
	res<-list(ind.NZ=ind.NZ,lab=lab,
		countTable=countTable,
		GTF.noz=GTF.noz,
		res.test=res.test,
		sizes=sizes) 
	res
}
library(scales)
VOCALNO_deseq<-function(res.deseq,title="",CORNER="top",onlyshow="",cexin=1/3,GENECOL=5){
	annoti<-res.deseq$GTF.noz
	testmat<-res.deseq$res.test

	LOGFC<-testmat$log2FoldChange
	LOGFC[is.infinite(LOGFC)]<-log2((testmat$baseMeanB+1)/(testmat$baseMeanA+1))[is.infinite(LOGFC)]

	PADJ<-testmat$padj
	GENES<-as.character(annoti[,GENECOL])
	indRBC<-grepl("^Hb[ab]",GENES)

	print(head(GENES))
	plot(LOGFC[!indRBC],-log10(testmat$pval[!indRBC]),
		xlab="log2(FC)",ylab="-log10(P)",pch=".",cex=2,
		col="grey",main=title)
	abline(v=0,lty=2)
	indUp<-which(PADJ<0.05&LOGFC>0&!indRBC)
	indDown<-which(PADJ<0.05&LOGFC< 0&!indRBC)

	if(onlyshow[1]==""){
		ind.show<-which(GENES%in%onlyshow)
		indUpShow<-indUp
		indDownShow<-indDown
		print(table(PADJ<0.05))

	}else{
		indUpShow<-which(PADJ<0.05&LOGFC>2&GENES%in%onlyshow)
		indDownShow<-which(PADJ<0.05&LOGFC< -2&GENES%in%onlyshow)
		print(indUpShow)
		print(indDownShow)
	}

	if(length(indUpShow)>0){
		points(LOGFC[indUpShow],-log10(testmat$pval)[indUpShow],pch=16,col="#F8766D")
		text(LOGFC[indUpShow],-log10(testmat$pval)[indUpShow],GENES[indUpShow],
			pos=rep(c(2,4),ceiling(length(indUpShow)/2)),cex=cexin,xpd=T)
	}
	if(length(indDownShow)>0){
		points(LOGFC[indDownShow],-log10(testmat$pval)[indDownShow],pch=16,col="#00BFC4")
		text(LOGFC[indDownShow],-log10(testmat$pval)[indDownShow],GENES[indDownShow],
			pos=rep(c(2,4),ceiling(length(indDownShow)/2)),cex=cexin,xpd=T)
	}
	print(levels(res.deseq$lab))
	str1<-paste0("Up in ",levels(res.deseq$lab),
		" (n=",c(length(indDown),length(indUp)),")")
	print(str1)
	legend(CORNER,
		pch=16,
		legend=str1,
		col=c("#00BFC4","#F8766D"))
	L1<-list(up=GENES[indUpShow],
		down=GENES[indDownShow])
	return(L1)
}

##############################################
##############################################
load("../htseq-count.gtf.outputs-grch38-bulk.RData")


colnames(dat.counts)<-gsub("_1.fast.*","",basename(list1))
dat.counts<-dat.counts[,-4]

shortIDs<-colnames(dat.counts)

GROUP<-gsub("_.*","",shortIDs)


################################
library(gplots)
################################

indCDS<-which(GROUP=="CDS")
indDS<-which(GROUP=="DS")

RES.DESeq<-DESEQ(indCDS,indDS,
			paste0("all_CDS"),
			paste0("all_DS"))

fnout<-"step5.DEseq2.CDS-DS.no147.May14.RData"
save(RES.DESeq,file=fnout)

##############################################
pdf("VOLCANO.DEseq2.CDS-DS.May14.pdf",width=12,height=6)

		layout(t(matrix(seq(2))),widths=c(6,6))
		res<-VOCALNO_deseq(RES.DESeq,"DS-vs-CDS",
			CORNER="topleft",cexin=1/2,GENECOL=4)
		##################
		source("G:/my-func-lib/CAT_DEGs.R")
		CAT_DEG2(res$up,res$down,BSIZE=3,WIDTH=90,Xpos=0.72)
		save(res,file="DEseq2.CDS-DS.no147.May14.RData")

dev.off()


library(ggplot2)
x<-dat.counts[which(annot.gtf[,4]=="TAGLN"),]

dat34<-data.frame(Taln=log2(x),GROUP=GROUP)

pdf("TAGLN.violin.pdf")
	p34<- ggplot(dat34, aes(x=GROUP, y=Taln)) + 
		geom_violin( scale = "width")+ 
		geom_boxplot(width=0.4,outlier.shape =NA)+ 
		geom_jitter(shape=16,size=4, position=position_jitter(0.2))+
		ggtitle("Tagln")
	plot(p34)

dev.off()


