#!/usr/bin/env Rscript

library('getopt');
#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
		'help' , 'h', 0, "logical",
		"outpath",'o',1,"character",
		"key",'k',1,"character"
		), byrow=TRUE, ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
		Usage example: 
		1) Rscript topGO.R --outpath /out/path --key Citrullus_lanatus_97103v1

		Options: 
		--help          -h  NULL        get this help
		--outpath       -o  character   the path for output file [optional,current working directory]
		--key           -k  character   the key words for output file [optional,default:'demo'] 
		\n")
	q(status=1);
}
library(topGO)
library(ALL)
data(ALL)

geneID2GO<-readMappings(paste(opt$outpath,"topGO.map",sep="/"))
#data<-read.table("opt$outpath/topGO.list", row.names = 1, header=TRUE)
data<-read.table(paste(opt$outpath,"topGO.list",sep="/"), row.names = 1, header=TRUE,check.names =F)
geneList<-data[,1]
names(geneList) <- rownames(data)
topDiffGenes<-function(allScore){return(allScore<0.05)} #返回值是TRUE or FALSE

go<-c('BP','MF','CC')
for(i in go)
{
	sampleGOdata <- new("topGOdata",nodeSize = 6,ontology=i, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO,geneSel=topDiffGenes)
#	resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
#	allRes <- GenTable(sampleGOdata,KS = resultKS.elim,ranksOf = "classic", topNodes = attributes(resultKS.elim)$geneData[4], numChar=200)
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
	allRes <- GenTable(sampleGOdata,Pvalue = resultFisher,ranksOf = "classicFisher", topNodes = attributes(resultFisher)$geneData[4], numChar=200)

	#write.table(allRes, file="opt$outpath/opt$key.topGO_BP.xls", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	write.table(allRes, paste(opt$outpath,"/",opt$key,".topGO_",i,'.xls',sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	################################# add gene list to topGO.*,xls
	sig_genes<-sapply(1:nrow(allRes),function(x){intersect(sigGenes(sampleGOdata),genesInTerm(sampleGOdata,allRes[x,1])[[1]])})
	allRes<-cbind(allRes,sapply(1:length(sig_genes),function(x){paste(sig_genes[[x]],collapse=',')}))
	colnames(allRes)[7]<-'Significant.genes'
	write.table(allRes, file=paste(opt$outpath,"/",opt$key,".topGO_",i,"_gene.xls",sep=""), sep="\t", quote=FALSE,col.names=TRUE, row.names=FALSE)

	###  score(resultKS.elim) can not be 0. (cause ERROR)
#	pvalue <- score(resultKS.elim)
    pvalue <- score(resultFisher)
	pvalue[pvalue==0]= 1e-30

	#pdf("opt$outpath/opt$key.topGO_BP.pdf")
	pdf(paste(opt$outpath,"/",opt$key,".topGO_",i,".pdf",sep=""))
	showSigOfNodes(sampleGOdata, pvalue, firstSigNodes = 10, useInfo = "all")
	dev.off()
	#png("opt$outpath/opt$key.topGO_BP.png")
	png(paste(opt$outpath,"/",opt$key,".topGO_",i,".png",sep=""))
	showSigOfNodes(sampleGOdata, pvalue, firstSigNodes = 10, useInfo = "all")
	dev.off()

}
