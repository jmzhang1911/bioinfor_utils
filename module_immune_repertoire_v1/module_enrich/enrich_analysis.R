#!/usr/bin/env Rscript
library(getopt)
# parameter setting
spec <- matrix(c(
  'help'        , 'h', 0, "logical"  , "print help usage",
  'go'          , 'g', 1, "character", "character	GO info contained three column, eg:go_id	GO_term name	gene_id",
  'kegg'        , 'k', 1, "character", "character	KEGG info contained three column, eg:ko_id	pathway_name	gene_id",
  'deg'         , 'd', 1, "character", "character	input , cluster diff marker gene file",
  'od'          , 'o', 1, "character", "enrich.out dir, default current path",
  'prefix'      , 'p', 1, "character", "character	enrich.output key,default cluster",
  'enrichn'     , 'n', 1, "integer"  , "int the first enrichn term would be plotted, default 10",
  "color"       , 'c', 1, "character", "character	One of pvalue p.adjust qvalue, default p.adjust",
  "column"      , 'C', 1, "integer"  , "interger	for deg file which column is gene id, default 1",
  "len"         , "l", 1, "integer"  , "int	The length of GO/KEGG term would be Truncated. default 50",
  'label.len'   , 'N', 1, "integer"  , "label length, default:100",
  'all'         , 'A', 0, 'logical'  , "-all, show full go term name",
  'x.lab'       , 'X', 1, "character", "dot plot x lab",
  'y.lab'       , 'Y', 1, "character", "dot plot y lab",
  'title.lab'   , 'T', 1, "character", "dot plot title",
  'lab.size'    , 'L', 1, "integer"  , "lab size",
  'axis.size'   , 'a', 1, "integer"  , "axis size",
  'legend.size' , 'D', 1, "integer"  , "legend size"
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:enrich_analysis.R
  Version: Version v1.0
  Description: this program is used to do enrich analysis ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help))	 PrintUsage()
if(is.null(opt$deg))	   PrintUsage()
if(is.null(opt$go) && is.null(opt$kegg)) PrintUsage()
if(is.null(opt$od))	     { opt$od <- "./" }
if(is.null(opt$prefix))  { opt$prefix <- "cluster" }
if(is.null(opt$color))	 { opt$color <- "qvalue" }
if(is.null(opt$column))	 { opt$column <- 1 }
if(is.null(opt$enrichn)) { opt$enrichn <- 20 }
if(is.null(opt$label.len)) { opt$label.len <- 100 }
#if(is.null(opt$len))	   opt$len <- 50
if ( is.null(opt$y.lab) )	{ opt$y.lab <- "" }
if ( is.null(opt$x.lab) )	{ opt$x.lab <- 'Enrichment Factor' }
if ( is.null(opt$height ) )		{ opt$height <- 3000 }
if ( is.null(opt$width ) )		{ opt$width <- 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size <- 14 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size <- 9 }
if ( is.null(opt$legend.size ) )	{ opt$legend.size <- 12 }
if(!dir.exists(opt$od))  {dir.create(opt$od)}

# load packages
library(farver)
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(ggplot2))
library(RColorBrewer)

# read data
GO <- read.csv(opt$go, sep="\t", header=F)
KEGG <- read.csv(opt$kegg, sep="\t", header=F)
marker.gene <- read.csv(opt$deg, sep="\t", header=T)
genes <- unique(sort(as.character(marker.gene[,opt$column])))
if(length(genes) < 1){
  print("No marker gene exists!")
  q()
}
my.width <- 6
my.height <- 4
# save plot
SavePlot <- function(od, filename, data, width = my.width, height = my.height, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}
SavePlot1 <- function(od, filename, data, width = my.width, height = my.height, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width-1, height = height-1, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width-1, height = height-1, scale = scale)
}


# barplot
PrivateBarplot <- function(enrich.data, top){
	bar.data <- enrich.data %>% as_tibble()
	bar.data <- bar.data[match(unique(bar.data$Description), bar.data$Description),]
	plot.data <- bar.data %>% group_by(Ontology) %>% arrange(qvalue) %>% slice(1:top)
	if(!is.null(opt$all)){
          Description <- plot.data %>% .[, 3, drop = TRUE]
        }else{
          Description <- plot.data %>% .[, 3, drop = TRUE] %>% substring(., 1, opt$len) %>% paste(., "...", sep = "")
        }
	plot.data$Description <- Description
	term <- factor(plot.data$Description, levels = plot.data$Description)
	log10Pvalue <- -(log(plot.data$pvalue)/log(10))
	bar.data <- data.frame(Description = term, log10Pvalue = log10Pvalue, Ontology = plot.data$Ontology, gene_number = plot.data$gene_number)
	p <- ggplot(bar.data , aes(Description, log10Pvalue, fill = Ontology)) + geom_col(position = 'dodge', width = 0.5) + scale_x_discrete(labels = function(x) str_wrap(x, width = opt$label.len)) + geom_text(aes(label = gene_number), size = 2, color = "gray20", position = position_dodge(0.9), hjust = 0.5, vjust = -0.5) 
	p1 <- p + scale_fill_manual(values = brewer.pal(8, "Set2")[1:3])
	# p2 <- p1 + theme(axis.text.x = element_text(size = 6, vjust = 1, hjust = 1, angle = 45), axis.text.y = element_text(size = 4), axis.title.x = element_blank()) + ylab("-log10pvalue") 
	p2 <- p1 + theme(axis.text.x = element_text(family="serif", size = opt$axis.size-2, angle = 315, hjust=0,vjust=1), axis.title.x = element_blank(), axis.text.y = element_text(family="serif", size = opt$axis.size), legend.title = element_text(family = "serif"), legend.text = element_text(family = "serif"), axis.title.y = element_text(family = "serif", size = opt$legend.size)) + ylab("-log10pvalue")
	file.name <- paste(opt$prefix, "go_enrich_barplot", sep = "_")
	SavePlot(od = opt$od, filename = file.name, data = p2, width = 8)
}


# enrich analysis
EnrichAnalysis<-function(genes, type, info){
  type <- chartr(" ","_",type)
  
  TERM2GENE <- data.frame(info[,1], info[,3])
  TERM2NAME <- data.frame(info[,1], info[,2])
  TERM2NAME <- unique(TERM2NAME[order(TERM2NAME[,1]),])
  enrich <- enricher(genes,TERM2GENE=TERM2GENE,TERM2NAME=TERM2NAME,pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "fdr")
  if(is.null(enrich)){
    print("null")
    return(1)
  }
  enrich <- data.frame(enrich)
  enrich_factor <- enrich %>% select(3,4) %>% as.matrix %>% 
    apply(., 1, function(x){
      GeneRatio <- strsplit(x[1], split = "/") %>% unlist %>% as.numeric 
      BgRatio <- strsplit(x[2], split = "/") %>% unlist %>% as.numeric
      rich <- as.numeric((GeneRatio[1]/GeneRatio[2])/(BgRatio[1]/BgRatio[2]))
      return(rich)}) %>% round(., 2)
  gene.number <- enrich %>% select(8) %>% .[, 1, drop = TRUE] %>% as.character %>% 
    sapply(., function(x){num <- strsplit(x, split = "/") %>% unlist() %>% length ; return(num)}) %>% unlist
  enrich.out <- data.frame(enrich$ID, enrich$Description, enrich$GeneRatio, enrich$BgRatio, 
                           enrich_factor, enrich$pvalue, enrich$qvalue, enrich$geneID, gene_number = gene.number)
  enrich.outfile <- paste(opt$od, "/", opt$prefix, "_", type, "_enrich.list", sep="")
  header <- c("ID", "Description", "GeneRatio", "BgRatio", "enrich_factor", "pvalue", "qvalue", "geneID", "gene_number")
  write.table(t(header), enrich.outfile, row.names=F, col.names=F, quote=F, sep="\t")
  enrich.out <- enrich.out[order(enrich.out[,7]),]
  write.table(enrich.out, enrich.outfile, row.names=F, col.names=F, quote=F,sep="\t", append=T)
  
  if(type != "KEGG_pathway"){
    Ontology <- rep(type, nrow(enrich.out))
    enrich.data <- data.frame(Ontology = Ontology, enrich.out) 
	write.table(enrich.data, all.outfile, row.names=F, col.names=F, quote=F, sep="\t", append=T)
  }
  TERM2NAME <- data.frame(info[,1], info[,2])
  TERM2NAME <- unique(TERM2NAME[order(TERM2NAME[,1]),])
  enrich.new <- enricher(genes, TERM2GENE=TERM2GENE, TERM2NAME=TERM2NAME, pvalueCutoff=1, qvalueCutoff=1, pAdjustMethod = "fdr")
  if(is.null(opt$all)){
   enrich.new@result$Description <- substring(enrich.new@result$Description, 1, opt$len) %>% paste(., "...", sep = "")
  }
  # enrich <- data.frame(enrich.new)
  enrich <- enrich.new@result
  enrich_factor <- enrich %>% select(3,4) %>% as.matrix %>% 
    apply(., 1, function(x){
      GeneRatio <- strsplit(x[1], split = "/") %>% unlist %>% as.numeric 
      BgRatio <- strsplit(x[2], split = "/") %>% unlist %>% as.numeric
      rich <- as.numeric((GeneRatio[1]/GeneRatio[2])/(BgRatio[1]/BgRatio[2]))
      return(rich)}) %>% round(., 2)
  gene.number <- enrich %>% select(8) %>% .[, 1, drop = TRUE] %>% as.character %>% 
    sapply(., function(x){num <- strsplit(x, split = "/") %>% unlist() %>% length ; return(num)}) %>% unlist
  enrich.plot <- data.frame(enrich$ID, enrich$Description, enrich$GeneRatio, enrich$BgRatio, 
                           enrich_factor, enrich$pvalue, enrich$qvalue, enrich$geneID, gene_number = gene.number)
  names(enrich.plot) <- c("ID", "Description", "GeneRatio", "BgRatio", "enrich_factor", "pvalue", "qvalue", "geneID", "gene_number")
  if(nrow(enrich.plot) < opt$enrichn){
	enrich.plot <- enrich.plot[order(enrich.plot[,7]),]
  }else{
	enrich.plot <- enrich.plot[order(enrich.plot[,7]),] %>% .[1:opt$enrichn, ]
  }
 
  if(opt$label.len > 80){
    cnet <- cnetplot(enrich.new, showCategory = 5, vertex.label.cex = 0.1, colorEdge = TRUE, circular = TRUE, cex_label_category = 0.5, cex_label_gene = 0.6)+guides(size = guide_legend(order = 1)) + ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, ceiling(opt$label.len/2)))
  }else{
    cnet <- cnetplot(enrich.new, showCategory = 5, vertex.label.cex = 0.1, colorEdge = TRUE, circular = TRUE, cex_label_category = 0.5, cex_label_gene = 0.6)+guides(size = guide_legend(order = 1)) + ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, opt$label.len))
  } 
  cnet.file <- paste(opt$prefix, "_", type,"_enrich_cnetplot",sep="")
  SavePlot(od = opt$od, filename = cnet.file, data = cnet)
  
   #print(enrich$Description)
  # print(enrich.plot$Description)
  # print(enrich.plot$qvalue)
#  print(enrich.plot$enrich_factor)
#  write.table(enrich.plot, "enrich.plot", row.names=F, col.names=F, quote=F,sep="\t", append=T)
#  dot <- ggplot(enrich.plot, aes(enrich_factor, y=reorder(Description, pvalue))) + 
#    geom_point(aes(colour = pvalue, size = gene_number)) + 
#    scale_colour_gradientn(colours=rainbow(4),guide = "colourbar") + 
#    scale_y_discrete(labels = function(x) str_wrap(x, width = 50) ) +
#    expand_limits(color=seq(0, 0.05, by=0.001)) + 
#    ggtitle(type) + xlab(opt$x.lab) + ylab(opt$y.lab) + 
#    theme_bw() + theme(panel.border = element_rect(colour = "black"), 
#                       plot.title = element_text(vjust = 1), 
#                       legend.key = element_blank(), 
#                       title = element_text(face = "bold", size = opt$lab.size), 
#                       axis.text.x = element_text(face = "bold", size = opt$axis.size-2), 
#                       axis.text.y = element_text(face = "bold", size = opt$axis.size), 
#                       legend.title = element_text(face = "bold", size = opt$legend.size),
#                       legend.text = element_text(size = opt$legend.size))+
#					   guides(size = guide_legend(order = 1))
  dot <- ggplot(enrich.plot, aes(x=reorder(Description,pvalue), y=enrich_factor)) + 
    geom_point(aes(colour = pvalue, size = gene_number)) +
    scale_radius(range = c(1, 4))+ 
    #scale_size(range = c(2, 10))+
    # scale_colour_viridis_c()+
    scale_colour_gradient(low = "OrangeRed", high = "MidnightBlue") + 
    #scale_color_distiller(palette ="PRGn")+
    #scale_colour_gradientn(colours=rainbow(4),guide = "colourbar") + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = opt$label.len)) +
    expand_limits(color=seq(0, 0.05, by=0.001)) + 
    ggtitle(type) + ylab(opt$x.lab) + xlab(opt$y.lab) + 
    theme_bw() + theme(panel.border = element_rect(colour = "black"), 
                       plot.title = element_text(vjust = 1), 
                       legend.key = element_blank(), 
                       title = element_text(face = "bold", size = opt$lab.size-5,family="serif"), 
                       axis.text.x = element_text(family="serif",angle=315,hjust=0,vjust=1,face = "plain", color="black",size = opt$axis.size-2.5), 
                       axis.text.y = element_text(family="serif",face = "plain", color="black",size = opt$axis.size-2), 
                       legend.title = element_text(face = "bold", size = opt$legend.size-6),
					   legend.key.size = unit(8, "pt"),
                       legend.text = element_text(family="serif",size = opt$legend.size-6))+
					   guides(size = guide_legend(order = 1))
  dot.file <- paste(opt$prefix, "_", type, "_enrich_dotplot", sep="")
  SavePlot(od = opt$od, filename = dot.file, data = dot)
}

if(!is.null(opt$kegg)){
  types <- c("Biological Process","Molecular Function","Cellular Component")
  all.outfile <- paste(opt$od, "/", opt$prefix, "_", "GO_enrich.list", sep="")
  header <- c("Ontology", "ID", "Description", "GeneRatio", "BgRatio", "enrich_factor", "pvalue", "qvalue", "geneID", "gene_number")
  write.table(t(header), all.outfile, row.names=F, col.names=F, quote=F, sep="\t")
  sapply(types, function(x){
    index <- which(GO[,4] == x)
    if(length(index) > 0){
      info <- GO[index,1:3]
      EnrichAnalysis(genes, x, info)
    }
	return("Done")
  })

  plot.data <- read.csv(all.outfile, header = TRUE, sep = "\t")
  PrivateBarplot(plot.data, top = opt$enrichn)
  
}
if(!is.null(opt$kegg)){
  EnrichAnalysis(genes, "KEGG pathway", KEGG)
}


