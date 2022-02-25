#!/usr/bin/env Rscript
library(getopt)
# parameter setting
spec <- matrix(c(
  'help'        , 'h', 0, "logical"  , "print help usage",
  'deg'         , 'd', 1, "character", "character  input , cluster diff marker gene file",
  'od'          , 'o', 1, "character", "enrich.out dir, default current path",
  'prefix'      , 'P', 1, "character", "character  enrich.output key,default cluster",
  'species'     , 's', 1, "character", "Species, exmple: human、mouse、rat、fly、celegans、zebrafish、arabidopsis",
  'pvalue'      , 'p', 1, "numeric"  , "pvalue cut off, default: 0.05",
  'qvalue'      , 'q', 1, "numeric"  , "qvalue cut off, defalut: 0.2",
  'enrichn'     , 'n', 1, "integer"  , "int the first enrichn term would be plotted, default 20",
  "color"       , 'c', 1, "character", "character       One of pvalue p.adjust qvalue, default pvalue",
  'label.len'   , 'N', 1, "integer"  , "label length, default:100",
  'size'   , 'a', 1, "integer"  , "axis size"
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:reactome_analysis.R
  Version: Version v1.0
  Description: this program is used to do enrich analysis ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)
}

if(!is.null(opt$help))   PrintUsage()
if(is.null(opt$deg) | is.null(opt$species))       PrintUsage()
if(is.null(opt$od))          { opt$od <- "./" }
if(is.null(opt$prefix))  { opt$prefix <- "cluster" }
if(is.null(opt$color))   { opt$color <- "pvalue" }
if(is.null(opt$pvalue)) { opt$pvalue <- 0.05 }
if(is.null(opt$qvalue)) { opt$qvalue <- 0.2 }
if(is.null(opt$enrichn)) { opt$enrichn <- 20 }
if(is.null(opt$label.len)) { opt$label.len <- 80 }
if ( is.null(opt$size ) )          { opt$size <- 7 }
if(!dir.exists(opt$od))  {dir.create(opt$od)}
library(tidyverse)
library(ReactomePA)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.At.tair.db)
library(org.Dr.eg.db)


# diff.exp <- read.table("statistic/cluster0.diff_featuregene.xls", sep = "\t", header = TRUE)
od <- opt$od

if(opt$species == "Homo_sapiens"){
  g2Ensembl <- toTable(org.Hs.egENSEMBL)
  g2Symbol <- toTable(org.Hs.egSYMBOL)
  species <- "human"
}else if(opt$species == "Mus_musculus"){
  g2Ensembl <- toTable(org.Mm.egENSEMBL)
  g2Symbol <- toTable(org.Mm.egSYMBOL)
  species <- "mouse"
}else if(opt$species == "Rattus_norvegicus"){
  g2Ensembl <- toTable(org.Rn.egENSEMBL)
  g2Symbol <- toTable(org.Rn.egSYMBOL)
  species <- "rat"
}else if(opt$species == "Danio_rerio"){
  g2Ensembl <- toTable(org.Dr.egENSEMBL)
  g2Symbol <- toTable(org.Dr.egSYMBOL)
  species <- "zebrafish"
}else if(opt$species == "Caenorhabditis_elegans"){
  g2Ensembl <- toTable(org.Ce.egENSEMBL)
  g2Symbol <- toTable(org.Ce.egSYMBOL)
  species <- "celegans"
}else if(opt$species == "Drosophila_melanogaster"){
  g2Ensembl <- toTable(org.Dm.egENSEMBL)
  g2Symbol <- toTable(org.Dm.egSYMBOL)
  species <- "fly"
}else if(opt$species == "Arabidopsis_thaliana"){
  g2Ensembl <- toTable(org.At.tairENSEMBL)
  g2Symbol <- toTable(org.At.tairSYMBOL)
  species <- "arabidopsis"
}else{
  q()
}
EG2List <- inner_join(g2Ensembl, g2Symbol, by = "gene_id")

GetGene <- function(vectors, data){
  conv.gene <- vector()
  conv.id <- vector()
  for(i in 1:length(vectors)) {
    geneID <- strsplit(vectors[i], split = "/") %>% unlist
    genes <- data[data$gene_id %in% geneID, "symbol"]
    ids <- data[data$gene_id %in% geneID, "ensembl_id"]
    Gene <- paste(genes, collapse = "/")
    Ids <- paste(ids, collapse = "/")
    conv.gene <- c(conv.gene, Gene)
    conv.id <- c(conv.id, Ids)
  }
  re <- data.frame(geneID = conv.id, geneSymbol = conv.gene)
  return(re)
}

data <- read.table(opt$deg, sep = "\t", header = TRUE)
# mouse: org.Mm.eg.db ; human:org.Hs.eg; rat: org.Rn.eg.db
# deal with data
data1 <- inner_join(data, EG2List[,1:2], by = c("ID" = "ensembl_id")) 
data1 <- data1[order(data1$log2FC, decreasing = TRUE),]
# analysis
enrich <- enrichPathway(data1$gene_id, organism = species, pvalueCutoff = opt$pvalue, qvalueCutoff = opt$qvalue, minGSSize = 2, maxGSSize = 3000)
enrich.re <- enrich@result
mygene <- GetGene(enrich.re$geneID, EG2List)
my.re <- data.frame(enrich.re[1:7], mygene, Count = enrich.re$Count)
colnames(my.re) <- c("ReactomeID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "geneSymbol", "Count")
write.table(file = file.path(od, paste(opt$prefix, "_reactome_enrich.list", sep="")), my.re, quote = FALSE, sep = "\t", row.names = FALSE)
enrich.plot <- my.re[order(my.re[, opt$color]),] %>% .[1:opt$enrichn, ]
dot <- ggplot(enrich.plot, aes(x=reorder(Description, pvalue), y = GeneRatio)) +
  geom_point(aes(colour = pvalue, size = Count)) +
  scale_radius(range = c(1, 4))+
  scale_colour_gradient(low = "OrangeRed", high = "MidnightBlue") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = opt$label.len)) +
  expand_limits(color=seq(0, 0.05, by=0.001)) +
  ggtitle("") + ylab("GeneRatio") + xlab("") +
  theme_bw() + theme(panel.border = element_rect(colour = "black"),
                     plot.title = element_text(vjust = 1),
                     legend.key = element_blank(),
                     title = element_text(face = "bold", size = opt$size,family="serif"),
                     axis.text.x = element_text(family="serif",angle=315,hjust=0,vjust=1,face = "plain", color="black",size = opt$size),
                     axis.text.y = element_text(family="serif",face = "plain", color="black",size = opt$size),
                     legend.title = element_text(face = "bold", size = opt$size + 4),
                     legend.key.size = unit(8, "pt"),
                     legend.text = element_text(family="serif",size = opt$size + 2)) + guides(size = guide_legend(order = 1))
ggsave(file = file.path(od, paste(opt$prefix, "_reactome_enrich_dotplot.png", sep="")), dot, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste(opt$prefix, "_reactome_enrich_dotplot.pdf", sep="")), dot, width = 6, height = 4, scale = 1.3)
log10Pvalue <- -(log(enrich.plot$pvalue)/log(10))
bar.data <- data.frame(Description = enrich.plot$Description, log10Pvalue = log10Pvalue, gene_number = enrich.plot$Count)
p2 <- ggplot(bar.data , aes(Description, log10Pvalue)) + geom_col(fill = "#8DA0CB", position = 'dodge', width = 0.6, show.legend = FALSE) + scale_x_discrete(labels = function(x) str_wrap(x, width = opt$label.len)) + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(size = 0.4), axis.line.y = element_line(size = 0.4), axis.text.x = element_text(family="serif", size = 6, angle = 300, hjust=0,vjust=1), axis.title.x = element_blank(), axis.text.y = element_text(family="serif", size = opt$size), legend.title = element_blank(), axis.title.y = element_text(family = "serif", size = opt$size + 2)) + ylab("-log10pvalue")
ggsave(file = file.path(od, paste(opt$prefix, "_reactome_enrich_barplot.png", sep="")), p2, width = 6, height = 4, scale = 1.3)
ggsave(file = file.path(od, paste(opt$prefix, "_reactome_enrich_barplot.pdf", sep="")), p2, width = 6, height = 4, scale = 1.3)

# gsea
#mlist <- data1$log2FC
#names(mlist) <- data1$gene_id
#gsea <- gsePathway(mlist, organism = "human")
#gsea.re <- gsea@result
#mygene <- GetGene(enrich.re$geneID, EG2List)
#my.re <- data.frame(enrich.re, mygene)


