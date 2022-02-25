#!/usr/bin/env Rscript

library(data.table)
library(igraph)
library(getopt)
#library(ggsci)

spec = matrix(c(
        'help' , 'h', 0, "logical","help description",
        'outdir', 'o', 2, "character", "graph file output dir, default [./]",
        'indir', 'i', 1, "character","input file eg node target weigth",
        'cfg', 'c', '1', 'character', 'configure file',
        'layout', 'l', '2', 'numeric', 'layout type' 
    ), byrow=TRUE, ncol=5)
args <-  getopt(spec)

if (!is.null(args$help) || is.null(args$indir)) {
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

if(is.null(args$outdir)) args$outdir <- '.'
if(is.null(args$layout)) args$layout <- 1
layout_types <- list(layout.auto, layout_with_fr, layout_with_kk, layout_as_star, 
                  layout_in_circle, layout_on_sphere, layout_randomly
                  )

plot_ppi <- function(diffgenes, ppi){
    #deg <- fread("M1.statistic/M1.cluster11.diff_featuregene.xls", sep='\t', select=1:4)
    pdf_file <- file.path(outdir, sub('\\.diff_featuregene.*$', '\\.ppi\\.cytoscapeInput\\.pdf', basename(diffgenes)))
    png_file <- sub('\\.pdf', '\\.png', pdf_file)
    deg <- fread(diffgenes, sep='\t', select=1:4)
    deg <- deg[!duplicated(deg[,2]),]
    setkeyv(ppi, '#Query_id1')
    ppi <- ppi[deg$ID, nomatch=NULL]
    setkeyv(ppi, 'Query_id2')
    ppi <- ppi[deg$ID, nomatch=NULL]
    keep_cols <- c('symbol1', 'symbol2', 'Type', 'Mode', 'Score')
    if(dim(ppi)[1]<1){
        write.table(ppi[, ..keep_cols], file=sub('\\.pdf', '\\.sig', pdf_file), sep='\t', col.names=T, row.names=F, quote=F) 
        return(0)
    }
    ppi <- ppi[!duplicated(ppi[, c('#Query_id1', 'Query_id2', 'Mode')]), ]
    ppi <- ppi[!duplicated(t(apply(as.matrix(ppi[, c("#Query_id1", "Query_id2", "Mode")]), 1, sort)))]
    setkey(deg, ID)
    nodes <- deg[unique(c(ppi$`#Query_id1`, ppi$Query_id2)), .(symbol, Pvalue, log2FC)]
    ppi_net <- graph_from_data_frame(ppi[, ..keep_cols], directed = F, vertices = nodes)
    
    n_mode <- sort(unique(edge_attr(ppi_net)$Mode))
    #n <- pal_igv()(10)[3:(3+length(n_mode))]
    n <- c('#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF') #[1:length(n_mode)]
    names(n) <- c('activation', 'binding', 'catalysis', 'expression', 'inhibition', 'ptmod', 'reaction')
    E(ppi_net)$color <- n[edge_attr(ppi_net)$Mode]

    E(ppi_net)$width <- edge_attr(ppi_net)$Score/max(edge_attr(ppi_net)$Score)
    #V(ppi_net)$size <- degree(ppi_net)^0.9+2 #^(max(0.8, 1-ecount(ppi_net)/1000)) + 0.5
    V(ppi_net)$size <- 12*(degree(ppi_net)/max(degree(ppi_net)))^0.9+5
    V(ppi_net)$frame.color <- NA
    V(ppi_net)$label.font <- 1
    V(ppi_net)$label.cex <- 0.3
    V(ppi_net)$label.color <- 'black'
    #V(ppi_net)$color <- cols_scale(vertex_attr(ppi_net)$log2FC)
    #V(ppi_net)$color <- ifelse(vertex_attr(ppi_net)$log2FC>0, "#FF7F0EFF", "#1F77B4FF") #"red", "RoyalBlue")
    V(ppi_net)$color <- ifelse(vertex_attr(ppi_net)$log2FC>0, "salmon", "deepskyblue")

    pdf(pdf_file)
    if(args$layout == 0){
      weights <- if(ecount(ppi_net)>150) (1/edge_attr(ppi_net)$Score)^0.5 else NULL
      #weights <- NULL
      plot(ppi_net, layout=layout_with_fr(ppi_net, niter=500, weights=weights), mark.col=NA, mark.border=NA)
    }else{
        plot(ppi_net, layout=layout_types[[args$layout]], mark.col=NA, mark.border=NA)
    }
    legend(x=-1.4, y=1.3, c("down", "up"), pch=21,
           pt.bg=c("deepskyblue", "salmon"), pt.cex=1.6, box.lwd=0, pt.lwd=0, bty='n', cex=0.7)
    legend(x=1., y=1.3, c('activation', 'binding', 'catalysis', 'expression', 'inhibition', 'ptmod', 'reaction'), lty='solid',
            col=c('#F0E685FF','#466983FF','#BA6338FF','#5DB1DDFF','#802268FF','#6BD76BFF','#D595A7FF'), cex=0.6, box.lwd=0, box.lty='blank', border=NA)
    dev.off()
    write.table(ppi[, ..keep_cols], file=sub('\\.pdf', '\\.sig', pdf_file), sep='\t', col.names=T, row.names=F, quote=F)
    #system(sprintf("convert  -pointsize 72 -density 600 %s -quality 90 %s",pdf_file,png_file))
    system(sprintf("convert -density 300 -alpha remove %s  %s", pdf_file, png_file))
}

cols_scale = circlize::colorRamp2(c(-2, 0, 2), c("RoyalBlue", "white", "red"))
ppi_file <- unlist(strsplit(grep("^PPI\t", readLines(args$cfg), value=T), '\t'))[2]
ppi_data <- fread(ppi_file, sep='\t')

diffGenesFiles <- list.files(path=args$indir, pattern='*diff_featuregene.xls', full.names=T)
outdir <- file.path(args$outdir, sub("\\.statistic", '\\.ppi_result', basename(args$indir)))
if(!dir.exists(outdir)) dir.create(outdir, recursive=T)
for(s in diffGenesFiles){
  plot_ppi(s, ppi_data)
}
