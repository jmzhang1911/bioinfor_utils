#!/usr/bin/env Rscript
# #############
library(getopt)
library(R.utils)
library(RColorBrewer)
spec = matrix(c(
        'help' , 'h', 0, "logical","help description",
        'output', 'o', 2, "character", "graph file output dir, default [./]",
        'infile', 'i', 1, "character","input file eg node target weigth",
        'layout' , 'l', 2, "character","the graph layout ,layout_randomly  layout_as_star layout_nicely ,default  [layout.auto]",
        'size' , 's', 2, "integer","node size default [degree]",
        'title', 't', 1, "character", "graph title, must be given",
        'node.col', 'n', 2, "integer", "graph title, default [1]",
        'target.col', 'g', 2, "integer", "graph title, default [3]",
        'weight', 'w', 0, "logical","set edge size depend on the weight,default [FALSE]",
        'color.type' , 'r', 2, "integer","the type of ppi network pic color [optional, default: 1]",
        'color.node' , 'a', 2, "character","the color of node from regulate network pic  [optional, default: red]",
        'color.target' , 'b', 2, "character","the color of target from regulate network pic [optional, default: skyblue]",
        'label.cex', 'x', 2, "numeric", "node label font size ,default [0.2]"
		),byrow=TRUE, ncol=5)
opt = getopt(spec)

# define usage function
print_usage <- function(spec=NULL){
	#cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
1) Rscript plot_ppi_network.r --infile ppi_cytoscapeInput.sif --title PPI_network --weight
2) Rscript plot_ppi_network.r --infile lncRNA_mRNA_regulate.cytoscapeInput.sif --output ./ \\
 --title lncRNA_mRNA_regulate --color.node red --color.target yellow
3) Rscript plot_ppi_network.r --infile lncRNA_mRNA_regulate.cytoscapeInput.sif --output ./  \\
--title lncRNA_mRNA_regulate --color.node red --color.target yellow --label.cex 0.5


Options:
    -h|--help            help description
    -o|--output          graph file output dir, default [./]
    -i|--infile          input file eg node target weigth, must be given
    -t|--title           graph title, must be given
    -l|--layout          the graph layout ,layout_randomly  layout_as_star layout_nicely ,default  [layout.auto]
    -s|--size            node size default degree
    -w|--weight          set edge size depend on the weight,default [TRUE]
    -x|--label.cex       node label font size ,default [0.2]
    -n|--node.col        node col  title,  default [1]
    -g|--target.col      target col, default [2]

    -r|--color.type      the type of ppi network pic color ,only for ppi network [optional, default: 1]
    -a|--color.node      the color of node from regulate network pic , [optional, default: red]
    -b|--color.target    the color of target from regulate network pic [optional, default: skyblue]
    ##########attention
    the color.type  for ppi network,if draw regulate network please give the color.node and color.target choised
\n")
	q(status=1);
}


#############
if ( !is.null(opt$help) | is.null(opt$infile) | is.null(opt$title) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
input<-as.vector(opt$infile)
if(file.exists(input)){
  input=getAbsolutePath(input)
}else{
  stop("Error:input file not exist!")
}
od<-getAbsolutePath(opt$output)
#####################
#+--------------------
# some default options
#+--------------------
if ( is.null(opt$node.col) ) opt$node.col <- 1
if ( is.null(opt$target.col) ) opt$target.col <- 2
if ( is.null(opt$layout) ) opt$layout <- 'layout.auto'
if ( is.null(opt$weight) ) opt$weight <- 'FALSE'
if ( is.null(opt$label.cex)) opt$label.cex <- 0.2
if ( is.null(opt$output) ) opt$output <- "./"
type<-"ppi"
if ( (is.null(opt$color.type)) &&  (is.null(opt$color.node)) && (is.null(opt$color.target))){type <- "ppi"}
if ( (!is.null(opt$color.type)) &&  (is.null(opt$color.node)) && (is.null(opt$color.target))){type <- "ppi"}
if ( (is.null(opt$color.type)) &&  (!is.null(opt$color.node) || !is.null(opt$color.target))){type <- "regulate"}
if ( (!is.null(opt$color.type) && !is.null(opt$color.target)) ||  (!is.null(opt$color.type) && !is.null(opt$color.target))){
    cat("make sure to draw ppi network pic or regulate pic \n")
    print_usage(spec)
    stop("check the para !!!\n")
}
if (type =="ppi"){
    if ( is.null(opt$color.type) ) opt$color.type <- 1
    if ( (opt$color.type < 1) || (opt$color.type > 6) ) {
        cat("Final Error: color.type must be 1-6\n")
        print_usage(spec)
    }
}
if (type =="regulate"){
    if ( is.null(opt$color.node) ) opt$color.node <- "red"
    if ( is.null(opt$color.target) ) opt$color.target <- "skyblue"
    if (opt$color.node==opt$color.target){
        print_usage(spec)
        #stop("the node color cannot be same as the target color  !!!\n")
    }
}
####################
pdf_to_png<-function(pdf_file,png_file){
	pdf_to_png_cmd<-sprintf("convert  -pointsize 72 -density 600 %s -quality 90 %s",pdf_file,png_file)
	print(pdf_to_png_cmd)
	system(pdf_to_png_cmd)
}
###################
color.1 <- colorRampPalette(rev(c("#00ff00", "#FFFF00FF", "#ff0000")))
color.2 <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4", rep("#4575B4",2))))
color.3 <- colorRampPalette(brewer.pal(3,"Set1"))
color.4 <- colorRampPalette(brewer.pal(3,"Set2"))
color.5 <- colorRampPalette(brewer.pal(3,"Set3"))
color.6 <- colorRampPalette(rev(c("#ff0000", "#ffffff",  "#0000ff")))
color.set <- list(color.1, color.2, color.3, color.4, color.5, color.6)
###################
library(igraph)
ppi<-read.table(input,sep='\t',check.names=F,header=F,stringsAsFactors=F,colClasses = "character")
ppi_net<-graph_from_edgelist(as.matrix(ppi[,c(opt$node.col,opt$target.col)]),directed = F)
node_names<-V(ppi_net)$name
genenames<-as.character(ppi[,opt$node.col])
print(length(node_names))
####计算权重，用于边的粗细
weight <- seq(1,3,length.out=ecount(ppi_net))
line <- length(weight)
print(line)
#####计算连接度，用于节点大小
d <- degree(ppi_net)
print(node_names[which(d==max(d))])
#####节点标签字体大小
cex <- opt$label.cex
print(cex)
V(ppi_net)$label.cex <- cex
if(opt$weight){E(ppi_net)$width<- weight}
V(ppi_net)$frame.color <- NA
 V(ppi_net)$label.color <-"black"
#####set color
if (type =="regulate" ){
    types<-rep('1',length(node_names))
    types[node_names %in% genenames]<-'2'
    types
    col<-c(opt$color.node,opt$color.target)
    V(ppi_net)$color=col[as.numeric(types)]
}
if (type=="ppi"){
    print(opt$color.type)
    V(ppi_net)$color=color.set[[opt$color.type]](max(d))[d]
}

if ( !is.null(opt$size) ){size <- V(ppi_net)$size<-opt$size}else{V(ppi_net)$size<- d}

pdffile<-paste(od,"/",opt$title,".pdf",sep="")
pngfile<-paste(od,"/",opt$title,".png",sep="")
pdf(file=pdffile)
if(opt$layout=="layout.auto"){
    plot(ppi_net,main =opt$title ,layout=layout.auto)
}else if(opt$layout=="layout_randomly"){
    plot(ppi_net,main =opt$title ,layout=layout_randomly)
}else if(opt$layout=="layout_as_star"){
    plot(ppi_net,main =opt$title ,layout=layout_as_star)
}else if(opt$layout=="layout_nicely"){
    plot(ppi_net,main =opt$title ,layout=layout_nicely)
}else{
    plot(ppi_net,main =opt$title ,layout=layout.auto)
}

dev.off()

pdf_to_png(pdf_file=pdffile,png_file=pngfile)
