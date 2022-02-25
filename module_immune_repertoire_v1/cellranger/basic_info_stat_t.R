#!/usr/bin/env Rscript

#import packages
library(getopt)
suppressMessages(library(tidyverse))

# get options
spec <- matrix(c(
		'files' , 'f', 1, 'character',
		'od'  , 'o', 1, 'character',
		'help', 'h', 0, 'logical'
		), byrow = TRUE, ncol = 4
)

opt <- getopt(spec)
print.usage <- function(){
	cat("ProgramName:basic_info_stat.R
Version: Version v1.0
Contact: zhengly <zhengly@biomarker.com.cn> 
Program Date:   2020.12.3
Description: this program is used to statistic information ;\n")
	cat(getopt(spec, usage = TRUE))
	q(status = 1)
}

if(!is.null(opt$help)){
	print.uasge()
}

if(is.null(opt$file)){
	print.usage()
}

if(is.null(opt$od)) {opt$od <- "./"}

if(!file.exists(opt$od)){
	dir.create(opt$od)
}

files <- unlist(strsplit(opt$files, split = ","))
header1 <- matrix(c("sampleID", "Number of Read", "Valid Barcodes", "Q30 Bases in Barcode", "Q30 Bases in RNA Read 1",  "Q30 Bases in UMI"), nrow = 1)
header2 <- matrix(c("sampleID", "Estimated Number of Cells", "Number of Cells With Productive V-J Spanning Pair", "Mean Reads per Cell", "Mean Used Reads per Cell","Fraction Reads in Cells"), nrow = 1)
header3 <- matrix(c("sampleID", "Reads Mapped to Any V(D)J Gene", "Reads Mapped to TRA","Reads Mapped to TRB","Median TRA UMIs per Cell","Median TRB UMIs per Cell","Cells With Productive V-J Spanning (TRA, TRB) Pair","Cells With TRA Contig","Cells With TRB Contig","Cells With CDR3-annotated TRA Contig","Cells With CDR3-annotated TRB Contig","Cells With V-J Spanning TRA Contig","Cells With V-J Spanning TRB Contig"), nrow = 1)
write.table(header1, file = paste(opt$od, "t_total_seqence_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(header2, file = paste(opt$od, "t_total_cell_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F,col.names = F)
write.table(header3, file = paste(opt$od, "t_total_vdj_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F,col.names = F)
sapply(files, function(x){
	sample.name  <- str_split(basename(x), pattern = '\\.')[[1]][1]
	info <- read.csv(x, header = T, check.name = F)
	seq.data <- info %>% select(c("Number of Read Pairs", 
                              "Valid Barcodes",  
                              "Q30 Bases in Barcode", 
                              "Q30 Bases in RNA Read 1", 
                              "Q30 Bases in UMI"))
	seq.data <- cbind(sample.name, seq.data)
	cell.data <- info %>% select(c("Estimated Number of Cells", 
                               "Number of Cells With Productive V-J Spanning Pair",
                               "Mean Read Pairs per Cell", 
                               "Mean Used Read Pairs per Cell",
                               "Fraction Reads in Cells"))
	cell.data <- cbind(sample.name, cell.data)   
	mapped.data <- info %>% select(c("Reads Mapped to Any V(D)J Gene", 
                                 "Reads Mapped to TRA",
                                 "Reads Mapped to TRB",
                                 "Median TRA UMIs per Cell",
                                 "Median TRB UMIs per Cell",
                                 "Cells With Productive V-J Spanning (TRA, TRB) Pair",
                                 "Cells With TRA Contig",
                                 "Cells With TRB Contig",
                                 "Cells With CDR3-annotated TRA Contig",
                                 "Cells With CDR3-annotated TRB Contig",
                                 "Cells With V-J Spanning TRA Contig",
                                 "Cells With V-J Spanning TRB Contig"
                                 ))
	mapped.data <- cbind(sample.name, mapped.data)
	write.table(seq.data, file = paste(opt$od, "t_total_seqence_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F, append = TRUE)
	write.table(cell.data, file = paste(opt$od, "t_total_cell_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F, append = TRUE)
	write.table(mapped.data, file = paste(opt$od, "t_total_vdj_info_stat.xls", sep = "/"), sep = "\t", quote = F, row.names = F, col.names = F, append = TRUE)
})




















