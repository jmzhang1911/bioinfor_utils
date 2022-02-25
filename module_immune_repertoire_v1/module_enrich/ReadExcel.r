#!/usr/bin/env Rscript
#-----------------------------------------------------------------
# getting script path
#-----------------------------------------------------------------
#getScriptPath <- function(){
#    cmd.args <- commandArgs()
#    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#    script.dir <- dirname(regmatches(cmd.args, m))
#    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
#    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
#    return(script.dir)
#}
#scr_path <- getScriptPath()
scr_path<-dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# get env path
.libPaths(paste0(scr_path,'/readxl_old_lib'))


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript ReadExcel.r input.xls(x) output/dir/ sheetname1 sheetname2 ...")
	print("1) input.xls(x): Forced, intput excel file")
	print("2) output/dir/ : Forced, output directory")
	print("3) sheetname : Optional, wanted sheet name(s); output all sheets if sheetname not given.")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) <2 ) {
	print(args)
	usage()
	stop("the length of args < 2")
}

infile <- args[1]
outdir <- args[2]
#-----------------------------------------------------------------
# main
#-----------------------------------------------------------------
library(readxl)
# read excel
print("read excel ...")

sheets <- excel_sheets(path=infile)

if ( length(args) >2 ){
	wantedSheets <- args[c(-1,-2)]
	outsheets <- intersect(sheets,wantedSheets)
}else {
	outsheets <- sheets
}

if ( length(outsheets)==0 ){
	stop("can not find wanted sheet(s)")
}

getOutName <- function(str=NULL,num=NULL){
	# check NULL
	if( is.null(str) ) stop("str is NULL")
	if( is.null(num) ) stop("num is NULL")
	
	# encoding str
	oldenc <- Encoding(str)
	if ( oldenc=='unknown' ){
		Encoding(str)<-'UTF-8'
	}else if ( oldenc != 'UTF-8' ){
		str <- iconv(str, oldenc, 'UTF-8')
	}
	Encoding(str) <- 'UTF-8'

	out <- gsub(' ','_',str)
	# is str chinese?
	# regular expression for chinese word: [\u4e00-\u9fa5]
	if ( grepl('[^a-zA-Z0-9._-]',out) ){ # chinses or special character
		out <- paste0('sheet',num)
	}

	return(out)
}

for ( i in 1:length(sheets) ){
	if ( sheets[i] %in% outsheets ) {
		dat <- read_excel(path=infile,sheet=i,col_names=F,col_types='text')
		outname <- getOutName(sheets[i],i)
		write.table(dat,file=paste0(outdir,'/',outname),quote=F,sep='\t',eol='\n',row.names=F,na='',col.names=F,fileEncoding='UTF-8')
		print(paste0("write sheet [",sheets[i],"] to [",outname,"] done."))
	}
}

print("read excel done.")
