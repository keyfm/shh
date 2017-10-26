#!/usr/bin/env Rscript

## USAGE:
## ./plot_srss.r /my/folder/with/snptable spl1,spl2,spl3
## Script requires the folder that contains the snptable.tsv file. No trailing '/'. Second argument is the individual ID(s) that shall be plotted. Several IDs should be separated by ','.
## Output will be written to input folder: site_read_support_spectrum_ID.pdf

args = commandArgs(trailingOnly=TRUE)
folder=args[1]
colIDs=args[2]

d=read.table(paste(folder,'snpTable.tsv',sep="/"),header=T,stringsAsFactors=F,sep="\t",check.names = FALSE)
colIDs <- strsplit(colIDs,",")[[1]]
for (colID in colIDs){
    d2 <- as.numeric(gsub("\\(|\\)","",substring(grep(" ",d[,colID],invert=F,value=T),3)))/100
    d3 <- d2[ d2 != 1 ]
    pdf(paste(folder,'/site_read_support_spectrum_',colID,".pdf",sep=""))
    hist(d3,xlim=c(0,1),main=paste(colID,": SNPs w/o 100% read support",sep="") ,xlab="read support",col='grey',ylab='# SNPs')
    legend('topleft',legend=c(paste('All Sites:',length(d2)),paste('Var Sites:',length(d3))))
    dev.off()
}

