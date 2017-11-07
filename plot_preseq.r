#!/usr/bin/env Rscript

## USAGE:
## ./plot_preseq.r /my/path/file.lcextrap outfolder/ sampleID MaxReads
## SampleID: used for header in pdf and file name.
## MaxReads: Max number of reads (sequencing effort) to be plotted.
## EAGER output summary incl. sample ID (.csv file). Required.
## NOTE: For each CF I plot the amount of distinct reads expected at this CF, and the amount of extra sequencing necessary to get there (considering %endo and already performed sequencing)
## NOTE2: Present data-state shown as red dot. Derived from EAGER out table as well, not from preseq results --> leads to minor deviations sometimes.

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfolder <- args[2]
sampleID <- args[3]
if( !is.na(args[4]) ){maxR <- as.integer(args[4])} else {maxR <- 20000000}
if( !is.na(args[5]) ){ endo <- read.csv(args[5],check.names=F);rownames(endo) <- endo[,'Sample Name'] }
## Test if procvided endo file has actual sample name
if( sampleID %in% rownames(endo) ){
	current.reads <- endo[ sampleID, 'Mapped Reads after RMDup' ]
	endo.per <- endo[ sampleID, 'Endogenous DNA (%)' ]
	current.cf <- round(endo[sampleID,'# mapped reads prior RMDup']/endo[sampleID,'Mapped Reads after RMDup'],2)
	sequencing.done <- endo[sampleID,'# reads after C&M prior mapping'] ## use this number to calculate how much seq has to be done. Note: removal of reads due to Q cut off not considered
} else {
	current.reads <- NA
	endo.per <- NA
	warning('SampleID not found in EAGER summary. Current read count and % endogenous cannot be shown in plot.')
}

## format data
d1=read.table(infile, header=T,sep="\t")
d2=d1[d1[,1]<=maxR,]

## plot
maxY <- max(d2[,'UPPER_0.95CI'])+(0.1*max(d2[,'UPPER_0.95CI']))
pdf(paste(outfolder,sampleID,"_preseq.pdf",sep=""))
plot("",ylim=c(0,maxY),xlim=c(0,maxR),xlab="Mapped reads prior RMDup",ylab="Mapped Reads after RMDup (w/ 95CI)",main=sampleID,cex.main=0.95)
mtext(paste("% Endo:",endo.per))
polygon(c( d2[,'TOTAL_READS'] , rev(d2[,'TOTAL_READS']) ) , c( d2[,'LOWER_0.95CI'] , rev(d2[,'UPPER_0.95CI']) ) ,col="lightblue1",border="lightskyblue4",lwd=1.5)
grid(NULL,NULL,lty=2,col="grey")
lines(d2[,'TOTAL_READS'],d2[,'EXPECTED_DISTINCT'],col='blue',lwd=2)

## extract/plot CFs
cluster.fac <- d2[-1,'TOTAL_READS']/d2[-1,'EXPECTED_DISTINCT'];my.idx <-vector(length=4)
for (i in c(1:4)){
	c = c(2.5,3,3.5,4)[i]
	my.idx[i] <- which(abs(cluster.fac-c)==min(abs(cluster.fac-c)))
}
points(d2[my.idx,'TOTAL_READS'],d2[my.idx,'EXPECTED_DISTINCT'],pch=21,col='black',bg=colorRampPalette(c('darkolivegreen1','darkgreen'))(4))
## current seq depth
points(endo[sampleID,'# mapped reads prior RMDup'],endo[sampleID,'Mapped Reads after RMDup'],pch=21,col='black',bg='red')
## legend
## NOTE: For each CF I plot the amount of distinct reads expected at this CF, and the amount of extra sequencing necessary to get there (considering %endo and already performed sequencing)
legend('bottomright',legend=paste('CF:',c('2.5','3.0','3.5','4.0'),'; #Dist.:',round(d2[my.idx,'EXPECTED_DISTINCT']/1000000,2),'M; #SeqLibEx:',round(((d2[my.idx,'TOTAL_READS']*(100/endo.per))-sequencing.done)/1000000,2),'M',sep=""),fill=colorRampPalette(c('darkolivegreen1','darkgreen'))(4))
dev.off()
