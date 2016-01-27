#R script for downloading and analyzing ChIP data
#This will be used to compare ChIP binding between YY1 and EZH2 in B-cells

#load packages needed
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#load YY1 ChIP-Seq data and load peaks
YY1File<-readPeakFile("wgEncodeSydhTfbsGm12878Yy1StdPk.narrowPeak.bed")
#annotate YY1 file
YY1Anno<-annotatePeak(YY1File, tssRegion = c(-3000, 3000), 
                      TxDb =TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
YY1Names<-as.data.frame(YY1Anno)
head(YY1Names)
#download and load EZH2 file
#filter file by EZH2 peaks in CB but not ESCs, also remove those without peaks
#also filter file by EZH2 peaks in CB (all)


#compare gene names from YY1 annotation to EZH2 filtered files ($ORF)
#get number of gene names overlapping and a CSV or output list of names