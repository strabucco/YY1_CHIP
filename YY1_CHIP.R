#R script for downloading and analyzing ChIP data
#This will be used to compare ChIP binding between YY1 and EZH2 in B-cells

#load packages needed
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(dplyr)

#load YY1 ChIP-Seq data and load peaks
YY1File<-readPeakFile("wgEncodeSydhTfbsGm12878Yy1StdPk.narrowPeak.bed")
#annotate YY1 file
YY1Anno<-annotatePeak(YY1File, tssRegion = c(-3000, 3000), 
                      TxDb =TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
YY1Names<-as.data.frame(YY1Anno)

#download and load EZH2 file
EZH2File<-read.csv("EZH2Table.csv")
#filter file by EZH2 peaks in CB but not ESCs, also remove those without peaks
EZH2Only<-EZH2File[-c(1,3,5, 7:11)]
EZH2CB<-EZH2Only[grepl("CB,EZH2",EZH2Only$CB.EZH2.bound),] #list of all with EZH2 peaks in centroblasts
EZH2CBOnly<-EZH2CB[!grepl("hESC,EZH2", EZH2CB$hESC.EZH2.bound),] #list of genes with EZH2 peaks ONLY
#in centroblasts-so this excludes genes that also appear in hESCs

head(EZH2CBOnly)
#also filter file by EZH2 peaks in CB (all)


#compare gene names from YY1 annotation to EZH2 filtered files ($ORF)
#get number of gene names overlapping and a CSV or output list of names