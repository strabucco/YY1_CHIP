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
#filter YY1 list by only promoters for comparison because EZH2 list is based on ChIP-CHIP,
#which only used probes in the promoter region, so we really are only able to compare what ChIPs
#to promoter regions anyway.
YY1Promoter<-YY1Names[grepl("Promoter", YY1Names$annotation),]
YY1PromoterClean<-YY1Promoter[!duplicated(YY1Promoter$SYMBOL),]
YY1noNA<-YY1PromoterClean[!(is.na(YY1PromoterClean$SYMBOL)), ]


#download and load EZH2 file
EZH2File<-read.csv("EZH2Table.csv")
#filter file by EZH2 peaks in CB but not ESCs, also remove those without peaks
EZH2Only<-EZH2File[-c(1,3,5, 7:11)]
#also filter file by EZH2 peaks in CB (all)
EZH2CB<-EZH2Only[grepl("CB,EZH2",EZH2Only$CB.EZH2.bound),] #list of all with EZH2 peaks in centroblasts
EZH2CBOnly<-EZH2CB[!grepl("hESC,EZH2", EZH2CB$hESC.EZH2.bound),] #list of genes with EZH2 peaks ONLY
#in centroblasts-so this excludes genes that also appear in hESCs
EZH2CBClean<-EZH2CB[!duplicated(EZH2CB$ORF),]
EZH2CBnoNA<-EZH2CBClean[!(is.na(EZH2CBClean$ORF)),]
EZH2CBOnlyClean<-EZH2CBOnly[!duplicated(EZH2CBOnly$ORF),]
EZH2CBOnlynoNA<-EZH2CBOnlyClean[!(is.na(EZH2CBOnlyClean$ORF)),]

#compare gene names from YY1 annotation to EZH2 filtered files ($ORF)
bothCBOnly<-intersect(EZH2CBOnlynoNA$ORF, YY1noNA$SYMBOL)
str(bothCBOnly)
bothAll<-intersect(EZH2CBnoNA$ORF,YY1noNA$SYMBOL)
str(bothAll)
head(bothAll)
#get number of gene names overlapping and a CSV or output list of names