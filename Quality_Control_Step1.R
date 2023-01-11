
library(edgeR)
library(DESeq2)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(GEOquery)
library(tibble)
library(tidyverse)
library(Rsubread)
library(airway)
library(Rsamtools)
library(GenomicFeatures)
library(BiocParallel)
library(rtracklayer)
library(Rgff)
library(tximport)
library(tximeta)
library(RColorBrewer)
library(RNAseqQC)
library(purrr)
library(ensembldb)
library(fastqcr)
library(Rqc)
library(BiocParallel)

Illumina_read_data_ <- read_csv("Illumina read data .csv")
View(Illumina_read_data_)

Il<-Illumina_read_data_%>%
  mutate(Sample_Name = factor(Sample_Name))

Total_Read_Pairs 

ggplot(data=Il, aes(x=Sample_Name, y=Total_Read_Pairs))+
  geom_bar(stat="identity", color = "black", fill="lightblue")+
  ylab("Total # of Read Pairs")+
  geom_text(aes(label = Total_Read_Pairs, vjust = 2))+
  xlab("Samples")+
  theme_linedraw()+
  theme(text=element_text(size=16,  family="Times New Roman"))

 ggplot(data=Il, aes(x=Sample_Name, y=Total_Reads))+
  geom_bar(stat="identity", color = "black", fill="red")+
  ylab("Total Reads")+
  geom_text(aes(label = Total_Reads, vjust = 2))+
  xlab("Samples")+
  theme_linedraw()+
  theme(text=element_text(size=16,  family="Times New Roman"))
  
 ggplot(data=Il, aes(x=Sample_Name, y=BP_30))+
   geom_bar(stat="identity", color = "black", fill="Purple")+
   ylab(" % of Base Pairs in Q3")+
   geom_text(aes(label = BP_30, vjust = 2))+
   xlab("Samples")+
   theme_linedraw()+
   theme(text=element_text(size=16,  family="Times New Roman"))






