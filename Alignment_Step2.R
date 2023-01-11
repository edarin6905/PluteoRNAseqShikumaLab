BiocInstaller::biocValid()
BiocManager::install("tximeta", force=TRUE)
aBiocManager::install(c("airway", "Rsamtools", "GenomicFeatures", "BiocParallel"))
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
library(fastqcr)

####Alignment and Counting####
#Do not run following lines of code, only needed to be run once (december 18,2022) all BAM files should be 
#accounted for. 

#Create metadata file 
FileName<-c("r1_S28_R1_001","r1_S29_R1_001","r1_S30_R1_001","r1_S31_R1_001",
            "wt1_S24_R1_001","wt1_S25_R1_001","wt1_S26_R1_001","wt1_S27_R1_001")
SampleName<-c("r1S28","r2S29","r3S30","r4S31","wt1S24","wt1S25","wt1S26","wt1S27")
Condition<-c("Experimental","Experimental","Experimental","Experimental","Wildtype","Wildtype","Wildtype","Wildtype")
metadata<-as.data.frame(cbind(FileName,SampleName,Condition))

#Load in Illumina sequencing data  

fastq.files_1 <- list.files(path = "./data", pattern = "R1_001.fastq$", full.names = TRUE)

fastq.files_2 <- list.files(path = "./data", pattern = "R2_001.fastq$", full.names = TRUE)




#run fastqc on our unprocessed reads, and save 
#the results in a new directory within our base directory

fastqc_dir <- "~Users/emilydarin/Desktop/PludioRNAseq/data"

fastqc(fq.dir = base_dir, qc.dir = fastqc_dir, threads = 4)

qc <- qc_read(file.path(fq_dir, "file_fastqc.zip"))
names(qc)



# 
# #Load in reference genome for P.Luteo taken from NCBI 
# 
buildindex(basename="p.luteo.ref",reference="GCF_000814765.1_ASM81476v1_genomic.fna")
# 
# #Align the reads to the reference genome 
# #The default setting for align is that it only keeps reads that 
# #uniquely map to the reference genome. For testing differential expression of genes, this is what we want
# #Takes 2 hours 
# 
Rsubread::align(index = "p.luteo.ref",
                readfile1 = fastq.files_1,
                readfile2 = fastq.files_2,
                type = "rna",
                input_format = "FASTQ",
                output_format = "BAM",
                PE_orientation = "rf",
                nsubreads = 10,
                TH1 = 3,
                TH2 = 1,
                maxMismatches = 3,
                
                useAnnotation = TRUE,
                annot.ext='GCF_000814765.1_ASM81476v1_genomic.gtf',
                isGTF = TRUE,
                GTF.featureType="gene",
                GTF.attrType="gene_id",
                chrAliases=NULL,
                #level of summarization
                useMetaFeatures=TRUE)




# 
# #Pull curated BAM files formed during alignment
#  #Pull curated BAM files formed during alignment

bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE)
bam.files


# Make proportions 
test<- propmapped(
  
  bam.files,
  countFragments = FALSE,
  properlyPaired = FALSE,
  verbose = FALSE)

test

test<-as.data.frame(test)
view(test)

Sample<-c("r1_S28","r1_S29","r1_S30","r1_S31","wt_S24","wt_S25","wt_S26","wt_S27")
Prop<-cbind(test,Sample)

library(ggplot2)


p<-ggplot(data=Prop, aes(x=Sample, y=PropMapped)) +
  geom_bar(stat="identity", color = "black", fill="green")+
  ylab("Proportion of Mapped Reads")+
  xlab("Sample Name")+
  theme_linedraw()+
  geom_text(aes(label=PropMapped,vjust = 2 ))+
  theme(text=element_text(size=16,  family="Times New Roman"))
  
p



# 