
BiocInstaller::biocValid()
BiocManager::install("Rqc", force=TRUE)
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
library(ggplot2) 
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(extrafont)

font_import()
loadfonts(device = "win")
 
#Pull curated BAM files formed during alignment
#Compiles any file that has a "BAM" in it. 
bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE) 
bam.files #Check that all of your BAM files are accounted for 

#Counting reads 
#annotation text is your GTF file 
#select true for is GTF annotation file 
#select GTF feature type which is gene (bacteria don't have exons)
#GTF attrType identify ninth column by gene_id 
#You want paired end data with meta features 
 
  fc<-featureCounts(bam.files,
              annot.ext='GCF_000814765.1_ASM81476v1_genomic.gtf',#annotation text is your GTF file 
              isGTFAnnotationFile=TRUE, #select true for is GTF annotation file 
               GTF.featureType="gene",#select GTF feature type which is gene (bacteria don't have exons)
               GTF.attrType="gene_id", #GTF attrType identify ninth column by gene_id
             chrAliases=NULL,isPairedEnd=TRUE,#You want paired end data with meta features 
  useMetaFeatures=TRUE) #Data summary level 

  counters<-fc$counts #Pull just the counts
  colnames(counters) #Check column names of the dataframe

 
 #Create metadata file 
 
 FileName<-c("r1_S28_R1_001","r1_S29_R1_001","r1_S30_R1_001","r1_S31_R1_001",
             "wt1_S24_R1_001","wt1_S25_R1_001","wt1_S26_R1_001","wt1_S27_R1_001") #Names of the files you are loading in
 SampleName<-c("r1S28","r2S29","r3S30","r4S31","wt1S24","wt1S25","wt1S26","wt1S27") #Names you would like them to be called 
 Condition<-c("Experimental","Experimental","Experimental","Experimental","Wildtype","Wildtype","Wildtype","Wildtype") #Treatments
 metadata<-as.data.frame(cbind(FileName,SampleName,Condition)) #create data frame 

 metadata$Condition <- as.factor(metadata$Condition) #Factor condition for DEseq
 
 count_data<-as.data.frame(counters)%>% #Renaming all column heads
   dplyr::rename(r1_S28 = r1_S28_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(r2_S29 = r2_S29_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(r3_S30 = r3_S30_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(r4_S31 = r4_S31_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(wt1_S24 = wt1_S24_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(wt2_S25 = wt2_S25_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(wt3_S26 = wt3_S26_R1_001.fastq.subread.BAM)%>%
   dplyr::rename(wt4_S27 = wt4_S27_R1_001.fastq.subread.BAM)
 View(count_data)#Views data to check if column headers really switched 



#Check out the data 
 rownames(metadata)
 View(metadata)
 
 colnames(count_data)
 all(rownames(metadata) == colnames(count_data))

 #Creates the design count matrix from future counts 
 dds<- DESeqDataSetFromMatrix(countData= count_data,
                               colData = metadata,
                               design = ~ Condition)
 
 #we can use plotCounts fxn to compare the normalized counts
 #between treated and control groups for our top 6 genes

 
 #6 genes from the gene cluster that have Gene ID's 
 plotCounts(dds, gene="JF50_RS12590", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12595", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12600", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12605", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12610", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12615", intgroup="Condition")
 
 
 #Top 6 genes by Pvalue 
 plotCounts(dds, gene="JF50_RS00810", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12675", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12720", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS06640", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS12700", intgroup="Condition")
 plotCounts(dds, gene="JF50_RS16100", intgroup="Condition")
 
 
 
 #guts
 dds <- DESeq(dds)
 res <- results(dds)
 
 #Not subsetting the reads 
 
 res_unsub<-res
 
 resdata_unsub <- merge(as.data.frame(res_unsub), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
 names(resdata_unsub)[1] <- 'gene'
 write.csv(resdata_unsub, file = paste0( "-results-with-normalized_unsub.csv"))
 
 #Subsetting P.adjusts by p<0.05
 
 res= subset(res, padj<0.05) # order results by padj value (most significant to least)
 res <- res[order(res$padj),] # should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
 
 
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
 names(resdata)[1] <- 'gene'
 write.csv(resdata, file = paste0( "-results-with-normalized_subset.csv"))
 
 
 
 
 # save data results and normalized reads to csv
 resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
 names(resdata)[1] <- 'gene'
 write.csv(resdata, file = paste0( "-results-with-normalized.csv"))
 
 
 
 
 # send normalized counts to tab delimited file for GSEA, etc.
 write.table(as.data.frame(counts(dds),normalized=T), file = paste0("_normalized_counts.txt"), sep = '\t')
 
 # produce DataFrame of results of statistical tests
 mcols(res, use.names = T)
 write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0("-test-conditions.csv"))
 
 # replacing outlier value with estimated value as predicted by distrubution using
 # "trimmed mean" approach. recommended if you have several replicates per treatment
 # DESeq2 will automatically do this if you have 7 or more replicates
 
 ddsClean <- replaceOutliersWithTrimmedMean(dds)
 ddsClean <- DESeq(ddsClean)
 tab <- table(initial = results(dds)$padj < 0.05,
              cleaned = results(ddsClean)$padj < 0.05)
 addmargins(tab)
 write.csv(as.data.frame(tab),file = paste0("-replaceoutliers.csv"))
 resClean <- results(ddsClean)
 resClean = subset(res, padj<0.05)
 resClean <- resClean[order(resClean$padj),]
 write.csv(as.data.frame(resClean),file = paste0("-replaceoutliers-results.csv"))
 
 ####################################################################################
 # Exploratory data analysis of RNAseq data with DESeq2
 #
 # these next R scripts are for a variety of visualization, QC and other plots to
 # get a sense of what the RNAseq data looks like based on DESEq2 analysis
 # 
 # 
 # 1) Differential heat expression maps 
 # 2) MA plot
 # 3) rlog stabilization and variance stabiliazation
 # 4) variance stabilization plot
 # 5) heatmap of clustering analysis
 # 6) PCA plot
 #
 #
 ####################################################################################
 
 #Differential expression gene heat map 

 #148 SIGNIFICANT differential expression genes 
 
 results_with_normalized <- read_csv("-results-with-normalized_subset.csv")
 View(results_with_normalized_subset) #Results that were only significant n=148
 
 results_with_normalized_unsub <- read_csv("-results-with-normalized_unsub.csv")
 View(results_with_normalized_unsub) # All Results 
 
 results_with_normalized_50<- head( results_with_normalized,50)
 View(results_with_normalized_50) #Top 50 Significat genes 
 
 results_with_normalized_20<- head( results_with_normalized,20)
 View(results_with_normalized_20) #Top 20 Significant genes 
 
 
#Top 50 genes data manipulation
 results_with_normalized_sum_50 <-   results_with_normalized_50 %>%
    select(c(gene,r1_S28,r2_S29,r3_S30,r4_S31,wt1_S24,wt2_S25,wt3_S26,wt4_S27))%>%
    pivot_longer(!gene, names_to = "Sample", values_to = "DE_values")%>%
    mutate(Log_DE = log(DE_values))%>%
    mutate(gene=factor(gene))
 

 #Top 20 genes data manipulation 
 results_with_normalized_sum_20 <-   results_with_normalized_20 %>%
   select(c(gene,r1_S28,r2_S29,r3_S30,r4_S31,wt1_S24,wt2_S25,wt3_S26,wt4_S27))%>%
   pivot_longer(!gene, names_to = "Sample", values_to = "DE_values")%>%
   mutate(Log_DE = log(DE_values))%>%
   mutate(gene=factor(gene))
 
 
 #Top 20 genes ggplot 
 
 ggplot(data =  results_with_normalized_sum_20, mapping = aes(x = Sample,
       y = gene, fill = Log_DE)) +
    scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"),breaks=seq(-2,10,1)) +
   xlab("Samples")+
   ylab("Gene")+
   ggtitle("20 Differentially Expressed Genes (P value < 0.05)", position = "center")+
    geom_tile(color = "white",
              lwd = 1.5,
              linetype = 1)+
    theme_linedraw()+
    theme(text=element_text(size=16,  family="Times New Roman"))
 
 ggsave("Figures/Top20DE.jpg")
 
 #Top 50 genes ggplot 
 ggplot(data =  results_with_normalized_sum_50, mapping = aes(x = Sample,
                                                           y = gene, fill = Log_DE)) +
   scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"),breaks=seq(-2,10,1)) +
   guides(color = guide_legend(title = "Users By guides"))+
   xlab("Samples")+
   ylab("Gene")+
   ggtitle("50 Differentially Expressed Genes (P value < 0.05)")+
   geom_tile(color = "white",
             lwd = 1.5,
             linetype = 1)+
    theme_linedraw()+
    theme(text=element_text(size=16,  family="Times New Roman"))
 
 ggsave("Figures/Top50DE.jpg")
 
 
 
 #Subsetted genes for P.luteo mac cluster 
 
 View(results_with_normalized_unsub)

  #Select rows 2532 to 2554
 
 p_luteo_cluster<-results_with_normalized_unsub[2532:2574, 1:16]
View( p_luteo_cluster)
 

#P.luteo cluster data manipulation 
p_luteo_cluster_results  <-   p_luteo_cluster %>%
   dplyr::filter(gene != "JF50_RS26085")%>%
   dplyr::filter(gene != "JF50_RS26080")%>%
   dplyr::filter(gene != "JF50_RS12570")%>%
   dplyr::select(c(gene,r1_S28,r2_S29,r3_S30,r4_S31,wt1_S24,wt2_S25,wt3_S26,wt4_S27))%>%
   pivot_longer(!gene, names_to = "Sample", values_to = "DE_values")%>%
   mutate(Log_DE = log(DE_values))%>%
   mutate(gene=factor(gene))


#P.luteo cluster ggplot 

ggplot(data =  p_luteo_cluster_results, mapping = aes(x = Sample,
                                                             y = gene, fill = Log_DE)) +
   scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"),breaks=seq(-2,10,1)) +
   xlab("Samples")+
   ylab("Gene")+
   ggtitle("P.Luteo MAC cluster")+
   geom_tile(color = "white",
             lwd = 1.5,
             linetype = 1)+
   theme_linedraw()+
   theme(text=element_text(size=16,  family="Times New Roman"))

ggsave("Figures/P.Luteo.jpg")


#SIGNIFICANT P.luteo cluster data manipulation 
p_luteo_cluster_results  <-   p_luteo_cluster %>%
   dplyr::filter(padj < 0.05)%>%
   dplyr::filter(gene != "JF50_RS26085")%>%
   dplyr::filter(gene != "JF50_RS12570")%>%
   dplyr::select(c(gene,r1_S28,r2_S29,r3_S30,r4_S31,wt1_S24,wt2_S25,wt3_S26,wt4_S27))%>%
   pivot_longer(!gene, names_to = "Sample", values_to = "DE_values")%>%
   mutate(Log_DE = log(DE_values))%>%
   mutate(gene=factor(gene))


#SIGNIFICANT P.luteo cluster ggplot 

ggplot(data =  p_luteo_cluster_results, mapping = aes(x = Sample,
                                                      y = gene, fill = Log_DE)) +
   scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"),breaks=seq(-2,10,1)) +
   xlab("Samples")+
   ylab("Gene")+
   ggtitle("P.Luteo MAC cluster (P<0.05)")+
   geom_tile(color = "white",
             lwd = 1.5,
             linetype = 1)+
   theme_linedraw()+
   theme(text=element_text(size=16,  family="Times New Roman"))

ggsave("Figures/P.Luteo_SIG.jpg")

 
 
 #we can use plotCounts fxn to compare the normalized counts
 #between treated and control groups for our top 6 genes
 par(mfrow=c(2,3))
 
 plotCounts(dds, gene="JF50_RS00810", intgroup="dex")
 plotCounts(dds, gene="ENSG00000179094", intgroup="condition")
 plotCounts(dds, gene="ENSG00000116584", intgroup="condition")
 plotCounts(dds, gene="ENSG00000189221", intgroup="condition")
 plotCounts(dds, gene="ENSG00000120129", intgroup="condition")
 plotCounts(dds, gene="ENSG00000148175", intgroup="condition")
 
 
 # MA plot of RNAseq data for entire dataset
 # http://en.wikipedia.org/wiki/MA_plot
 # genes with padj < 0.1 are colored Red
 plotMA(dds, ylim=c(-8,8),main = "Wildtype vs. Experimental")
 dev.copy(png, paste0("MAplot_initial_analysis.png"))
 dev.off()
 
 # transform raw counts into normalized values
 # DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
 # variance stabilization is very good for heatmaps, etc.
 rld <- rlogTransformation(dds, blind=T)
 vsd <- varianceStabilizingTransformation(dds, blind=T)
 
 
 # save normalized values
 write.table(as.data.frame(assay(rld),file = paste0("-rlog-transformed-counts.txt"), sep = '\t'))
 write.table(as.data.frame(assay(vsd),file = paste0("-vst-transformed-counts.txt"), sep = '\t'))
 
 # plot to show effect of transformation
 # axis is square root of variance over the mean for all samples
 par(mai = ifelse(1:4 <= 2, par('mai'),0))
 px <- counts(dds)[,1] / sizeFactors(dds)[1]
 ord <- order(px)
 ord <- ord[px[ord] < 150]
 ord <- ord[seq(1,length(ord),length=50)]
 last <- ord[length(ord)]
 vstcol <- c('blue','black')
 matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
 legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
 dev.copy(png,paste0("variance_stabilizing.png"))
 dev.off()
 
 # clustering analysis
 # excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
 library("RColorBrewer")
 library("gplots")
 distsRL <- dist(t(assay(rld)))
 mat <- as.matrix(distsRL)
 rownames(mat) <- colnames(mat) <- with(colData(dds), paste(Condition,sep=" : "))
 #Or if you want conditions use:
 #rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
 hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
 dev.copy(png, paste0("clustering.png"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
 dev.off()
 
 
 #Principal components plot shows additional but rough clustering of samples
 library("genefilter")
 library("ggplot2")
 library("grDevices")
 
 rv <- rowVars(assay(rld))
 select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
 pc <- prcomp(t(assay(vsd)[select,]))
 pc
 # set condition
 condition <- Condition
 scores <- data.frame(pc$x, condition)
 
 (pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
   + geom_point(size = 5)
   + ggtitle("Principal Components")
   + scale_colour_brewer(name = " ", palette = "Set1")
   + theme(
     plot.title = element_text(face = 'bold'),
     legend.position = c(.9,.2),
     legend.key = element_rect(fill = 'NA'),
     legend.text = element_text(size = 10, face = "bold"),
     axis.text.y = element_text(colour = "Black"),
     axis.text.x = element_text(colour = "Black"),
     axis.title.x = element_text(face = "bold"),
     axis.title.y = element_text(face = 'bold'),
     panel.grid.major.x = element_blank(),
     panel.grid.major.y = element_blank(),
     panel.grid.minor.x = element_blank(),
     panel.grid.minor.y = element_blank(),
     panel.background = element_rect(color = 'black',fill = NA)
   ))
 
 ggsave(pcaplot,file=paste0("PCA-ggplot2.jpg"))
 
 # scatter plot of rlog transformations between Sample conditions
 # nice way to compare control and experimental samples
 head(assay(rld))
 # plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
 plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
 plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
 plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
 plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
 
 
 
 # heatmap of data
 library("RColorBrewer")
 library("gplots")
 # 1000 top expressed genes with heatmap.2
 select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:1000]
 select
 my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
 heatmap.2(assay(vsd)[select,], col=my_palette,
           scale="row", key=T, keysize=1, symkey=T,
           density.info="none", trace="none",
           cexCol=0.6, labRow=F,
           main="1000 Top Expressed Genes Heatmap")
 dev.copy(png, paste0("1000_top_expressed_genes-HEATMAP.png"))
 dev.off()
 
 # 100 top expressed genes with heatmap.2
 select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:100]
 my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
 heatmap.2(assay(vsd)[select,], col=my_palette,
           scale="row", key=T, keysize=1, symkey=T,
           density.info="none", trace="none",
           cexCol=0.6, labRow=F,
           main="100 Top Expressed Genes Heatmap")
 dev.copy(png, paste0("100_top_expressed_genes-HEATMAP.png"))
 dev.off()
 
 
 # 20 top expressed genes with heatmap.2
 select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:20]
 my_palette <- colorRampPalette(c("blue",'white','red'))(n=20)
 heatmap.2(assay(vsd)[select,], col=my_palette,
           scale="row", key=T, keysize=1, symkey=T,
           density.info="none", trace="none",
           cexCol=0.6, labRow=F,
           main="100 Top Expressed Genes Heatmap")
 dev.copy(png, paste0("10_top_expressed_genes-HEATMAP.png"))
 dev.off()
 
 #Sizefactor and abline plots 
 
 dds <- estimateSizeFactors(dds)# Size factor estimation
 
 sizeFactors(dds) #Sizefactor estimation check point
 
 plot(sizeFactors(dds), colSums(counts(dds))) #Plot size factor (QA check for samples)
 abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0)) #Plot abline on the graph
 dev.copy(png, paste0("abline.png"))
 
 #Dendrogram cluster 
 
 normlzd_dds <- counts(dds, normalized=T) #normalize your counts 
 
 head(normlzd_dds) #check first five lines of code for normalized counts 
 
 plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$protocol) #plot normalized counts 
 
 
 plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$Condition) #Plot dendrogram of normalized reads against conditions
 
 plot(log(normlzd_dds[,1])+1, log(normlzd_dds[,2])+1, cex =.1) #Normalize dendrogram
 dev.copy(png, paste0("hclust.png"))
 #Sample Correlation Heat Map
 
 vsd <- vst(dds, blind = T)  # Varaiance Stabilizing transformation
 
 
 vsd_mat <- assay(vsd)  # extract the vst matris from the object
 
 
 vsd_cor <- cor(vsd_mat)  # compute pairwise correlation values
 
 vsd_cor #Check the call line 
 
 pheatmap(vsd_cor) #Sample heatmap comparing sample correlations 
 dev.copy(png, paste0("Correlation-HEATMAP.png"))
 
 ##DESeq2 model - Dispersion
 
 mean_readCounts <- apply(count_data[,1:3], 1, mean) #calculating mean for each gene
 
 
 var_readCounts <- apply(count_data[,1:3], 1, var)  # Calculating variance for each gene
 
 df <- data.frame(mean_readCounts, var_readCounts) #Turning mean and variance into a data frame (df)

 #ggplot 
 ggplot(df) + 
   geom_point(aes(x=mean_readCounts, y= var_readCounts)) + #ggplot of mean and variance 
   scale_y_log10() + #Convert y axis to log scale 
   scale_x_log10() + #convert x axis to log scale 
   xlab("Mean counts per gene") + #label x axis 
   ylab("Variance per gene") + #label y axis 
   labs(title = "DESeq2 model - Dispersion") #Title the graph 
 
 ggsave("Figures/DEseq2Model.jpg")
 