# RUVseq data normalizartion
#https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html
# Analyzing RNA-seq data with DESeq2
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# load libraries * --------------------------------------------------------
library(RUVSeq)
library("tximeta")
library(DESeq2)
library(tximport)
library(RColorBrewer)
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("dplyr")
library("ggplot2")

# Load data * ------------------------------------------------------------
dir = #directory containing files
list.files(dir)

file = file.path(dir, "SampleInfo.csv") #contains sampleName,fileName, Condition and Replicate 
file.exists(file)
sampleTable <- read.csv(file, stringsAsFactors=FALSE)
sampleTable$Condition <- factor(sampleTable$Condition)

# build the DESeqDataSet * --------------------------------------------------

quants = file.path(dir, "/Quants") # directory that contains read counts files 
list.files(quants)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = quants,
                                  design= ~ Condition)
dds

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
#dds <- dds[rowSums(counts(dds)) < 1000000, ]

dds <- estimateSizeFactors(dds)

filtered <- counts(dds)
genes <- rownames(counts(dds))

#x <- as.factor(rep(c("Mc", "Ec"), each = 4))
#x = as.factor(c("Mc", "Mc", "Mc", "Mc", "Mm", "Mm", "Mm", "Mm","Ec", "Ec", "Ec","Em", "Em", "Em"))
x = as.factor(c("Ec", "Ec", "Ec","Em", "Em", "Em"))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set

# Raw data exploration *-------------------------------------------------
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=TRUE, ylim=c(-2.1,2.1), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2,  ylim=c(-1,1), xlim=c(-0.8,0.8))

# Upper-quartile (UQ) normalization * ---------------------------------------
# adjust for sequencing depth 
set <- betweenLaneNormalization(set, which="upper")
# The boxplots of relative log expression (RLE = log-ratio of read count to median read count across sample) and plots of principal components (PC) reveal a clear need for betwen-sample normalization.
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=TRUE, ylim=c(-2.1,2.1), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2,  ylim=c(-1,1), xlim=c(-0.8,0.8))

# RUVs normalization *------------------------------------------------------
## Estimating the factors of unwanted variation using replicate samples
# e.g., batch, library preparation, etc.
differences <- makeGroups(x)
differences
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)
colors <- brewer.pal(3, "Set2")
plotRLE(set3, outline=FALSE, ylim=c(-0.5,0.5), col=colors[x])
plotPCA(set3, col=colors[x], cex=1.2, ylim=c(-1,1), xlim=c(-0.8,0.8))


# Exploración de calidad datos ----------------------------------------------------
#LOVE tutorial
normalized_counts  <- normCounts(set3)
log.norm.counts <- log2(normalized_counts + 1)
rlog.norm.counts <- rlogTransformation(normalized_counts, blind = FALSE)

# Effects of transformations on the variance 
## log2 transform
meanSdPlot(log.norm.counts)

## rlogtranform, with variance stabilization
meanSdPlot(rlog.norm.counts)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(rlog.norm.counts))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sampleNames(set3)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Heatmap of sample-to-sample distances using the Poisson Distance.
library("PoiClaClu")
poisd <- PoissonDistance(t(rlog.norm.counts ))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- sampleNames(set3)
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# dendogram
# cor () calculates the correlation between columns of a matrix
distance.m_rlog <- as.dist (1 - cor(rlog.norm.counts , method = "pearson"))

# plot () can directly interpret the output of hclust ()
plot(hclust(distance.m_rlog, method = "complete"),
     labels = colnames(rlog.norm.counts),
     main = "rlog transformed read counts\ndistance:Pearson correlation")

# heat map per gene (20 genes)
select <- order(rowMeans(normalized_counts),
                decreasing=TRUE)[1:20]

pheatmap(rlog.norm.counts[select,], cluster_rows=FALSE, show_rownames=TRUE,
        cluster_cols=FALSE) #annotation_col=df 


# Differential expression analysis Wald test of significance *--------------
dds = DESeqDataSetFromMatrix(countData = counts(set3),
                             colData = pData(set3),
                             design = ~ W_1 + x)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds = DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
res = results(dds,  alpha=0.05)
res
summary(res)
resultsNames(dds)

# Export data *----------------------------------------------------
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# AllData
resOrdered <- res[order(res$pvalue),] #ordenar por pvalue
resOrderedDF <- as.data.frame(resOrdered)
#write.csv(resOrderedDF, file = "ResultsDEanalysis.csv")

# only the results which pass an adjusted p value threshold 
resSig <- subset(resOrdered, padj < 0.05)
write.csv(resSig, file = "ResultsDEanalysis_SigGenes.csv")

# Session information *------------------------------------------------
sessionInfo()
