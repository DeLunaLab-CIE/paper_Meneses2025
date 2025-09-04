                    # GO enrichment analysis #
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
################################################################################
# The GO term analysis was performed using differentially expressed genes obtained 
# with a cutoff of p > 0.1

# load libraries 
# BiocManager::install("clusterProfiler")
#BiocManager::install("sos")
library(clusterProfiler)
library(org.Sc.sgd.db)
library(sos)
library(stringr)
library("AnnotationHub")
library(ggplot2)

# Load Data ---------------------------------------------------------------
# All detected genes in the transcriptome. CSV file contain sgd id, systemic name and standard name 
fileAll <- #directory/file.csv
file.exists(fileAll)
Bkground = read.csv(fileAll, row.names=1, stringsAsFactors=FALSE)

# Upregulated genes. CSV file contain sgd id, systemic name and standard name
fileUp <- #directory/file.csv
file.exists(fileUp)
UP = read.csv(fileUp, row.names=1, stringsAsFactors=FALSE) 

# file containing downregulated genes. CSV file contain sgd id, systemic name and standard name
fileDown <- #directory/file.csv
file.exists(fileDown)
Down= read.csv(fileDown, row.names=1, stringsAsFactors=FALSE) #


# Create annotation -------------------------------------------------------
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Saccharomyces cerevisiae")[[1]]
orgdb

# enrichGO ----------------------------------------------------------------
### GO Enrichment Analysis of a gene set

UP = str_trim(rownames(UP))
Down = str_trim(rownames(Down))

gene = UP #UP (upregulates) or DOWN (downregulated) 
geneList = str_trim(rownames(Bkground))
#org.Sc.sgd.db
#keytypes(orgdb)
EGO <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = orgdb,
                keyType       = "SGD",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 5,
                maxGSSize     = 300, 
                readable      = FALSE)
head(EGO)
EGO[,2]
nrow(EGO) #numero de go enriquecidos

# Reduce GO term redundancy
SimEGO = simplify(EGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                   semData = NULL)
nrow(SimEGO)
SimEGO[,2]


# resulttable
tab = SimEGO@result
namefile = "EGO_BP_UP"; #cambiar nombre de archivo
write.table(tab, file = paste0(namefile, ".txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)
# plot
goplot(EGO)


# session information
sessionInfo()
