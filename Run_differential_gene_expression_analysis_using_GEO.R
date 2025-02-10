## Study: Rare Variant Analyses in Ancestrally Diverse Cohorts Reveal Novel ADHD Risk Genes
## Analysis: Differential expression analysis
## Author: Seulgi Jung


BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("umap")


# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

library(dplyr)
library(XML)
library(edgeR)
library(org.Hs.eg.db)
library(annotate)
library(GO.db)

setwd("C:\\Users\\jsg79\\Dropbox\\MSSM\\Seulgi\\ADHD\\GEO") ## set work space


# load series and platform data from GEO
gset <- getGEO("GSE159104", GSEMatrix =TRUE, AnnotGPL=TRUE)
pheno <- read.delim("GSE159104_samples.txt.gz", sep = "\t", header = TRUE)
head(pheno)
write.csv(pheno, "GSE159104_samples.csv", row.names = TRUE)


# Check supplementary file links
sup_files <- experimentData(gset[[1]])@other$supplementary_file
print(sup_files)

# Define the URL for the supplementary file
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159104/suppl/GSE159104_DGE_raw_154samples.txt.gz"

# Download the file
download.file(url, destfile = "GSE159104_DGE_raw_154samples.txt.gz")


# Read the file into R
expression_data <- read.delim("GSE159104_DGE_raw_154samples.txt.gz", sep = "\t", header = TRUE)

# Convert to DataFrame
expression_df <- as.data.frame(expression_data)

# Inspect the first few rows
head(expression_df)

# Optionally save to a CSV file
write.csv(expression_df, "GSE159104_expression_data.csv", row.names = FALSE)


input <- expression_df %>%
  dplyr::select(Gene.name, P14.2016.03.31T12_40_25.0220160771004.fc1.ch8, P14.2016.03.31T12_40_25.0220160771004.fc1.ch10, P14.2016.03.31T12_40_25.0220160771004.fc1.ch12, P14.2016.03.31T12_40_25.0220160771004.fc1.ch19, P14.2016.03.31T12_40_25.0220160771004.fc1.ch23, P14.2016.03.31T12_40_25.0220160771004.fc1.ch25, P14.2016.03.31T12_40_25.0220160771005.fc2.ch1, P14.2016.03.31T12_40_25.0220160771005.fc2.ch2, P14.2016.03.31T12_40_25.0220160771005.fc2.ch4, P14.2016.03.31T12_40_25.0220160771005.fc2.ch9, P14.2016.03.31T12_40_25.0220160771005.fc2.ch12, P14.2016.03.31T12_40_25.0220160771005.fc2.ch18, P14.2016.03.31T12_40_25.0220160771005.fc2.ch19, P13.2016.04.08T16_18_11.0220160891003.fc1.ch1, P13.2016.04.08T16_18_11.0220160891003.fc1.ch3, P13.2016.04.08T16_18_11.0220160891003.fc1.ch5, P13.2016.04.08T16_18_11.0220160891003.fc1.ch19, P13.2016.04.08T16_18_11.0220160891003.fc1.ch20, P14.2016.04.15T16_10_44.0220160891004.fc1.ch22, P14.2016.04.15T16_10_44.0220160891005.fc2.ch1, P14.2016.04.15T16_10_44.0220160891005.fc2.ch20, P14.2016.04.15T16_10_44.0220160891005.fc2.ch24, P13.2016.04.22T14_16_47.0220160891007.fc1.ch14, P14.2016.03.31T12_40_25.0220160771004.fc1.ch7, P14.2016.03.31T12_40_25.0220160771004.fc1.ch17, P14.2016.03.31T12_40_25.0220160771004.fc1.ch24, P14.2016.03.31T12_40_25.0220160771005.fc2.ch3, P14.2016.03.31T12_40_25.0220160771005.fc2.ch5, P14.2016.03.31T12_40_25.0220160771005.fc2.ch6, P14.2016.03.31T12_40_25.0220160771005.fc2.ch11, P14.2016.03.31T12_40_25.0220160771005.fc2.ch15, P14.2016.03.31T12_40_25.0220160771005.fc2.ch20, P14.2016.03.31T12_40_25.0220160771005.fc2.ch23, P14.2016.03.31T12_40_25.0220160771005.fc2.ch25, P13.2016.04.08T16_18_11.0220160891003.fc1.ch2, P13.2016.04.08T16_18_11.0220160891003.fc1.ch4, P13.2016.04.08T16_18_11.0220160891003.fc1.ch10, P14.2016.04.15T16_10_44.0220160891004.fc1.ch20, P14.2016.04.15T16_10_44.0220160891004.fc1.ch24, P14.2016.04.15T16_10_44.0220160891005.fc2.ch14, P13.2016.04.22T14_16_47.0220160891007.fc1.ch9, P13.2016.04.22T14_16_47.0220160891008.fc2.ch9, P13.2016.05.31T15_11_07.1231461441002.fc1.ch1, P14.2016.06.07T15_56_34.0231461441003.fc1.ch19)

dim(input)

## generate a matrix of gene expression
y <- DGEList(counts=input[,2:45], genes=input[,1])

idfound <- y$genes$genes %in% mappedRkeys(org.Hs.egALIAS2EG)
y <- y[idfound,]
dim(y)
head(y)

egREFSEQ <- toTable(org.Hs.egALIAS2EG)
head(egREFSEQ)
m <- match(y$genes$genes, egREFSEQ$alias_symbol)
head(m)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
head(y$genes)

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$genes)
y <- y[!d,]
e <- duplicated(y$genes$EntrezGene)
y <- y[!e,]
nrow(y)

## check a MDS plot
y$samples$lib.size <- colSums(y$counts)
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL
y <- calcNormFactors(y)
y$samples
plotMDS(y) ## checking the outliers
dim(y)


# group membership for all samples
## case-control information is here: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE90598

## allocate case-control phenotypes
gsms <- paste0("11111111111111111111111000000000000000000000")
sml <- strsplit(gsms, split="")[[1]]
group <- factor(sml)

data.frame(Sample=colnames(y),group)

design <- model.matrix(~group)
rownames(design) <- colnames(y)
design


## MDS plot
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, contrast=c(1, -1))
lfc <- coef(fit)
plotMDS(y, labels=colnames(y))


## perform the differential expression analysis
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
dim(qlf) ## 22779
topTags(qlf, n=22779) -> p
write.table(p, "GSE90598_result.out", col.names=T, row.names = F, quote=F, sep="\t")
result <- read.table("GSE90598_result.out", header=T)
head(result)
result[which((result$logFC > 1.5) & (result$FDR < 0.05)),]


colnames(design)
o <- order(qlf$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(qlf))
plotMD(qlf) ## volcano plot
abline(h=c(-1, 1), col="green")

go <- goana(qlf)
topGO(go, ont="BP", sort="Up", n=5000)
write.table(go, "GSE90598_result_GO.out", col.names=T, row.names = F, quote=F, sep="\t")

