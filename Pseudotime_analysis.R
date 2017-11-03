#source("http://bioconductor.org/biocLite.R")
#biocLite("TSCAN")
#biocLite("M3Drop")
#biocLite("monocle")
#biocLite("destiny")
#biocLite("SLICER")
#biocLite("scater")


library(TSCAN)
library(M3Drop)
library(monocle)
library(destiny)
library(SLICER)

#install_github("davismcc/scater")
#biocLite("SingleCellExperiment")
#biocLite("scater")
#biocLite("knitr")
#biocLite("SingleCellExperiment")

library(SingleCellExperiment)
library(scater)
library(knitr)
options(stringsAsFactors = FALSE)

setwd("/media/vimal/10b357a2-8127-40d9-ab03-3e38ea800ce2/SC/scRNA.seq.course/")
#molecules <- read.table("tung/molecules.txt", sep = "\t")
#anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
#head(molecules[ , 1:3])

###  4 Biological Analysis

#biocLite("pcaMethods")
#install_github("JustinaZ/pcaReduce")
#install_github("drisso/SingleCellExperiment")
#biocLite("SC3")
#biocLite("pheatmap")
#biocLite("mclust")

library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
library(DESeq2)

set.seed(1234567)

deng <- readRDS("deng/deng-reads.rds")
table(colData(deng)$cell_type2)
scater::plotPCA(deng, colour_by = "cell_type2")


deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "cell_type1")

