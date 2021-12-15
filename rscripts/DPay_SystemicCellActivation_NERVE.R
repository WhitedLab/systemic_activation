# HEADER ---------------------------- #
# Contributor: Steven Blair
# Email:  steven_blair@fas.harvard.edu
# Date: 2021-12-03
# Script Name: DPay_SystemicCellActivation_NERVE.R
# Script Description: DGE analysis of denervated vs sham, sham representing baseline.
# i. PREPARE ENVIRONMENT --------------
# Clear and previous data or values as well as closing any active devices.
rm(list=ls())
dev.off()
# SET WORKING DIRECTORY ------------- #
setwd("/Users/sblair/Documents/DPay_SystemicCellActivation_20211201")
# SET OPTIONS ----------------------- #
# ii. INSTALL PACKAGES & LOAD LIBRARIES
# Libraries for creating DDS
# BiocManager::install("rhdf5")
library(tximport) # BiocManager::install("tximport")
library(DESeq2) # BiocManager::install("DESeq2")
library(tibble) # install.packages("tibble")
# For annotation
library(dplyr) # devtools::install_github("hadley/dplyr")
library(readr) # install.packages("readr")
library(tidyverse) # install.packages("tidyverse")
library(data.table)
library(pheatmap)
library(viridis)

# 1) TXI IMPORT -----------------------
# Before we begin the analysis we must import the count data created through kallisto.
# We first create a variable for the location of where the counts are stored on 
# the machine in relation to the working directory
dir <- "./kallisto_output"
# a simple samplesheet was created for the following samples
# samples:
# G01v1_7_S6 (dfn)  denervated forelimb
# D02v1_12_S11 (dfn)  denervated forelimb
# G02v1_15_S14 (dfn)  denervated forelimb
# G03v1_23_S22 (dfn)  denervated forelimb
# B04v2_26_S25 (sfn)  sham forelimb
# C06v1_43_S41 (sfn)  sham forelimb
# A05v1_33_S32 (sfn)  sham forelimb
# this was saved as a tsv file and loaded with the following command
samplesheet <- read.table("samplesheet.tsv", header = T, sep = '\t', row.names = NULL)

# note the column "ident" has the full sample name and is similar to the 
# corresponding filename.  The column sample_type identifies if the sample 
# is denervated or sham

# txi will be used to create a summary matrix of transcript lengths, and estimated counts
# for later downstream analysis. we can simply use the data in the samplesheet 
# to load the counts and other data in from the file location.
# We create a list of locations for where the counts are by file.path
# this effectively concatenate the file locations by the directory "dir", 
# subdirectory that the counts are in "counts", the samplename which is 
# part of the filename and append the rest of the file names which in this case is "_counts.genes.results".  
files <- file.path(dir, paste0(samplesheet$ident,"/","abundance.h5"))
# Now we identify the samples by the sample_type which was saved in the annotation file samplesheet.tsv
names(files) <- paste0(samplesheet$sample_type[1:7])

# The files Values is an array of lists. Each list has a sampletype 
# (in this case sham or denervated) and a location for the sample.  This is 
# important for later to keep track of the samples by read id, condition, well id, 
# and others, all through the saved samplesheet.

# Now to use the bioconductor tool tximport to create a matrix of counts by sampletype
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
# this has created al arge list with 5 elements
# clean up the environment of uneeded values
rm(dir, files)

# 2) CREATE THE DESeqDataSet ----------
# Import the transcript quantifications (txi) from kallisto into DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, samplesheet, ~ operation)

# Now we can run DESeq2 on this primed dataset using DESeq2 command
dds <- DESeq(dds)

# save the Deseq2 DataSet
saveRDS(dds,file = "DPay_SystemicCellActivation_nerve.rds")

# Extract a results table using sham as the reference level
res <- results(dds, contrast=c("operation","denervated","sham"))

# print all results with a pvalue
results <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% # tibble library
  as_tibble()
results <- results %>% filter(!is.na(pvalue))
# write.table(results, file='DPay_SystemicActivation_NerveAllResults.xls', quote=FALSE, sep='\t', row.names = T)
# Seperate the significant results using an adjusted p of 0.05 and a lfc of 0.58 (log1.5 base 2)
rm(results)
sig_genes <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% # tibble library
  filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
  as_tibble()

# Create a dataframe of the counts
counts <- data.frame(counts(dds))

# Clean up
rm(samplesheet, txi.kallisto)
# 3) ANNOTATION -----------------------

# The Bryant transcriptome is annotated by blast results and uniprot data.
# DESeq2 results are labeled by contig, using the contigs from the 
# DESeq2 results (res) allow us to gather these annotations.
# We change the label of column to contig, to allow the 
# matrix to be merged using these identifiers.

names(sig_genes)[1] <- "contig"

# load the provided putative annotation
anno <- read_tsv("DPay_SystemicCellActivationNerve_Annotation.txt")

# Join the DESeq2 results with the annotation
sig_genes <- left_join(sig_genes, anno)

# tidy the result matrix (remove unannotated contigs for the purpose of visualization)
sig_genes <- sig_genes %>% filter(!is.na(gene))

# later we will need to remember the gene name in reference to the contig
con2g <- sig_genes[,-2:-9]

# tidy up
rm(anno)
# 4) VISUALIZATION --------------------

# If starting here you will want to load the following data
dds <-readRDS(file = "DPay_SystemicCellActivation_nerve.rds")

# perform rlog transformation using the design formula (blind = F) and make the table reference
rld <- rlog(dds, blind = F)
rld_assay <- assay(rld)
rld_assay <- data.frame(rld_assay)
rld_assay$contig <- row.names(rld_assay)
# now subset these rlog transformed results with the sig_genes from before
rld_sig_counts <- subset(rld_assay, rld_assay$contig %in% sig_genes$contig)

# scale the gene data, we need to transpose the matrix to do this
rld_sig_gene <- t(rld_sig_counts[1:7])
rld_sig_gene <- scale(rld_sig_gene)
rld_sig_gene <- t(rld_sig_gene)
# n <- head(rld_sig_gene, n = 248L)

# Let us merge this list with the deseq2 results
sig_gene_results <- left_join(sig_genes, rld_sig_counts)
                              
# rename the rows using the putative gene names
sig_gene_results <- left_join(sig_gene_results, con2g)
rownames(sig_gene_results) <- make.names(sig_gene_results$gene, unique = TRUE)

# save the complete results matrix
write.table(sig_gene_results, file='DPay_SystemicCellActivation_NerveResults.xls', quote=FALSE, sep='\t', row.names = F)

# create a list of gene names that will be used to label the rows (italics)
heatmap_names <- lapply(
  rownames(sig_gene_results),
  function(x) bquote(italic(.(x))))

#draw heatmap allowing larger margins and adjusting row label font size
heatmap <- pheatmap(sig_gene_results[11:17], 
         scale="row",
         cluster_rows=T,
         cluster_cols=T,
         show_rownames = T,
         main = "Nerve Denervated vs. Sham DGE (row-based z-score)",
         color = viridis(16,direction = -1),
         border_color = NA,
         labels_row = as.expression(heatmap_names),
         treeheight_row = 25,
         treeheight_col =  15)
dev.off()

# save the heatmap in multiple formats
svg("DPay_heatmap_allNerve.svg",width=7, height=38)
heatmap
dev.off()

png("DPay_heatmap_allNerve.png",width=672, height=3456)
heatmap
dev.off()

pdf("DPay_heatmap_allNerve.pdf",width=7, height=38)
heatmap
dev.off()

rm(heatmap_names, heatmap)
# WE will create a subset of 50 genes for a smaller heatmap.  
# grabbing 50 rows, 25 up, 25 down, including the two genes Adra2a and lef1.
subset <- c("Col11a2","Creg2","Fbln5","Lmod1","Cyp20a1","Mfap4","Add3","Sptbn1","Synrg","Pofut1","Gtf3c3","Fbln5.1","Kat7","Lvrn","Cldn9","Bccip","Kcnmb2","Avd","Dpp4","Snca","Drp2","C4b","Nsl1","Magi1","Cryab","Adamts17","Thbs4","Crabp1","Kazald2","Sall3","Srpx2","Vgf","Tmpo","Kcp","Gpam","Sestd1","Shc4","Mmp13","Znf236","Lef1","Fosl2","Ralgds","Rhbdl2","Cd276","Swap70","Sepsecs","Cpm","Il11","Adra2a","Cyp27b1")
subset <- subset(sig_gene_results, sig_gene_results$gene %in% subset)

heatmap_names <- lapply(
  rownames(subset),
  function(x) bquote(italic(.(x))))

heatmap <- pheatmap(subset[11:17], 
                    scale="row",
                    cluster_rows=T,
                    cluster_cols=T,
                    show_rownames = T,
                    main = "Nerve Denervated vs. Sham DGE (row-based z-score)",
                    color = viridis(16,direction = -1),
                    border_color = NA,
                    labels_row = as.expression(heatmap_names),
                    treeheight_row = 25,
                    treeheight_col =  15)
dev.off()
svg("DPay_heatmap_nerveSubset.svg",width=7, height=9)
heatmap
dev.off()
pdf("DPay_heatmap_nerveSubset.pdf",width=7, height=9)
heatmap
dev.off()
png("DPay_heatmap_nerveSubset.png",width=672, height=710)
heatmap
dev.off()
