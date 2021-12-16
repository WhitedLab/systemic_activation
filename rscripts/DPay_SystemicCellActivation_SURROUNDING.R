# HEADER ---------------------------- #
# Contributor: Steven Blair
# Email:  steven_blair@fas.harvard.edu
# Date: 2021-12-03
# Script Name: DPay_SystemicCellActivation_SURROUNDING.R
# Script Description: DGE analysis of denervated vs sham, sham representing baseline in surrounding tissue.
# i. PREPARE ENVIRONMENT --------------
# Clear and previous data or values as well as closing any active devices.
rm(list=ls())
dev.off()
# SET WORKING DIRECTORY ------------- #
setwd("/Users/sblair/Documents/DPay_SystemicCellActivation_surroundingTissue")
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
# For visualization
library(pheatmap)
library(viridis)

# 1) TXI IMPORT -----------------------
# Before we begin the analysis we must import the count data created through kallisto.
# We first create a variable for the location of where the counts are stored on 
# the machine in relation to the working directory
dir <- "./kallisto_output"
# a simple samplesheet was created for the following samples
# samples:
# A02v1_9_S8  denervated	forelimb	tissue
# A03v1_17_S16  denervated	forelimb	tissue
# A04v1_25_S24  denervated	forelimb	tissue
# B01v1_2_S2  denervated	forelimb	tissue
# B02v1_10_S9 denervated	forelimb	tissue
# D03v1_20_S19  denervated	forelimb	tissue
# E03v1_21_S20  denervated	forelimb	tissue
# C06v1_41_S39  sham	forelimb	tissue
# E05v1_37_S35  sham	forelimb	tissue
# E06v2_45_S43  sham	forelimb	tissue
# F04v1_30_S29  sham	forelimb	tissue
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
names(files) <- paste0(samplesheet$sample_type[1:11])

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

# Now we can estimate size factors and gene-wise dispersion estimates using 
# DESeq2 on this primed dataset using DESeq2 command
dds <- DESeq(dds)

# save the Deseq2 DataSet
saveRDS(dds,file = "DPay_SystemicCellActivation_surroundingDDS.rds")

# Extract a results table using sham as the reference level
res <- results(dds, contrast=c("operation","denervated","sham"))

# print all results with a pvalue
results <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% # tibble library
  as_tibble()
results <- results %>% filter(!is.na(pvalue))
write.table(results, file='DPay_SystemicActivation_SurroundingAllResults.xls', quote=FALSE, sep='\t', row.names = T)
# Seperate the significant results using an adjusted p of 0.05 and a lfc of 0.58 (log1.5 base 2)
rm(results)
sig_genes <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% # tibble library
  filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
  as_tibble()

# Create a dataframe of the counts
counts <- data.frame(counts(dds))

# save the counts
write.table(counts, file="DPay_SystemicCellActivation_surroundingCounts.xls", quote=FALSE, sep='\t', row.names = F)

# Clean up
rm(samplesheet, txi.kallisto)

# 3) ANNOTATION -----------------------

# The Bryant transcriptome is annotated by blast results and uniprot data.
# DESeq2 results are labeled by contig, using the contigs from the 
# DESeq2 results (res) allow us to gather these annotations.
# We change the label of column to contig, to allow the 
# matrix to be merged using these identifiers.

names(sig_genes)[1] <- "contig"

tmp <- subset(anno, anno$contig %in% sig_genes$contig)
tmp <- tmp %>% distinct()
# load the provided putative annotation
anno <- read_tsv("DPay_SystemicCellActivation_surroundingAnnotation.txt")

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
dds <-readRDS(file = "DPay_SystemicCellActivation_surroundingDDS.rds")

# perform rlog transformation using the design formula (blind = F) and make the table reference
rld <- rlog(dds, blind = F)
rld_assay <- assay(rld)
rld_assay <- data.frame(rld_assay)
rld_assay$contig <- row.names(rld_assay)
# now subset these rlog transformed results with the sig_genes from before
rld_sig_counts <- subset(rld_assay, rld_assay$contig %in% sig_genes$contig)

# scale the gene data, we need to transpose the matrix to do this
rld_sig_gene <- t(rld_sig_counts[1:11])
rld_sig_gene <- scale(rld_sig_gene)
rld_sig_gene <- t(rld_sig_gene)

# Let us merge this list with the deseq2 results
sig_gene_results <- left_join(sig_genes, rld_sig_counts)

# rename the rows using the putative gene names
sig_gene_results <- left_join(sig_gene_results, con2g)
rownames(sig_gene_results) <- make.names(sig_gene_results$gene, unique = TRUE)

# save the complete results matrix
write.table(sig_gene_results, file='DPay_SystemicCellActivation_surroundingResults.xls', quote=FALSE, sep='\t', row.names = F)

# create a list of gene names that will be used to label the rows (italics)
heatmap_names <- lapply(
  rownames(sig_gene_results),
  function(x) bquote(italic(.(x))))

#draw heatmap allowing larger margins and adjusting row label font size
heatmap <- pheatmap(sig_gene_results[11:21], 
                    scale="row",
                    cluster_rows=T,
                    cluster_cols=T,
                    show_rownames = T,
                    main = "Surrounding Tissue Denervated vs. Sham DGE (row-based z-score)",
                    color = viridis(16,direction = -1),
                    border_color = NA,
                    labels_row = as.expression(heatmap_names),
                    treeheight_row = 25,
                    treeheight_col =  15)
dev.off()

# save the heatmap in multiple formats
svg("DPay_heatmap_surroundingTissue.svg", width=7, height=12)
heatmap
dev.off()

png("DPay_heatmap_surroundingTissue.png", width=672, height=700)
heatmap
dev.off()

pdf("DPay_heatmap_surroundingTissue.pdf", width=7, height=12)
heatmap
dev.off()

rm(heatmap_names, heatmap)
# WE will create a subset of 50 genes for a smaller heatmap.  
# grabbing 50 rows, 25 up, 25 down, including the two genes Adra2a and lef1.
subset <- c("Itga6","Rpp40","Pcna","Haus5","Wbp1l","Apt","Ryr3","Pigm","Pex10","Septin2","Rab11a","Mkrn1","Xpot","Rad54l","Capn8","Rtase","Lats1","Abhd12","Ttf1","Dtx3l","Macrod1","Atf7","Lgals8","Gng2","Abracl","Srsf10","Lclat1","Uba2","Cyp3a5","Tmem104","Sbf1","Recql4","Mmp13","Aebp1","Ano9","Col28a1","Dcn","Pou2f1","Dpp9","Smap1","Ppp1r9a","Csde1","Tp63","Mycbp2","Strn3","Nos2","Igf1r","Sbf1.1","Irf3")
subset <- subset(sig_gene_results, sig_gene_results$gene %in% subset)

heatmap_names <- lapply(
  subset$gene,
  function(x) bquote(italic(.(x))))

heatmap <- pheatmap(subset[11:21], 
                    scale="row",
                    cluster_rows=T,
                    cluster_cols=T,
                    show_rownames = T,
                    main = "Surrounding Tissue Denervated vs. Sham DGE (row-based z-score)",
                    color = viridis(16,direction = -1),
                    border_color = NA,
                    labels_row = as.expression(heatmap_names),
                    treeheight_row = 25,
                    treeheight_col =  15)
dev.off()
svg("DPay_heatmap_surroundingSubset.svg",width=7, height=9)
heatmap
dev.off()
pdf("DPay_heatmap_surroundingSubset.pdf",width=7, height=9)
heatmap
dev.off()
png("DPay_heatmap_surroundingSubset.png",width=672, height=710)
heatmap
dev.off()

save.image(file = "DPay_SystemicCellActivation_surrounding.RData")
