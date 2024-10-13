## Meta
##Date: September 2024
##Author: Emily
##Description: Download TCGA/cbioportal data and look at differential expression analysis between samples with different copy number status of SULT1A1 

## Data download
##Data was downloaded through cbioportal (https://www.cbioportal.org/datasets) there is also a github page (https://github.com/cBioPortal/datahub/tree/master) that should be useable, although it wasn't working today.


install.packages("BiocManager")
BiocManager::install("DESeq2")

## Data analysis
##```r
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

## Set working directory 
setwd("//student.otago.ac.nz/mds/Profiles-V2/y/youem737/RedirectedFolders/Documents/SULT1A1")

rna_rsem <- read.delim("brca_tcga_pan_can_atlas_2018data_mrna_seq_v2_rsem.txt")
cna <- read.delim("brca_tcga_pan_can_atlas_2018data_cna.txt")

# Check to see if the samples are in both data sets - also use on entrezid (check for dups)
gene_meta <- rna_rsem[,1:2]
table(duplicated(gene_meta$Entrez_Gene_Id))
## 26 duplicates - will need to handle these some how


rna_samples <- colnames(rna_rsem)[3:ncol(rna_rsem)]
cna_samples <- colnames(cna)[3:ncol(cna)]

table(rna_samples %in% cna_samples)
## 1068/1082 rna sample have cna data (14 missing)
table(cna_samples %in% rna_samples)
## 1068/1070 cna sample have cna data (2 missing)

# drop samples with incomplete data (e.g. missing RNA or CNA data) - check samples are in the same order
rna_rsem <- rna_rsem[,colnames(rna_rsem) %in% colnames(cna)]
cna <- cna[,colnames(cna) %in% colnames(rna_rsem)]
all.equal(colnames(cna), colnames(rna_rsem))


## Generate sample info matrix for CNV status for a give GoI (SULT1A1)
GoI <- "SULT1A1"
cna_status <- cna %>% filter(Hugo_Symbol == GoI) %>% select(-c(1:2)) %>% t() %>% as.data.frame() %>% rename(CNV =1) %>% 
	mutate(CNV =case_when(
		.[[1]] == 0 ~ "Diploid",
		.[[1]] %in% c(1,2) ~ "Duplication",
		.[[1]] %in% c(-1,-2) ~ "Deletion"))


## RSEM estimated counts need to be rounded to whole intergers for DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(rna_rsem[,-c(1:2)]),
				colData = cna_status,
				design = ~ CNV
)


### Not neccessary for deseq but can be useful to remove lowly expressed genes - many ways to this, this method is recommended by deseq
smallestGroupSize <- min(table(cna_status$CNV))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## move all duplicates up here 
filtered_gene <- gene_meta[duplicated(gene_meta$Hugo_Symbol),]
filtered_gene <- gene_meta[duplicated(gene_meta$Entrez_Gene_Id),]

# Check to see how many hugo symbol duplicates are 
table(duplicated(gene_meta$Hugo_Symbol))
## 19 duplicates

## Filter Duplicated Hugo Symbol 
x <- gene_meta %>% filter(Hugo_Symbol != "") %>% filter(duplicated(Hugo_Symbol))
x1<- rna_rsem %>% filter(Hugo_Symbol %in% x$Hugo_Symbol)

## Generate expression level of the duplicated genes 
data.frame(x1$Hugo_Symbol, mean_exprs=apply(x1[,-c(1:2)],1, mean))

## Expression level of duplicated genes with both Hugo Symbol and Gene Id                                                                                                 
data.frame(x1$Hugo_Symbol,x1$Entrez_Gene_Id, mean_exprs=apply(x1[,-c(1:2)],1, mean))

## Filter out low expressed genes and see how many duplicates remain:dds right expression level from above filtering

#Duplicated Hugo Symbols in dds (where no lowly expressed genes are) data
duplicated_genes <- dds$Hugo_Symbol[duplicated(dds$Hugo_Symbol)]

## Count how many duplicated Hugo Symbols remain 
num_duplicates <- length(unique(duplicated_genes))

## Display the results 
cat("Number of duplicated Hugo Symbols after filtering:", num_duplicates)
## 0 duplicates remain 

##Make row names Hugo_Symbols 
rownames(dds) <- rna_rsem[keep,1]

dds <- DESeq(dds)
res <- results(dds, contrast=c("CNV", "Diploid", "Deletion"))
res
res<-as.data.frame(res)

## Volcano Plot 
p1<- EnhancedVolcano (as.data.frame (res),
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                col=c('grey', 'grey', 'grey', 'red'), 
                colAlpha = 1)
p1

write.csv(as.data.frame(res),
file="SULT1A1 tcga.csv")


