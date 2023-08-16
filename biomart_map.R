#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART B:: DEG ##############################################
################################################################################
packages <- c("DESeq2", "biomaRt")
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE,
                         repos = "http://cran.r-project.org/")
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
################################################################################
library(DESeq2)
library(biomaRt)
library(stringr)
################################################################################
colnames(seq) <- gsub("\\.bam", "", colnames(seq))
seIdx <- match(colnames(seq), sampleTable$Run)
colData(seq) <- cbind(Run.1 = colData(seq), sampleTable[seIdx, ])
colData(seq) <- colData(seq)[6:11]
################################################################################
# ddsFull <- DESeqDataSet(se, design = ~ Treatment)
ddsSeq <- DESeqDataSet(seq, design = as.formula(~ Condition))
################################################################################
# View(rownames(ddsFull))
# ddsFull$Treatment <- relevel(ddsFull$Treatment, "control")
################################################################################
# dds <- DESeq(ddsFull)
dds <- DESeq(ddsSeq, betaPrior=TRUE)
# resultsNames(dds)
res1 <- results(dds, contrast = c("Condition", "WT_pos.abeta", "WT_neg.abeta"))
res2 <- results(dds, contrast = c("Condition", "exon2del_pos.abeta", "exon2del_neg.abeta"))
res3 <- results(dds, contrast = c("Condition", "exon2del_neg.abeta", "WT_neg.abeta"))
res4 <- results(dds, contrast = c("Condition", "exon2del_pos.abeta", "WT_pos.abeta"))
################################################################################
saveRDS(dds, file = "./results/dExp.rds")
dds <- readRDS(file = "./results/dExp.rds")
################################################################################
res1$ensembl <- sapply( strsplit( rownames(res1), split="\\+" ), "[", 1 )
res2$ensembl <- sapply( strsplit( rownames(res2), split="\\+" ), "[", 1 )
res3$ensembl <- sapply( strsplit( rownames(res3), split="\\+" ), "[", 1 )
res4$ensembl <- sapply( strsplit( rownames(res4), split="\\+" ), "[", 1 )
################################################################################
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl)
idx1 <- match(res1$ensembl, genemap$ensembl_gene_id)
res1$entrez <- genemap$entrezgene[idx1]
res1$hgnc_symbol <- genemap$hgnc_symbol[idx1]

idx2 <- match(res2$ensembl, genemap$ensembl_gene_id)
res2$entrez <- genemap$entrezgene[idx2]
res2$hgnc_symbol <- genemap$hgnc_symbol[idx2]

idx3 <- match(res3$ensembl, genemap$ensembl_gene_id)
res3$entrez <- genemap$entrezgene[idx3]
res3$hgnc_symbol <- genemap$hgnc_symbol[idx3]

idx4 <- match(res4$ensembl, genemap$ensembl_gene_id)
res4$entrez <- genemap$entrezgene[idx4]
res4$hgnc_symbol <- genemap$hgnc_symbol[idx4]

################################################################################
res1$Length <- str_count(res1$hgnc_symbol)
res1$hgnc_symbol <- ifelse(is.na(res1$hgnc_symbol) | res1$Length == 0,
                           res1$ensembl, res1$hgnc_symbol)

write.csv(as.data.frame(res1[order(res1$padj, decreasing = FALSE, na.last = TRUE), ]),
          file = "F:/data/RNAseq/P190850/analysis/redo/WT_pos.abeta_vs_WT_neg.abeta.csv")
###

res2$Length <- str_count(res2$hgnc_symbol)
res2$hgnc_symbol <- ifelse(is.na(res2$hgnc_symbol) | res2$Length == 0,
                           res2$ensembl, res2$hgnc_symbol)

write.csv(as.data.frame(res2[order(res2$padj, decreasing = FALSE, na.last = TRUE), ]),
          file = "F:/data/RNAseq/P190850/analysis/redo/exon2del_pos.abeta_vs_exon2del_neg.abeta.csv")

###

res3$Length <- str_count(res3$hgnc_symbol)
res3$hgnc_symbol <- ifelse(is.na(res3$hgnc_symbol) | res3$Length == 0,
                           res3$ensembl, res3$hgnc_symbol)

write.csv(as.data.frame(res3[order(res3$padj, decreasing = FALSE, na.last = TRUE), ]),
          file = "F:/data/RNAseq/P190850/analysis/redo/exon2del_neg.abeta_vs_WT_neg.abeta.csv")

###

res4$Length <- str_count(res4$hgnc_symbol)
res4$hgnc_symbol <- ifelse(is.na(res4$hgnc_symbol) | res4$Length == 0,
                           res4$ensembl, res4$hgnc_symbol)

write.csv(as.data.frame(res4[order(res4$padj, decreasing = FALSE, na.last = TRUE), ]),
          file = "F:/data/RNAseq/P190850/analysis/redo/exon2del_pos.abeta_vs_WT_pos.abeta.csv")
################################################################################
