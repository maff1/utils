#!/usr/bin/Rscript

################################################################################
# GSEA KEGG PATHWAYS  ###########################################################

rm(list = ls())
gc()

# ENVIRONMENT ---------------------------------------------------------------- #
library("psych")
library("data.table")
library("dplyr")
require("AnnotationDbi")
require("org.Hs.eg.db")
library("fgsea")
library("pheatmap")
library("parallel")
library("factoextra")
library("cluster")
library("OmnipathR")
library("readr")

source("./CODE/f_utils.r")

# INPUT ---------------------------------------------------------------------- #
# ADNINIGHTINGALELONG_05_24_21_bl_residuals
# ADMCGUTMETABOLITESLONG_12_13_21_bl_residuals
# ADMCLIPIDOMICSMEIKLELABLONG_bl_residuals
adniMet = "ADMCLIPIDOMICSMEIKLELABLONG_bl_residuals"
# PARAMETERS ----------------------------------------------------------------- #
dirPathways = "./CONFIG/msigdb_v7.2_GMTs"
adniGene = "residuals_ADNI_Gene_Expression_Profile.rds"
# OUTPUT --------------------------------------------------------------------- #
outDir = "./CLEAN_DATA/METABOLITES/SPEARMAN"
# ---------------------------------------------------------------------------- #

rho <- readRDS(file.path(outDir, paste("corRho", adniMet, adniGene, sep = "_")))
pval <- readRDS(file.path(outDir, paste("corPval", adniMet, adniGene, sep = "_")))

resRho <- melt(rho, id.vars=colnames(rho), 
            measure.vars=rownames(rho),
            value.name = "rho")
resPval <- melt(pval, id.vars=colnames(pval), 
               measure.vars=rownames(pval),
               value.name = "fdr")
tbl <- data.frame(resRho, fdr = resPval[, 3])

tbl %>%
  group_by(Var2) %>%
  slice_min(n=10, order_by=fdr) %>%
write_delim(., file = file.path(outDir, "NETWORK", paste0("LIPID", "_rho_pairs.tab")),
            delim = "\t")

imTbl <- data.frame(name = unique(tbl$Var2),
                    feature = "metabolite")
imTbl_2 <- data.frame(name = unique(tbl$Var1),
                    feature = "gene")
out <- rbind(imTbl, imTbl_2)
write_delim(out, file = file.path(outDir, "NETWORK", paste0("GUT", "ANNOT.tab")),
            delim = "\t")

# coefs <- read.csv(file = file.path("./CLEAN_DATA/METABOLITES/MET_LMM", "COEF", 
#                           paste0("NG", "_LMM_COEF.csv")), row.names = NULL)

# FIX GENE ENTRIES WITH SYMBOLS CHANGED TO DATE (EXCEL)
  gene_symbol <- data.frame(
    symbol = as.character(unique(tbl$Var1))
  )
toMatch <- gene_symbol[grep("^X[0-9]{1,2}\\.", gene_symbol$symbol), ]
gene_symbol$symbol_fix <- mgsub::mgsub(
  string = gene_symbol$symbol, 
  pattern = c(toMatch, "ENSG00000122432", "ENSG00000266826", "ENSG00000250374",
              "ENSG00000211893", "ENSG00000211896", "ENSG00000211897"),
  replacement = c("SEPT2", "SEPT9", "SEPT11", "SEPT8", 
                  "MARCH6", "SEPT10", "MARCH7", "SEPT4", "MARCH9", "SEPT6",
                  "SEPT3", "SEPT1", "MARCH2", "MARCH3", "MARCH4", "MARCH1", "MARCH1",
                  "SEPT14", "SEPT7", "MARCH8", "DEC1", "MARCH11", "MARCH5",
                  "SPATA1", "IGBP1P2", "TRIM75", "IGHG2", "IGHG1", "IGHG3")
  )
match_tbl <- merge(tbl, gene_symbol, by.x="Var1", by.y="symbol")
  tbl <- data.frame(
    Var1 = match_tbl[["symbol_fix"]],
    Var2 = match_tbl[["Var2"]],
    rho = match_tbl[["rho"]],
    fdr = match_tbl[["fdr"]]
  )
  
# ANY BIAS TOWARDS PARTICULAR PROTEIN FAMILIES
geneOlfactory <- unique(grep("^OR", tbl$Var1, value = TRUE))
geneOlfactory <- geneOlfactory[geneOlfactory %nin% c("ORAI1",
                                                       "ORAI2",
                                                       "ORAI3",
                                                       "ORC1",
                                                       "ORC2",
                                                       "ORC3",
                                                       "ORC4",
                                                       "ORC5",
                                                       "ORC6")]
set.seed(2020)
geneOlfactory <- geneOlfactory[geneOlfactory %nin% sample(geneOlfactory, size=1, replace=FALSE)]
tbl <- tbl[tbl$Var1 %nin% geneOlfactory, ]
# ---------------------------------------------------------------------------- #

lsSets <- lapply(unique(tbl$Var2), 
                 function(met) {
                   df <- subset(tbl, Var2 == met)[c("Var1", "rho")]
                     vDF <- as.numeric(df$rho)
                       names(vDF) <- df$Var1
                     vDF <- sort(vDF, decreasing = FALSE)
                    return(vDF)
                    }); names(lsSets) <- unique(tbl$Var2)

# ---------------------------------------------------------------------------- #
# MAP SYMBOLS TO ENTREZID
#   getMatrixWithSelectedIds <- function(df, type, keys) {
#     x <- data.frame(genes = as.character(unique(df)))
#     geneSymbols <- mapIds(org.Hs.eg.db, 
#                           keys=x$genes, 
#                           column=type,
#                           keytype=keys, 
#                           multiVals="first")
#     
#  
#     # inds <- which(!is.na(geneSymbols))
#     # found_genes <- geneSymbols[inds]
#     # 
#     # 
#     # df2 <- df[names(found_genes), ]
#     # rownames(df2) <- found_genes
#     return(geneSymbols)
#   }
# toEntrez <- getMatrixWithSelectedIds(tbl$Var1, type="ENTREZID", keys="SYMBOL")

# KEEP ONLY NON-DISEASE ASSOCIATED PATHWAY TERMS
keggPath <- gmtPathways(file.path(dirPathways, "c2.cp.kegg.v7.2.symbols.gmt"))
# filtered_paths <- names(keggPath)[names(keggPath) %nin% grep(
# "disease|cancer|ASTHMA|LEUKEMIA|DIABETES|CARCINOMA|MELANOMA|INFECTION|SCLEROSIS|CARDIOMYOPATHY|MYOCARDITIS|LUPUS|AUTOIMMUNE|GLIOMA|AMYOTROPHIC|IMMUNODEFICIENCY|REJECTION|ALLOGRAFT|DEPRESSION", 
#        names(keggPath), ignore.case = TRUE, value = TRUE)]
# keggPath <- keggPath[filtered_paths]
keggMetabolic <- subset(read.delim(
  file = file.path(outDir, "kegg_metabolic.txt"), 
  sep = "\t"), metabolic == 1)$kegg
keggPath <- keggPath[keggMetabolic]

keggPath <- gmtPathways(file.path(dirPathways, "c2.cp.kegg.v7.2.symbols.gmt"))
filtered_paths <- names(keggPath)[names(keggPath) %nin% grep(
  "disease|cancer|ASTHMA|LEUKEMIA|DIABETES|CARCINOMA|MELANOMA|INFECTION|SCLEROSIS|CARDIOMYOPATHY|MYOCARDITIS|LUPUS|AUTOIMMUNE|GLIOMA|AMYOTROPHIC|IMMUNODEFICIENCY|REJECTION|ALLOGRAFT|DEPRESSION", 
  names(keggPath), ignore.case = TRUE, value = TRUE)]
keggPath <- keggPath[filtered_paths]

# ---------------------------------------------------------------------------- #
# RUN GSEA
lsResKEGG_2 <- safe_mclapply(lsSets, mc.cores = 50, stop.on.error=FALSE,
                             function(vGenes) {
                               kegg <- fgseaMultilevel(
                                 pathways = keggPath, 
                                 stats = vGenes, 
                                 minSize = 3, 
                                 maxSize = Inf,
                                 eps = 1e-50,
                                 scoreType = c("std", "pos", "neg"),
                                 nPermSimple = 1000
                                 )
                               return(kegg)
                             })
stopifnot(names(lsResKEGG) == names(lsSets))
saveRDS(lsResKEGG, file = file.path(outDir, "GSEA_KEGG_FILTER", paste0("GSEA_", adniMet, ".rds")))
f_extract_NES <- function(X) {
  indx <- lsResKEGG[[1]][["pathway"]]
    df <- X[match(indx, X[["pathway"]])][, c("pathway", "NES")]
  return(df)
}

lsNES <- lapply(lsResKEGG, f_extract_NES)
mtx <- scale(
  matrix(data = unlist(lapply(lsNES, function(x) x[["NES"]])), 
                nrow = nrow(lsNES[[1]]), 
                ncol = length(lsNES),
                byrow = FALSE,
                dimnames = list(
                  lsNES[[1]][[1]],
                  names(lsNES)
                )
  ), center = TRUE, scale = TRUE
)
# ---------------------------------------------------------------------------- #

# p-value clustering from ORA
# metrics <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# data_mtx <- t(mtx)
# f_ag_clustering <- function(mtype) {
#   df <- scale(data_mtx)
#   agnes(df, method = mtype)$ac
# }
# sapply(metrics, f_ag_clustering)
#   hc_ward <- agnes(mtx, method = "ward")
# pltree(hc_ward, cex = 0.6, hang = -1, main = "Dendrogram")
# fviz_nbclust(data_mtx, FUNcluster = hcut, 
#                  method = "gap_stat",
#                  k.max = 10, nboot = 100)

f_get_cluster_pval_fisher <- function(dat, K) {
  set.seed(789)
  d <- dist(t(dat), method = "euclidean")
  final_clust <- hclust(d, method = "ward.D2" )
  groups <- cutree(final_clust, k=K)
  final_data <- cbind(cluster = groups, t(dat))
return(final_data)
}
nes_clust <- f_get_cluster_pval_fisher(dat=mtx, K=4)
  f_to_excel <- function(xls) {
    DF <- as.data.frame(nes_clust)
openxlsx::write.xlsx(x = DF, 
                       file = file.path(outDir, "GSEA_KEGG_FILTER", paste0(adniMet, "_GSEA_NES.xlsx")), 
                     rowNames = TRUE)
  }
f_to_excel()

# f_get_cluster_pval_freq <- function(dat, K, lsORA, pvalThr) {
#   set.seed(789)
#   d <- dist(t(dat), method = "euclidean")
#   final_clust <- hclust(d, method = "ward.D2" )
#   groups <- cutree(final_clust, k=K)
#   final_data <- cbind(t(dat), cluster = groups)
#   
#   pathCut <- sapply(lsORA, function(x) subset(x, pval < pvalThr)[,1])
#   names(pathCut) <- final_data[, c("cluster")]
#   ls2mtx <- ComplexHeatmap::list_to_matrix(pathCut)
#   df <- data.frame(
#     cluster = dimnames(ls2mtx)[[2]],
#     t(ls2mtx),
#     row.names = NULL
#   )
#   setDT(df)
#   res <- df[, lapply(.SD, sum), by="cluster"]
#   print(
#     colnames(res)[apply(res, 1, which.max)]
#   )
# return(res)
# }
# pval_freq <- f_get_cluster_pval_freq(dat=mtx, lsORA=lsPval, K=4, pvalThr=0.05)

# pheatmap ------------------------------------------------------------------- #

set.seed(1231)
rowAnno <- data.frame(cluster = as.factor(nes_clust[, "cluster"]), 
                      row.names = rownames(nes_clust)
                      )
colorPalette = colorRampPalette(c("blue", "white", "red"))(100)  
pdf(file = file.path(outDir, "GSEA_KEGG_FILTER", paste0(adniMet, "_heatmap_GSEA_NES.pdf")), height = 20, width = 30)
pheatmap::pheatmap(mat = t(mtx),
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE,
                   clustering_method = "ward.D2",
                   cutree_rows = 4,
                   color = colorPalette,
                   annotation_row = rowAnno,
                   fontsize_row = 4,
                   border_color = "black",
                   main = paste("GSEA", adniMet, sep = "|")
                   )
dev.off()
# ---------------------------------------------------------------------------- #

# lsPaths <- readRDS(file = file.path(outDir, 
#             "GSEA", "GSEA_ADNINIGHTINGALELONG_05_24_21_bl_residuals.rds"))
# clusGroup <- openxlsx::read.xlsx(
#   xlsxFile = file.path(outDir, "GSEA",
#   "ADNINIGHTINGALELONG_05_24_21_bl_residuals_GSEA_NES.xlsx"), sheet = 1, rowNames = TRUE)
# clusGroup["metabolite"] <- rownames(clusGroup)
# spadj <- 0.05
# res_padj <- lapply(lsPaths, function(x) subset(x, padj < spadj))
#     stopifnot(rownames(clusGroup) == names(res_padj))
#     res <- data.frame(metabolite = rep(names(res_padj), sapply(res_padj, nrow)),
#                       do.call(rbind, res_padj)
#     )
#     res2 <- merge(res, clusGroup[, c("cluster", "metabolite")], by = "metabolite")
#       clusterName <- data.frame(table(res2$pathway, res2$cluster))
#         clusterName["freq_perc"] <- ifelse(clusterName$Var2 == 1, clusterName$Freq/392*100, 
#                                         ifelse(clusterName$Var2 == 2, clusterName$Freq/570*100,
#                                             ifelse(clusterName$Var2 == 3, clusterName$Freq/406*100,
#                                                ifelse(clusterName$Var2 == 4, clusterName$Freq/254*100, NA)
#                                         )))
# # 1 - (std(x) / (max(x) - min(x))
#         f_extract_NES <- function(X) {
#           indx <- lsPaths[[1]][["pathway"]]
#           df <- X[match(indx, X[["pathway"]])][, c("pathway", "NES")]
#           return(df)
#         }
#         
#         lsNES <- lapply(lsPaths, f_extract_NES)
#         mtx <- scale(
#           matrix(data = unlist(lapply(lsNES, function(x) x[["NES"]])), 
#                  nrow = nrow(lsNES[[1]]), 
#                  ncol = length(lsNES),
#                  byrow = FALSE,
#                  dimnames = list(
#                    lsNES[[1]][[1]],
#                    names(lsNES)
#                  )
#           ), center = TRUE, scale = TRUE
#         )
# ---------------------------------------------------------------------------- #
# all(dimnames(mtx)[[2]] == rownames(clusGroup))
# # metabolite = dimnames(mtx)[[2]]
# df <- data.frame(cluster = clusGroup$cluster,
#                  t(mtx)
#                  )
# setDT(df)
#   res3 <- df[, lapply(.SD, function(x) 1 - (sd(x) / (max(x) - min(x)))), by=cluster]
#     colnames(res3)[apply(res3, 1, which.min)]
    
# RASTER PLOT ---------------------------------------------------------------- #

# ng <- read.delim(file = file.path(outDir, "ng_sig_kegg.txt"), sep = "\t")
#   rownames(ng) <- gsub("KEGG_", "", ng$cluster)
#     ng$cluster <- NULL
#       colnames(ng) <- paste0("cluster", 1:4)
#         ng_log <- -log10(ng)
#         ng_txt <- data.frame(sapply(ng, 
#                                     function(x) ifelse(x > 0.05, "", format(x, digits = 3))),
#                              row.names = gsub("KEGG_", "", rownames(ng))
#         )
#       
#       bk1 <- seq(0,1.3,by=0.1)
#       bk2 <- c(seq(1.31,19,by=2),20)
#       bk <- c(bk1,bk2)
#       
#       my_palette <- c(colorRampPalette(colors = c("gray38", "gray38"))(n = length(bk1)-1),
#                       c(colorRampPalette(colors = c("yellow", "orange"))(n = length(bk2)-1)))
#       
#       set.seed(1231)
#       colorPalette = colorRampPalette(c("blue", "white", "red"))(100)  
#       pdf(file = file.path(outDir, paste0("RASTER_NIGHTINGALE", "_heatmap_ORA.pdf")), height = 20, width = 30)
#       pheatmap::pheatmap(mat = ng_log,
#                          cluster_rows = T, 
#                          cluster_cols = T,
#                          scale = "none",
#                          color = my_palette,
#                          breaks = bk,
#                          legend_breaks = c(0, 5, 10, 15, 20),
#                          legend_labels = c("0", "5", "10", "15", "-log10\n20"),
#                          fontsize_row = 12,
#                          fontsize = 12,
#                          border_color = "black",
#                          display_numbers = ng_txt,
#                          main = paste("ORA", "KEGG", "NIGHTINGALE", sep = " | ")
#       )
#       dev.off()
# ---------------------------------------------------------------------------- #

# ng <- readRDS(file = file.path(outDir, "GSEA_KEGG_FILTER",
#   "GSEA_ADNINIGHTINGALELONG_05_24_21_bl_residuals.rds"))
# 
# lspadj <- lapply(ng, function(x) subset(x, padj <0.05))
# goObo <- obo_parser(path = file.path("CONFIG", "go.obo"))
#   goObo_f <- as.data.frame(goObo$names[grep("bolic", goObo$names$name),])
# goa <- data.table::fread(file = file.path("CONFIG", "goa_human.gaf.gz"), skip = 13)
#       goa_sub_p <- subset(goa, V9 == "P" & V15 == "UniProt")[, c(5,11)]
#       goa_sub_p$symbol <- ifelse(nchar(goa_sub_p$V11) <= 1, NA, goa_sub_p$V11)
#       goa_sub_p <- goa_sub_p[!is.na(goa_sub_p$symbol), ]
#       goa_sub_p <- data.frame(goid = goa_sub_p$V5, symbol = goa_sub_p$symbol)
#       
#       
# goa_sub_p1 <- goa_sub_p[grep(pattern = "\\|", goa_sub_p$symbol),]
#   ls_goa_sub_p1 <- lapply(
#     unique(goa_sub_p1$goid), function(x) subset(goa_sub_p1, goid == x))
#   names(ls_goa_sub_p1) <- unique(goa_sub_p1$goid)
#   f_clean_go <- function(x) {
#     vSymbl <- unlist(strsplit(x[["symbol"]], split = "\\|"))
#       vSymbl[!grepl("\\/", vSymbl)]
#   }
# goa_sub_p1_clean <- lapply(ls_goa_sub_p1, f_clean_go)
#     names(goa_sub_p1_clean) <- gsub("\\:", "_", names(goa_sub_p1_clean))
#       indx <- gsub("\\:", "_", unique(goObo_f$term))
#         remove_null_lst <- function(x)  x[!sapply(x, is.null)]
#       goa_sub_p1_out <- remove_null_lst(goa_sub_p1_clean[indx])
#       goa_sub_p1_out <- goa_sub_p1_out[sapply(goa_sub_p1_out, function(x) length(x) > 2)]
# saveRDS(goa_sub_p1_out, file = file.path("CONFIG", "GO_BP_FILTER_METABOLIC.rds"))
# namesGO <- goObo_f[gsub("\\:", "_", goObo_f$term) %in% names(goa_sub_p1_out),]
# saveRDS(namesGO, file = file.path("CONFIG", "GO_BP_FILTER_METABOLIC_NAMES.rds"))
# ---------------------------------------------------------------------------- #

# NETWORK
dt <- rbindlist(lsResKEGG_2)
  dt <- dt[, metabolite := rep(names(lsResKEGG_2), each=nrow(lsResKEGG_2[[1]]))]
  data.frame(dt$pathway[1], dt$metabolite[1], gene = unlist(dt$leadingEdge[1]),
             dt$pval[1], dt$padj[1], dt$NES[1])
  
  
  res <- lapply(1:nrow(dt), function(x) {
    data.frame(dt$pathway[[x]], 
               dt$metabolite[[x]], 
               gene = unlist(dt$leadingEdge[[x]]),
               dt$pval[[x]], 
               dt$padj[[x]], 
               dt$NES[[x]])
    })
  DTres <- rbindlist(res)
    colnames(DTres) <- gsub("dt\\.|\\.\\.x\\.\\.", "", colnames(DTres))
  saveRDS(DTres, file = file.path(outDir, "KEGG", "gsea_ng_input_network.rds"))
  
# HC CLUSTERING -------------------------------------------------------------- # 
# metabolites using pathway (GSEA) NES
  # "gsea_ng_input_network.rds"
  f_get_cluster_nes <- function(dat, K) {
    set.seed(789)
    m <- reshape2::dcast(dat, formula = metabolite ~ pathway, value.var = "NES", fun.aggregate = mean)
      m <- as.data.frame(m[, -1], row.names = m$metabolite)
    d <- dist(m, method = "euclidean")
    final_clust <- hclust(d, method = "ward.D2")
    groups <- cutree(final_clust, k=K)
    final_data <- data.frame(cluster = groups, metabolite = rownames(m))
    final <- merge(dat, final_data, by = "metabolite", all.x = TRUE)
    return(final)
  }
  nes_clust <- f_get_cluster_nes(dat=DTres, K=4)
saveRDS(nes_clust, file = file.path(outDir, "KEGG", "gsea_ng_input_network_cluster.rds"))
write_delim(nes_clust, file = file.path(outDir, "KEGG", "gsea_ng_input_network_cluster.tab"),
            delim = "\t")

unique(subset(nes_clust, pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION")[,3])
##
met_path <- data.frame(id = paste(nes_clust$metabolite, nes_clust$pathway, sep="_"),
                       nes_clust[, -3]
)
met_path <- met_path[!duplicated(met_path$id),]
  met_path <- subset(met_path, padj <0.1)
    met_path <- met_path[, -1]
    met_path <- data.frame(met_path, 
                           padjLog = -log10(met_path$padj), 
                           NESabs = abs(met_path$NES),
                           NESdirection = ifelse(met_path$NES <0, "negative", "positive")
    )
write_delim(met_path, file = file.path(outDir, "KEGG", "gsea_ng_input_network_cluster_met_path.tab"),
                delim = "\t")

# READ AND FILTER
ng <- readRDS(file = file.path(outDir, "KEGG", "gsea_ng_input_network_cluster.rds"))
fil <- unique(read.delim(file = file.path(outDir, "KEGG", "filter_paths"), header = F)$V1)
  ng$id <- paste(ng$metabolite, ng$pathway, sep="_")
    ng <- ng[ng$id %in% fil, ]
    write_delim(ng, file = file.path(outDir, "KEGG", "add_genes_network.tab"),
                delim = "\t")

# NG BL RESIDUALS
bl <- readRDS(file = file.path("CLEAN_DATA", "METABOLITES", "RESIDUALS",
"ADNINIGHTINGALELONG_05_24_21_bl_residuals.rds")
)

gene <- readRDS(file = file.path("CLEAN_DATA", "gene", 
                                 "residuals_ADNI_Gene_Expression_Profile.rds"))
labelDx <- readRDS(file = file.path("CLEAN_DATA", "gene",
                                    "microarray_gene_merged_cova.rds"))
cova <- subset(readRDS(file = file.path("CLEAN_DATA", "covariates",
                                 "covariates_data_matrix.rds")),
               VISCODE == "bl")[, c(1, 18)]
  cova <- cova[cova$RID %in% labelDx$RID, ]
  
# LIMMA object
library(limma)
library(Biobase)
  
  g <- data.frame(t(gene), row.names = colnames(gene))
    colnames(g) <- paste0("rid", cova$RID)
      g <- as.matrix(g)
    colData <- data.frame(cova, row.names = paste0("rid", cova$RID))
  stopifnot(rownames(colData) == colnames(g))
  phenoData <- new("AnnotatedDataFrame", data=colData)
eset <- ExpressionSet(assayData=as.matrix(g), phenoData=phenoData)
subEset <- eset[, eset$DX %in% c("ALZ","CTR")]
subEset$DX <- droplevels(subEset$DX)

design <- model.matrix(~0 + DX, data=subEset)
  colnames(design) = gsub("DX","",colnames(design))
    colnames(design) = make.names(colnames(design))
contrast <- makeContrasts(contrasts="ALZ-CTR", levels=design)
fit1 <- lmFit(exprs(subEset), design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)

res = topTable(fit3, coef=c("ALZ-CTR"), n=Inf)
write.csv(res, file = file.path(outDir, "KEGG", "deg_ALZ-CTR.csv"), 
          row.names = T)
fc <- data.frame(gene = rownames(res), fc = res$logFC)
write_delim(fc, file = file.path(outDir, "KEGG", "fc_genes_network.tab"),
            delim = "\t")

# METABOLITE
met <- readRDS(file = file.path("CLEAN_DATA", "METABOLITES", "MET_LMM", "NG",
               "NG_lxm_model_time_interaction_met_DX.rds")
)

beta <- lapply(met, function(x) subset(x, coef.term == "DXALZ")[,c(2,8)])
  beta_df <- do.call(rbind, beta)
    beta_df <- data.frame(metabolite = rownames(beta_df), beta = beta_df[,1],
                          pval = -log10(beta_df[,2]))
write_delim(beta_df, file = file.path(outDir, "KEGG", "fc_metabolite_network_v2.tab"),
                delim = "\t")
  
# GUT ------------------------------------------------------------------------ #
gut <- readRDS(file = file.path(outDir, "GSEA_v2", 
                                "GSEA_ADMCGUTMETABOLITESLONG_12_13_21_bl_residuals.rds"))

# NETWORK
f_to_network <- function(lsResKEGG_2, dataset) {
dt <- rbindlist(lsResKEGG_2)
dt <- dt[, metabolite := rep(names(lsResKEGG_2), each=nrow(lsResKEGG_2[[1]]))]
data.frame(dt$pathway[1], dt$metabolite[1], gene = unlist(dt$leadingEdge[1]),
           dt$pval[1], dt$padj[1], dt$NES[1])


res <- lapply(1:nrow(dt), function(x) {
  data.frame(dt$pathway[[x]], 
             dt$metabolite[[x]], 
             gene = unlist(dt$leadingEdge[[x]]),
             dt$pval[[x]], 
             dt$padj[[x]], 
             dt$NES[[x]])
})
DTres <- rbindlist(res)
colnames(DTres) <- gsub("dt\\.|\\.\\.x\\.\\.", "", colnames(DTres))
saveRDS(DTres, file = file.path(outDir, "KEGG", paste0("gsea_", dataset, "_input_network.rds")))
return(DTres)
}

DTres <- f_to_network(gut, "gut")

# HC CLUSTERING -------------------------------------------------------------- # 
# metabolites using pathway (GSEA) NES
# "gsea_gut_input_network.rds"
f_get_cluster_nes <- function(dat, K) {
  set.seed(789)
  m <- reshape2::dcast(dat, formula = metabolite ~ pathway, value.var = "NES", fun.aggregate = mean)
  m <- as.data.frame(m[, -1], row.names = m$metabolite)
  d <- dist(m, method = "euclidean")
  final_clust <- hclust(d, method = "ward.D2")
  groups <- cutree(final_clust, k=K)
  final_data <- data.frame(cluster = groups, metabolite = rownames(m))
  final <- merge(dat, final_data, by = "metabolite", all.x = TRUE)
  return(final)
}
nes_clust <- f_get_cluster_nes(dat=DTres, K=4)
saveRDS(nes_clust, file = file.path(outDir, "KEGG", "gsea_gut_input_network_cluster.rds"))
write_delim(nes_clust, file = file.path(outDir, "KEGG", "gsea_gut_input_network_cluster.tab"),
            delim = "\t")

##
met_path <- data.frame(id = paste(nes_clust$metabolite, nes_clust$pathway, sep="_"),
                       nes_clust[, -3]
)
met_path <- met_path[!duplicated(met_path$id),]
met_path <- subset(met_path, padj <0.1)
met_path <- met_path[, -1]
met_path <- data.frame(met_path, 
                       padjLog = -log10(met_path$padj), 
                       NESabs = abs(met_path$NES),
                       NESdirection = ifelse(met_path$NES <0, "negative", "positive")
)
write_delim(met_path, file = file.path(outDir, "KEGG", "gsea_gut_input_network_cluster_met_path.tab"),
            delim = "\t")


# METABOLITE EXPRESSION
f_met_expr <- function(dataset) {
  met <- readRDS(file = file.path("CLEAN_DATA", "METABOLITES", "MET_LMM", dataset,
                                paste0(dataset, "_lxm_model_time_interaction_met_DX.rds"))
)
beta <- lapply(met, function(x) subset(x, coef.term == "DXALZ")[,c(2,8)])
beta_df <- do.call(rbind, beta)
beta_df <- data.frame(metabolite = rownames(beta_df), beta = beta_df[,1],
                      pval = -log10(beta_df[,2]))
write_delim(beta_df, file = file.path(outDir, "KEGG", paste0(dataset, "_fc_metabolite_network_v2.tab")),
            delim = "\t")
}
f_met_expr("GUT")

# ADD genes to network

f_add_gene_to_network <- function(cluster, dataset) {
  ng <- cluster
fil <- unique(read.delim(file = file.path(outDir, "KEGG", paste0(dataset,"_filter_paths")), header = F)$V1)
ng$id <- paste(ng$metabolite, ng$pathway, sep="_")
ng <- ng[ng$id %in% fil, ]
write_delim(ng, file = file.path(outDir, "KEGG", paste0(dataset, "_add_genes_network.tab")),
            delim = "\t")
}

f_add_gene_to_network(cluster=nes_clust, dataset="gut")
