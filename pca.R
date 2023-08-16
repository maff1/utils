################################################################################
# DIM REDUCTION PCA ## MAFF ## V001 ############################################

rm(list = ls())

dir <- "/mnt/marcodata/PROJECTS/ADNI/LONGITUDINAL/LONI-ADNI"
setwd(dir)

library(dplyr)
library(factoextra)
library(cowplot)
source("./CODE/f_utils.r")

# INPUT ---------------------------------------------------------------------- #
lsFiles <- list.files(path = "./REVISION/METABOLITES/MERGED", 
                      pattern = ".rds", 
                      full.names = T, ignore.case = T, recursive = T)
# exclusions
cancerRID <- "./REVISION/SENSITIVITY/cancer_exclusions_rid.tab"
sensiFile <- "./REVISION/SENSITIVITY/medication_disease_binary.rds"
# PARAMETERS ----------------------------------------------------------------- #
# OUTPUT --------------------------------------------------------------------- #
outDIR = "./REVISION/METABOLITES/PCA/MERGED"
# ---------------------------------------------------------------------------- #

f_pca_longitudinal <- function(n=n, excludeFile=sensiFile) {
  adniMet <- lsFiles[n]
  #sRid <- read.delim(file = exclude, col.names = "rid", header = F)
  sSensi <- readRDS(file = excludeFile)
    df <- readRDS(adniMet)
  #df <- merge(df, sSensi, by="RID_VISIT", all.x=TRUE)
}
  if(excludeCancer == "no") {
    df <- readRDS(adniMet)
      df <- merge(df, sSensi, by="RID_VISIT", all.x=TRUE)
  } else {
    df <- subset(readRDS(adniMet), cancer != 1)
  }
    
    stopifnot(is.data.frame(df))
    sMet <- colnames(df)[seq(grep("^VISIT$", colnames(df))+1, 
                              grep("^VISCODE2$", colnames(df))-1,by=1)]
    sCova <- colnames(df)[colnames(df) %nin% sMet]
    resPca <- cbind.data.frame(df[sCova],
                               df[sMet],
                               prcomp(df[sMet], 
                                      center = TRUE, 
                                      scale. = TRUE)$x[,1:40])
return(resPca)
}
res <- f_pca_longitudinal(n=3)
# ---------------------------------------------------------------------------- #



# permuted PCA bl ------------------------------------------------------------ #

f_wrap_pca_perm <- function(dat) {
lsVisit <- lapply(unique(dat$VISIT), function(x) dat[dat$VISIT %in% x,])
names(lsVisit) <- unique(dat$VISIT)
sMet <- colnames(dat)[seq(grep("^EXAMDATE.x$", colnames(dat))+1, 
                          grep("^EXAMDATE.y$", colnames(dat))-1,by=1)]

f_timer_permPCA <- function (x) {
  start <- Sys.time()
  permPca <- pca_eigenperm(x)
  end <- Sys.time()
total_time <- as.numeric(end - start, units = "mins")
  print(total_time)
return(permPca)
}

bl_pca <- prcomp(lsVisit$bl[sMet], center = TRUE, scale. = TRUE)
# permuation pca 1000 times
permPca <- f_timer_permPCA(lsVisit$bl[sMet])
out <- list(defaPCA=bl_pca, permPCA=permPca)
return(out)
}

lsDat <- list(ng=datNG, gut=datGUT, lipid=datLIPID)
resPCA <- lapply(lsDat, f_wrap_pca_perm)
saveRDS(resPCA, file = file.path(outDIR, "PCA", "pca_bl_metDatasets.rds"))

f_barPlotPermPCA<-function(permPCA, defaPCA, npc=9, ...){
  fa_pca_rand95<- 
    data.frame(Random_Eigenvalues = sapply(permPCA, quantile, 0.95)) %>%
    mutate(PC = colnames(defaPCA$rotation)) %>%
    cbind(Eigenvalues = defaPCA$sdev^2)
  fa_pca_rand95_long<-
    gather(fa_pca_rand95[1:npc, ], key = Variable, value = Value, -PC)
  fa_pca_rand95_long$PC <- factor(fa_pca_rand95_long$PC,
                                  levels = paste0("PC", seq(1, npc, by=1)))
  barPlot<- fa_pca_rand95_long %>%
    mutate(PC = as.numeric(gsub(pattern = "PC", "", PC))) %>%
    ggplot(aes(PC, Value, fill = Variable)) +
    scale_x_continuous(breaks=seq(1, npc, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    geom_bar(stat = "identity", position = position_dodge())+
    labs(y=bquote(bold(Eigenvalue~(sdev^2))), x="", fill= "") +
    theme_classic() + 
    theme(axis.title.y=element_text(face="bold"),
          axis.title.x=element_text(face="bold"),
          legend.title=element_text(face="bold"),
          axis.text.x=element_text(angle = 45, size=6),
          legend.position = c(0.6, .8)) +
    labs(x = "PCs")
  return(barPlot)
}
# ---------------------------------------------------------------------------- #

ls_plot_pca_perm <- lapply(names(lsDat), 
                           function(x) {
                             if(x == "ng") {
                               out <- f_barPlotPermPCA(permPCA=resPCA[[x]]$permPCA, 
                                                       defaPCA=resPCA[[x]]$defaPCA, 
                                                       npc=12)
                             } else if(x == "gut") {
                               out <- f_barPlotPermPCA(permPCA=resPCA[[x]]$permPCA, 
                                                       defaPCA=resPCA[[x]]$defaPCA, 
                                                       npc=15)
                             } else {
                               out <- f_barPlotPermPCA(permPCA=resPCA[[x]]$permPCA, 
                                                       defaPCA=resPCA[[x]]$defaPCA, 
                                                       npc=35)
                               }
                             })

plot_grid(plotlist = ls_plot_pca_perm, 
          labels = "auto", 
          nrow = 1, 
          ncol = 3,
          rel_widths = c(1,1,2))

f_cumVarPlot <- function(bl_pca, npc) {
plot_pca <- get_eigenvalue(bl_pca) %>%
  mutate(PCs = rownames(.)) %>%
  mutate(PCs = as.numeric(gsub(pattern = "Dim.", "", PCs)), 
         cumVar = paste0(round(cumulative.variance.percent, digits = 1), "%")) %>%
  filter(PCs <= npc) %>%
  ggplot(., aes(x = PCs, y = variance.percent, label = cumVar)) +
  scale_x_continuous(breaks=seq(1, npc, 1)) +
  geom_point() +
  geom_line() +
  geom_text(nudge_x = 0.40, nudge_y = 0.75, check_overlap = T, size = 3) +
  theme_classic() +
  theme(axis.title = element_text(face ="bold"),
        axis.text.x=element_text(angle = 45, size=6),
        plot.title = element_blank()) +
  labs(x = "PCs", y = "explained variance (%)")
return(plot_pca)
}

ls_cumVarPlot <- lapply(names(lsDat), 
                        function(x) {
                          if(x == "ng") {
                            out <- f_cumVarPlot(resPCA[[x]]$defaPCA, 
                                                    npc=12)
                          } else if(x == "gut") {
                            out <- f_cumVarPlot(resPCA[[x]]$defaPCA, 
                                                    npc=15)
                          } else {
                            out <- f_cumVarPlot(resPCA[[x]]$defaPCA, 
                                                    npc=35)
                          }
                        })

p1 = plot_grid(plotlist = ls_plot_pca_perm, 
          labels = "auto", 
          nrow = 1, 
          ncol = 3,
          rel_widths = c(1,1,1.25))
p2 = plot_grid(plotlist = ls_cumVarPlot, 
          labels = c("d", "e", "f"), 
          nrow = 1, 
          ncol = 3,
          rel_widths = c(1,1,1.25))


pPCA = plot_grid(p1, p2, labels = NULL, ncol = 1)
ggsave(filename = file.path(outDIR, "PCA", "plot_pca_bl_metDatasets.png"), 
       plot = pPCA, 
       dpi = 300,
       height = 7,
       width = 20)
saveRDS(pPCA, file = file.path(outDIR, "PCA", "plot_pca_bl_metDatasets.rds"))

# PCA PER VISIT (BL, M12, M24) ----------------------------------------------- #

f_wrap_pca <- function(dat) {
  lsVisit <- lapply(unique(dat$VISIT), function(x) dat[dat$VISIT %in% x,])
  names(lsVisit) <- unique(dat$VISIT)
  sMet <- colnames(dat)[seq(grep("^EXAMDATE.x$", colnames(dat))+1, 
                            grep("^EXAMDATE.y$", colnames(dat))-1,by=1)]
  vCova <- colnames(dat)[colnames(dat) %nin% sMet]
  v_pca <- lapply(names(lsVisit), 
                   function(v) {
                     cbind.data.frame(lsVisit[[v]][vCova],
                                      lsVisit[[v]][sMet],
                                      prcomp(lsVisit[[v]][sMet], 
                                      center = TRUE, 
                                      scale. = TRUE)$x[,1:40]
                     )
                   }
  ); names(v_pca) <- names(lsVisit)
return(v_pca)
}

lsDat <- list(ng=datNG, gut=datGUT, lipid=datLIPID)
resPCAvisit <- lapply(lsDat, f_wrap_pca)

f_merge_pca_covar <- function(x) {
mPCA <- do.call(rbind.data.frame, 
                lapply(c("bl", "m12", "m24"), function(v) resPCAvisit[[x]][[v]])
                )
  rownames(mPCA) <- NULL
return(mPCA)
}
mPCA <- lapply(names(lsDat), f_merge_pca_covar)
names(mPCA) <- names(lsDat)
saveRDS(mPCA, file = file.path(outDIR, "PCA", "pca_merged_metDatasets.rds"))
#summary(resPCAvisit$ng$bl)$importance[2,]*100
# ---------------------------------------------------------------------------- #
