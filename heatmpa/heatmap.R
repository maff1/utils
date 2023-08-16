################################################################################
# Save top20 from fig4b as obj

library(data.table)
library(DESeq2)
library(ComplexHeatmap)

load("top20.RData")

# filled circles and triangles (16,17)
nshapes = c(16, 17) 
nshapes <- nshapes[as.numeric(DF$abeta)]

colData(vsp) <- DataFrame(colData(vsp), DF[,19:21])
vsp$ipsc <- gsub(".CLU.exon2.hom.KO|.M3.36S", "", vsp$ipsc)
# ---------------------------------------------------------------------------- #
f_wrap_heatmap<-function(dat,ntop,show_row_names){
  
  topVarGenes <- head(order(-rowVars(assay(dat))),ntop)
  
  ha_top = HeatmapAnnotation(
    abeta = anno_simple(dat$ipsc, 
                        col = c("CTR" = "blue", "A4" = "red", "D1" = "orange"),
                        pch = nshapes),
    gp = gpar(col = "black"),
    show_legend = FALSE,
    annotation_label = ""
  )
  
  mat = scale(assay(dat)[topVarGenes,])
  
  set.seed(123)
  Heatmap(mat, 
          name="z-score \n counts", 
          border = T, 
          show_column_names = FALSE,
          column_names_gp = gpar(fontsize = 9), 
          top_annotation = ha_top,
          show_heatmap_legend = TRUE,
          row_labels=rowData(dat)[topVarGenes,]$SYMBOL,
          show_row_names = show_row_names
          
  )
}
# ---------------------------------------------------------------------------- #

f_wrap_heatmap2<-function(dat,ntop,show_row_names){
  
  topVarGenes <- head(order(-rowVars(assay(dat))),ntop)
  
  ha_top = HeatmapAnnotation(abeta = anno_simple(
    dat$IPSC, 
    col = c("CTR" = "blue", "A4" = "red", "D1" = "orange"),
    pch = nshapes),
    gp = gpar(col = "black"),
    show_legend = FALSE,
    annotation_label = ""
  )
  
  mat = scale(assay(dat)[topVarGenes,])
  
  set.seed(123)
  Heatmap(mat, 
          name="z-score \n counts", 
          border = T, 
          show_column_names = FALSE,
          column_names_gp = gpar(fontsize = 9), 
          top_annotation = ha_top,
          left_annotation = ha_top,
          show_heatmap_legend = TRUE,
          row_labels=rowData(dat)[topVarGenes,]$SYMBOL,
          show_row_names = show_row_names
          
  )
}
# ---------------------------------------------------------------------------- #

heat_20 = grid::grid.grabExpr(draw(f_wrap_heatmap(dat=vsp,ntop=20,show_row_names=TRUE)))
# ---------------------------------------------------------------------------- #
library(circlize)

## legend abeta

lg1 <-Legend(title = "",  pch = c(16, 17), size = unit(6, "mm"),
         labels = gt_render(c("A&beta;(-)" , "A&beta;(+)")),
         at = c(0, 1),
         type = "points",
         direction = "vertical",
         labels_gp = gpar(fontsize = 14, fontface = "bold"),
         break_dist = 2,
         grid_height = unit(10, "mm"),
         grid_width = unit(10, "mm")
  )
## A4, CTR, D1
lg2 <-
  Legend(title = "",  
         #pch = c(15,15,15), 
         size = unit(6, "mm"),
         labels = c("A4", "CTR", "D1"),
         #type = "points",
         direction = "vertical",
         legend_gp = gpar(fill = c("red", "blue", "orange")),
         labels_gp = gpar(fontsize = 14, fontface = "bold"),
         #break_dist = 2,
         #grid_height = unit(10, "mm"),
         #grid_width = unit(10, "mm")
  )
# ---------------------------------------------------------------------------- #

heatMap20 = grid::grid.grabExpr(
  draw(f_wrap_heatmap(dat=vsp,ntop=20,show_row_names=TRUE), 
     heatmap_legend_list = list(lg1, lg2))
)
