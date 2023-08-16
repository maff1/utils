################################################################################
# summary enrichment plots (revision SreportsEvie)

library(ggplot2)
library(dplyr)
library(cowplot)

# ---------------------------------------------------------------------------- #
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# ---------------------------------------------------------------------------- #

# Figure A
go <- openxlsx::read.xlsx("export-GO-tabl.xlsx", sheet = 1)
  go_bp = subset(go, GO == "BP")
gr_bp = c(158,107,126,54,74,86,83) # DE with GO terms for BP
go_bp$ratio = ifelse(go_bp$Cluster == "res1", go_bp$Count/158,
         ifelse(go_bp$Cluster == "res2", go_bp$Count/107,
           ifelse(go_bp$Cluster == "res3", go_bp$Count/126,
             ifelse(go_bp$Cluster == "res4", go_bp$Count/54,
               ifelse(go_bp$Cluster == "res5", go_bp$Count/74,
                 ifelse(go_bp$Cluster == "res6", go_bp$Count/86,
                   ifelse(go_bp$Cluster == "res7", go_bp$Count/83, NA)
                 ))))))

go_sub <- rev(c("myelination", "cell−substrate adhesion", 
  "endochondral bone morphogenesis", "growth plate cartilage morphogenesis", 
  "chondrocyte morphogenesis", "growth plate cartilage chondrocyte morphogenesis", 
  "chondrocyte morphogenesis involved in endochondral bone morphogenesis", 
  "cell−cell adhesion via plasma−membrane adhesion molecules", 
  "homophilic cell adhesion via plasma membrane adhesion molecules", "learning", 
  "extracellular structure organization", "extracellular matrix organization", 
  "negative regulation of protein serine/threonine kinase activity", 
  "ERK1 and ERK2 cascade", "negative regulation of MAPK cascade", 
  "response to corticosterone", "regulation of ERK1 and ERK2 cascade", 
  "heart morphogenesis", "epithelial cell proliferation", 
  "negative regulation of ERK1 and ERK2 cascade", "regulation of MAP kinase activity", 
  "epithelial tube morphogenesis","negative regulation of MAP kinase activity")
)
go_bp2 <- go_bp %>%
  filter(Description %in% go_sub) %>%
  mutate(Description = firstup(Description))
  
plot_go = ggplot(go_bp2, aes(x=Cluster, y=Description)) +
  geom_point(aes(size = ratio, color = p.adjust)) +
  scale_color_gradientn(
    trans = "log10",
    colours = c("red", "yellow", "blue"),
    values = c(0, scales::rescale(log10(0.01), from = log10(range(go_bp2$p.adjust))), 1)
  ) +
    scale_y_discrete(limits = firstup(rev(go_sub)),
                     labels=c("Extracellular structure organization"=expression(bold(Extracellular~structure~organization)), 
                              "Extracellular matrix organization"=expression(bold(Extracellular~matrix~organization)),parse=TRUE)) +
    scale_x_discrete(labels=paste0("res", 1:7, "\n", "(", gr_bp, ")")) +
    xlab("") + ylab("") + labs(color = "adj.\npvalue", size = "gene ratio", title = "GO::Biological Process (BP)") +
    theme(plot.title = element_text(hjust = 0.5))
# ---------------------------------------------------------------------------- #

# GSEA WITH KEGG
gs_kegg <- openxlsx::read.xlsx("GSEA_KEGG.xlsx", sheet = 1)
  numCols <- c("pval","padj","NES","size")
  gs_kegg[numCols] <- lapply(gs_kegg[numCols], as.numeric)
  gs_kegg$status <- ifelse(gs_kegg$NES <0 & gs_kegg$padj <0.1, "supressed", 
                             ifelse(gs_kegg$NES >0 & gs_kegg$padj <0.1, "activated", "ns"))
  gs_kegg$status <- factor(gs_kegg$status, levels = c("activated","supressed","ns"), labels = c("activated","supressed","ns"))
  gs_kegg$pathway <- gsub("\\."," ", tolower(gs_kegg$pathway))
  
  
pathsK  <- gs_kegg %>% filter(padj < 0.1) %>% select(pathway)
plot_gsea = gs_kegg %>%
  filter(pathway %in% unique(pathsK$pathway)) %>%
  mutate(pathway = firstup(pathway)) %>%
  mutate(pathway = gsub("Ecm receptor interaction", "ECM-receptor interaction", pathway)) %>%
  ggplot(., aes(x=NES, y=pathway, color=status, xmin=0, xmax=NES)) +
  scale_color_manual(values = c("red","blue","grey")) + 
    geom_pointrange() +
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.75) +
    scale_y_discrete(
                     labels=c("ECM-receptor interaction"=expression(bold(ECM-receptor~interaction)), 
                              "Focal adhesion"=expression(bold(Focal~adhesion)),parse=TRUE)) +
    facet_grid(~group) +
    theme_bw() +
    ylab("") + labs(title = "GSEA::KEGG") +
    theme(plot.title = element_text(hjust = 0.5))
# ---------------------------------------------------------------------------- #

# SPIA-KEGG
# Topological - gene set enrichment analysis across groups comparisons. pSize is
# the number of genes on the pathway; NDE is the number of DE genes per pathway; 
# tA is the observed total perturbation accumulation in the pathway; pNDE is the 
# probability to observe at least NDE genes on the pathway using a hypergeometric model; 
# pPERT is the probability to observe a total accumulation more extreme than tA only by chance; 
# pG is the p-value obtained by combining pNDE and pPERT; pGFdr and pGFWER are the False Discovery Rate 
# and respectively Bonferroni adjusted global p-values; and the Status gives the direction in which the pathway 
# is perturbed (activated or inhibited).

spia <- openxlsx::read.xlsx("SPIA_KEGG.xlsx", sheet = 1)
  spia$status <- ifelse(spia$tA <0 & spia$pGFdr <0.1, "supressed", 
                        ifelse(spia$tA >0 & spia$pGFdr <0.1, "activated", "ns"))
  spia$status <- factor(spia$status, levels = c("activated","supressed","ns"), labels = c("activated","supressed","ns"))
  
  pathsKspia  <- spia %>% filter(pGFdr < 0.1) %>% select(Name)
  intP <- c(unique(spia$Name)[unlist(
    sapply(unique(pathsK$pathway), function(x) grep(x, unique(spia$Name), ignore.case = T))
  )], "ECM-receptor interaction")
  
  plot_spia = spia %>%
    filter(Name %in% intP) %>%
    ggplot(., aes(x=tA, y=Name, color=status, xmin=0, xmax=tA)) +
    scale_color_manual(values = c("red","blue","grey")) + 
    geom_pointrange() +
    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.75) +
    scale_y_discrete(
      labels=c("ECM-receptor interaction"=expression(bold(ECM-receptor~interaction)), 
               "Focal adhesion"=expression(bold(Focal~adhesion)),parse=TRUE)) +
    facet_grid(~Cluster) +
    theme_bw() +
    ylab("") + labs(title = "SPIA::KEGG") +
    theme(plot.title = element_text(hjust = 0.5))
# ---------------------------------------------------------------------------- #

pa <- plot_grid(plot_go, heatMap20, labels = c("a","b"), ncol=1, nrow=2, rel_heights = c(1,1))
pb <- plot_grid(plot_gsea, plot_spia, labels = c("c","d"), ncol=1, nrow=2, rel_heights = c(1,1))
pdf(file = "figure4_rev.pdf", height = 12, width = 20)
plot_grid(
    pa,
    pb, ncol=2, nrow=1
  )
dev.off()
# ---------------------------------------------------------------------------- #

p4 = plot_grid(
  pa,
  pb, ncol=2, nrow=1
)

save_plot(filename = "figure4_rev.svg", p4, base_height = NULL, base_width = 20)
save_plot(filename = "figure4_rev.png", p4, base_height = NULL, base_width = 20)
save_plot(filename = "figure4_rev.pdf", p4, base_height = NULL, base_width = 20)
# ---------------------------------------------------------------------------- #

save.image(file = "figure4_rev.RData")

# Res1 	Aβ25-35 treated CTR vs. untreated CTR neurons 
# Res2 	Aβ25-35 treated A4 vs. untreated A4 neurons 
# Res3 	Aβ25-35 treated D1 vs. untreated D1 neurons 
# Res4 	Untreated CTR vs. untreated A4 neurons 
# Res5 	Untreated CTR vs. untreated D1 neurons 
# Res6 	Aβ25-35 treated CTR vs. treated A4 neurons 
# Res7 	Aβ25-35 treated CTR vs. treated D1 neurons 

##

lnames = data.frame( snames = paste(
  paste0("res",1:7),
  c("AB25-35 treated CTR vs. untreated CTR neurons",
    "AB25-35 treated A4 vs. untreated A4 neurons",
    "AB25-35 treated D1 vs. untreated D1 neurons",
    "Untreated CTR vs. untreated A4 neurons",
    "Untreated CTR vs. untreated D1 neurons",
    "AB25-35 treated CTR vs. treated A4 neurons",
    "AB25-35 treated CTR vs. treated D1 neurons"),
  sep = ": "),
  x = 1:7,
  y = 1:7
)
lnames$snames <- as.factor(lnames$snames)


plotLeg = ggplot(lnames, aes(x=x, y=y, color=snames)) +
  geom_point(color="white") +
  theme(legend.key = element_blank(), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_blank()) +
  scale_color_discrete(labels = c(expression(A[beta[25-35]]~treated~CTR~bold(.vs.)~untreated~CTR~neurons), 
                                  expression(A[beta[25-35]]~treated~A4~bold(.vs.)~untreated~A4~neurons),
                                  expression(A[beta[25-35]]~treated~D1~bold(.vs.)~untreated~D1~neurons),
                                  expression(Untreated~CTR~bold(.vs.)~untreated~A4~neurons),
                                  expression(Untreated~CTR~bold(.vs.)~untreated~D1~neurons),
                                  expression(A[beta[25-35]]~treated~CTR~bold(.vs.)~treated~A4~neurons),
                                  expression(A[beta[25-35]]~treated~CTR~bold(.vs.)~treated~D1~neurons)
                                  )
                       )


# Using the cowplot package
legend <- cowplot::get_legend(plotLeg)

grid.newpage()
grid.draw(legend)
