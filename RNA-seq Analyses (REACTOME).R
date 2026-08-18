#Paper written by Xiao et al. can be found with this link: https://www.nature.com/articles/s42255-018-0008-5

#Set working directory (WD) to the same folder containing this script (In the options above, click session, set WD, To Source File Location)
#Make sure your folder contains the R script, appropriate dataset, and downloaded background genesets (if required)

OUTDIR_3 <- "R_outputs_Reactome"

dir.create(OUTDIR_3, showWarnings = FALSE)
FIGDIR_3 <- file.path(OUTDIR_3, "figures_Reactome"); dir.create(FIGDIR_3, showWarnings = FALSE)



#GSEA 
#Running GSEA using Reactome
#Using clusterProfiler

library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)
library(enrichplot)
library(ReactomePA)


#Reading in dataset
RNA_data <- read.csv("RNA-seq Prototype Data.csv", header = TRUE)

#Setting healthy state as baseline and fibrosis as experimental group 
FIBROSIS_IS_NEGATIVE_LOG2FC <- TRUE

flip <- if (FIBROSIS_IS_NEGATIVE_LOG2FC) -1 else 1
RNA_data$stat_fib   <- flip * RNA_data$stat
RNA_data$log2FC_fib <- flip * RNA_data$log2FoldChange
message(sprintf("    %d genes total; %d with a valid adjusted p-value.",
                nrow(RNA_data), sum(!is.na(RNA_data$padj))))

RNA_data  <- RNA_data[!is.na(RNA_data$padj), ]
RNA_data <- RNA_data[!duplicated(RNA_data$gene_id), ]


#Setting pre-ranked gene list
#Always check if there are duplicated genes prior
react_vector <- RNA_data$stat_fib
names(react_vector) <- RNA_data$gene_id

final_react_vector <- na.omit(react_vector)

final_react_vector = sort(final_react_vector, decreasing = TRUE)

head(final_react_vector)

#Converting gene IDs for gsePathway function 
#Some genes will be lost (i.e., not all IDs convert)
rID <- bitr(names(react_vector), fromType = "SYMBOL",
            toType = "ENTREZID", OrgDb = org.Mm.eg.db)

#Removing duplicate IDs
dedup_rID = rID[!duplicated(rID[c("SYMBOL")]),]

#New dataframe with successfully mapped genes
rdf2 = RNA_data[RNA_data$gene_id %in% dedup_rID$SYMBOL,]

#Creating new column in rdf2 with corresponding ENTREZ IDs
rdf2$Z = dedup_rID$ENTREZID

#Creating gene vector from rdf2
react_gene_list <- rdf2$stat_fib
names(react_gene_list) <- rdf2$Z

react_gene_list <- na.omit(react_gene_list)

react_gene_list = sort(react_gene_list, decreasing = TRUE)

head(react_gene_list)


#Mark for reproducible results 
set.seed(2026)


#Running GSEA
#Setting background genes
react1 <- gsePathway(geneList = react_gene_list, 
                     organism = "mouse",
                     nPerm = 10000,
                     minGSSize = 15,
                     maxGSSize = 500,
                     pvalueCutoff = 0.02)

head(react1)

#Filtering out any row with NAs
react1@result <- react1@result[complete.cases(react1@result), ]

#Converting react1 to use Gene Symbols
react1_readable <- setReadable(react1, 
                            OrgDb = org.Mm.eg.db, 
                            keyType = "ENTREZID")

#Converting GSEA results to a data frame 
react1_df <- data.frame(react1_readable)

react1_df <- na.omit(react1_df)

write.csv(react1_df, file.path(OUTDIR_3, "gsea_Reactome.csv"), row.names = FALSE)



#Visualizing results 
#Top 3 gene sets with most positive NES 
react1_df %>%
  dplyr::slice_max(NES, n = 3)


#Plotting most positive NES gene set 
React_most_pos_nes_plot <- enrichplot::gseaplot(
  react1,
  geneSetID = "R-MMU-6809371",
  title = "Formation of the cornified envelope",
  color.line = "cadetblue2"
)

React_most_pos_nes_plot

#To save
ggplot2::ggsave(file.path(FIGDIR_3, "Reactome_gsea_enrich_positive_plot.png"), React_most_pos_nes_plot, width = 18, height = 15, dpi = 300)



#Top 3 gene sets with most negative NES
react1_df %>%
  dplyr::slice_min(NES, n = 3)


#Plotting most negative NES gene set 
React_most_neg_nes_plot <- enrichplot::gseaplot(
  react1,
  geneSetID = "R-MMU-211945",
  title = "Phase I - Functionalization of compounds",
  color.line = "cadetblue2"
)

React_most_neg_nes_plot

#To save
ggplot2::ggsave(file.path(FIGDIR_3, "Reactome_gsea_enrich_negative_plot.png"), React_most_neg_nes_plot, width = 18, height = 15, dpi = 300)



#Creating more complex plots

#cnetplot 

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)

#Top 10 significantly enriched pathways 
#Only labelled shared leading edge genes among pathways

#Creating genelist to colour code genes based on their log2FC

genelist <- RNA_data$log2FC_fib
names(genelist) <- RNA_data$gene_id

genelist <- na.omit(genelist)

p13 <- cnetplot(react1_readable,
         showCategory = 10, 
         foldChange = genelist,
         fc_threshold = 1,
         categorySizeBy = ~p.adjust, 
         color_category = 'skyblue',
         node_label = "share"
         #, color_item = 'steelblue' (only run if foldChange does not exist)
)  +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick2",
    midpoint = 0
  ) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

print(p13)

ggsave(file.path(FIGDIR_3, "Top_10_sig_enrichpathways.png"), p13, bg = "white", width = 20, height = 18, dpi = 300)



#Top 3 positively and negatively enriched pathways (plots only 1 pathway for positively enriched condition due to fc_threshold)
#Labelled all genes from each pathway
#Genelist applied

top_up <- react1_df %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(p.adjust) %>%
  head(3) %>%
  dplyr::pull(ID)

top_down <- react1_df %>%
  dplyr::filter(NES < 0) %>%
  dplyr::arrange(p.adjust) %>%
  head(3) %>%
  dplyr::pull(ID)

target_pathways <- c(top_up, top_down)

reinserted_react1 <- react1_readable
reinserted_react1@result <- react1_readable@result %>%
  dplyr::filter(ID %in% target_pathways)

p14 <- cnetplot(reinserted_react1, 
         showCategory = length(target_pathways),
         foldChange = genelist,
         fc_threshold = 1,
         categorySizeBy = ~p.adjust,
         color_category = 'skyblue'
) +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick2",
    midpoint = 0
  ) +
theme(plot.background = element_rect(fill = "white", color = NA))

print(p14)

ggsave(file.path(FIGDIR_3, "Top_1_pos_3_neg_enrichpathways.png"), p14, bg = "white", width = 20, height = 15, dpi = 300)




#Use networkD3 and RCytoscape for larger gene sets
library(networkD3)
library(jsonlite)
library(tidyverse)

#Template code
# edge_list <- BP_gsea_df %>%
#   dplyr::select(Description, core_enrichment) %>%
#   tidyr::separate_rows(core_enrichment, sep = "/") %>%
#   dplyr::rename(Source = core_enrichment, Target = Description)
# 
# nodes <- data.frame(
#   name = unique(c(edge_list$Source, edge_list$Target)),
#   stringsAsFactors = FALSE
#   )
#  
# nodes$group <- ifelse(nodes$name %in% edge_list$Source, "Gene", "Pathway")
# nodes$id <- 0:(nrow(nodes) - 1)
#  
# links <- edge_list %>%
#   dplyr::left_join(nodes, by = c("Source" = "name")) %>%
#   dplyr::rename(source_id = id) %>%
#   dplyr::select(Source, Target, source_id) %>%
#   dplyr::left_join(nodes, by = c("Target" = "name")) %>%
#   dplyr::rename(target_id = id)
#   dplyr::mutate(value_weight = 1)
# 
# networkD3_plot <- forceNetwork(
#   Links = links, 
#   Nodes = nodes,
#   Source = "source_id",
#   Target = "target_id", 
#   Value = "value_weight",
#   NodeID = "name", 
#   Group = "group",
#   opacity = 0.8,
#   zoom = TRUE,
#   legend = TRUE,
#   fontSize = 12,
#   charge = -120
#   )
# 
# #print plot
# networkD3_plot


#Interactive plot created by networkD3
#Separating pathways based on NES 
activated_pathways <- react1_df %>% 
  dplyr::filter(NES > 0) %>% 
  dplyr::pull(Description)

suppressed_pathways <- react1_df %>% 
  dplyr::filter(NES < 0) %>% 
  dplyr::pull(Description)

#Creating Node and Links dataframe  
edge_list <- react1_df %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(Source = core_enrichment, Target = Description)

nodes <- data.frame(
  name = unique(c(edge_list$Source, edge_list$Target)),
  stringsAsFactors = FALSE
)

nodes <- nodes %>%
  dplyr::mutate(group = case_when(
    name %in% activated_pathways ~ "Pathway (Activated)",
    name %in% suppressed_pathways ~ "Pathway (Suppressed)",
    TRUE ~ "Gene"
  ))

nodes$id <- 0:(nrow(nodes) - 1)

links <- edge_list %>%
  dplyr::left_join(nodes, by = c("Source" = "name")) %>%
  dplyr::rename(source_id = id) %>%
  dplyr::select(Source, Target, source_id) %>%
  dplyr::left_join(nodes, by = c("Target" = "name")) %>%
  dplyr::rename(target_id = id) %>%
  dplyr::mutate(value_weight = 1)

#Adding colour to legend
three_color_scale <- 'd3.scaleOrdinal()
  .domain(["Pathway (Activated)", "Pathway (Suppressed)", "Gene"])
  .range(["#D32F2F", "#1976D2", "#4A4A4A"]);'

#Assembling plot
networkD3_plot <- forceNetwork(
  Links = links, 
  Nodes = nodes,
  Source = "source_id",
  Target = "target_id", 
  Value = "value_weight",
  NodeID = "name", 
  Group = "group",
  colourScale = htmlwidgets::JS(three_color_scale),
  opacity = 0.9,
  zoom = TRUE,
  legend = TRUE,
  opacityNoHover = TRUE,
  fontSize = 12,
  charge = -1000,
  bounded = FALSE
)

networkD3_plot$x$options$cbCustomJS <- htmlwidgets::JS(
  "function(el, x) {
     d3.select(el).selectAll('.legend text').each(function() {
       var textElement = d3.select(this);
       var originalText = textElement.text();
       if (originalText === 'Pathway (Activated)') {
         textElement.text('Pathway'); 
       } else if (originalText === 'Pathway (Suppressed)') {
         d3.select(this.parentNode).style('display', 'none');
       }
     });
   }"
)

#print plot
networkD3_plot





library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)

#Heatmap of top 50 leading edge genes shared among the top 20 significantly enriched pathways
heatmap_1 <- heatplot(react1_readable, 
                      foldChange = genelist, 
                      showCategory = 20, 
                      showTop = 50)

plot(heatmap_1)

ggsave(file.path(FIGDIR_3, "Top_50_genes_shared_between_top_20_pathways.png"), heatmap_1, width = 20, height = 15, dpi = 300)



#Heatmap of top 100 leading edge genes shared among the top 20 significantly enriched pathways
heatmap_2 <- heatplot(react1_readable, 
                      foldChange = genelist, 
                      showCategory = 20, 
                      showTop = 100)

plot(heatmap_2)

ggsave(file.path(FIGDIR_3, "Top_100_genes_shared_between_top_20_pathways.png"), heatmap_2, width = 20, height = 15, dpi = 300)



#Heatmap of select pathways (and their associated leading edge genes) from the paper
select_path <- c("Glycolysis", "Keratinization", "Formation of the cornified envelope", "Fatty acid metabolism", "Proteasome assembly", "ZNF598 and the Ribosome−associated Quality Trigger (RQT) complex dissociate a ribosome stalled on a no−go mRNA")

heatmap_3 <- heatplot(react1_readable,
                      foldChange = genelist,
                      showCategory = select_path,
                      showTop = 100)

plot(heatmap_3)

ggsave(file.path(FIGDIR_3, "Select_pathways_from_paper_htmp_Reactome.png"), heatmap_3, width = 20, height = 10, dpi = 300)



#Heatmap of Glycolysis pathway (extracted leading edge genes - some of which match those found in paper)
ext_path <- c("Glycolysis")

heatmap_4 <- heatplot(react1_readable,
                      foldChange = genelist,
                      showCategory = ext_path)

plot(heatmap_4)

ggsave(file.path(FIGDIR_3, "Glycolysis_genes.png"), heatmap_4, width = 15, height = 5, dpi = 300)



#Creating a Log2FC bar plot of leading edge genes from the Glycolysis pathway

#Extracting genes from heatmap_4 and creating a vector
plot_h4 <- heatmap_4$data
head(plot_h4)

h4_genes <- unique(plot_h4$Gene)

#Extracting rows from RNA_seq dataset containing only genes from heatmap_4
clean_RNAdf <- na.omit(RNA_data)

h4_RNAdf <- clean_RNAdf %>%
  dplyr::filter(gene_id %in% h4_genes) %>%
  dplyr::arrange(gene_id)

#Plotting bar plot 
gly_plot <- ggplot(h4_RNAdf, aes(x = gene_id, y = log2FC_fib)) +
  geom_col(fill = "brown2", color = "black", linewidth = 0.2) + 
  geom_errorbar(
    aes(
      ymin = log2FC_fib - lfcSE,
      ymax = log2FC_fib + lfcSE
    ),
    width = 0.2,
    color = "black",
    linewidth = 0.2) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  geom_text(
    aes(
      y = log2FC_fib + lfcSE,
      label = as.character(stats::symnum(padj,
                                         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                         symbols = c("*\n*\n*", "*\n*", "*", "")))
    ),
    vjust = -0.3,
    lineheight = 0.4,
    hjust = 0.5,
    size = 5,
    color = "black") +
  scale_y_continuous(
    limits = c(0, 3),              
    breaks = seq(0, 3, by = 1),
    expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Glycolysis Pathway",
       x = "Gene", 
       y = expression(Log[2]~Fold~Change)) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, face = "italic"), 
    axis.text.y = element_text(color = "black", size = 10),
    panel.grid.major = element_blank())

print(gly_plot)

ggsave(file.path(FIGDIR_3, "Glycolysis_barplot.png"), gly_plot, width = 15, height = 5, dpi = 300)



#Ploting refined heatmap
#Extracting shared genes from heatmap_2 and creating a vector
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)

plot_data <- heatmap_2$data 

head(plot_data)

#Select genes shared in at least 2 pathways
sgl <- plot_data %>%
  dplyr::select(Gene) %>%
  tidyr::separate_rows(Gene, sep = "/") %>%
  dplyr::count(Gene, name = "pathway_count") %>%
  dplyr::filter(pathway_count >= 2) %>%
  dplyr::arrange(desc(pathway_count)) 

print(sgl)

hmp2_genes <- sgl %>%
  dplyr::pull(Gene) %>%
  unique()


#Extracting rows from RNA_seq dataset containing only shared genes from heatmap_2
clean_RNAdf <- na.omit(RNA_data)

hmp2_RNAdf <- clean_RNAdf %>%
  dplyr::filter(gene_id %in% hmp2_genes)


hmp2_RNAdf <- hmp2_RNAdf[order(hmp2_RNAdf$log2FC_fib, decreasing = TRUE),]


#Creating a matrix from filtered dataset
target_col <- ncol(hmp2_RNAdf) - 2

hmp2_numeric <- hmp2_RNAdf[, -target_col]
hmp2_matrix <- as.matrix(hmp2_numeric$log2FC_fib)

rownames(hmp2_matrix) <- hmp2_RNAdf$gene_id
colnames(hmp2_matrix) <- "log2FC"


#Checking if matrix is completely numeric
is.numeric(hmp2_matrix)


#Creating additional matrices to add onto heatmap
#padj column
padj_val <- as.matrix(hmp2_RNAdf$padj)
rownames(padj_val) <- hmp2_RNAdf$gene_id
colnames(padj_val) <- "padj"
mode(padj_val) <- "numeric"

#avg. expression column
mean <- as.matrix(hmp2_RNAdf$baseMean)
rownames(mean) <- hmp2_RNAdf$gene_id
colnames(mean) <- "AvgExpr"


#Create colour map for heatmap
#L2FC values
col_l2 <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))

#padj values
col_padj <- colorRamp2(c(min(padj_val), max(padj_val)), c("steelblue1", "white"))

#avg. expression values (quantiles 0, 25, 50, 75 from mean matrix)
col_AvgExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))


#Assembling heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = "steelblue3"),
                                               height = unit(2, "cm")))

h1 <- Heatmap(hmp2_matrix,
              name = "Log2FC",
              col = col_l2,
              cluster_rows = T, 
              cluster_columns = F,
              row_dend_width = unit(15, "cm"),
              width = unit(3, "cm")
)

h2 <- Heatmap(padj_val,
              name = "padj",
              top_annotation = ha,
              col = col_padj,
              cluster_rows = F, 
              width = unit(2, "cm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                p <- padj_val[i, j]
                
                # Add asterisks based on threshold
                if (p < 0.001) {
                  grid.text("***", x = x, y = y, gp = gpar(fontsize = 12, fontface = "bold"))
                } else if (p < 0.05) {
                  grid.text("*", x = x, y = y, gp = gpar(fontsize = 12, fontface = "bold"))
                }
              }
)

h3 <- Heatmap(mean,
              name = "AvgExpr",
              col = col_AvgExpr,
              cluster_rows = F, 
              show_row_names = T,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 10, fontface = "italic"),
              width = unit(2, "cm"),
              cell_fun = function(j, i, x, y, width, height, col) {
                grid.text(round(mean[i, j], 2), x, y, gp = gpar(fontsize = 8))
              })

h <- h1+h2+h3
h



##Alternative GSEA using pre-downloaded dataset
#library(clusterProfiler)
#library(fgsea)
#   
#REACT_gmt <- read.gmt("c2.cp.reactome.v7.0.mouse_symbols.gmt")
#  
#REACT_gsea_results <- GSEA(
#geneList     = final_react_vector,  
#TERM2GENE    = REACT_gmt,
#pvalueCutoff = 0.05,
#pAdjustMethod = "BH",
#minGSSize = 15,
#maxGSSize = 500)
#  
# 
#REACT_df <- data.frame(REACT_gsea_results)
#  
#  
##Plotting heatmaps
#heatmap_2 <- heatplot(REACT_gsea_results, 
#                      foldChange = genelist, 
#                      showCategory = 20, 
#                      showTop = 100)
# 
#plot(heatmap_2)