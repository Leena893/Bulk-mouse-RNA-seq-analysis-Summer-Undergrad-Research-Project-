#Paper written by Xiao et al. can be found with this link: https://www.nature.com/articles/s42255-018-0008-5

#Set working directory (WD) to the same folder containing this script (In the options above, click session, set WD, To Source File Location)
#Make sure your folder contains the R script, appropriate dataset, and downloaded background genesets (if required)

OUTDIR_1 <- "R_outputs_GOBP&HALL"

dir.create(OUTDIR_1, showWarnings = FALSE)
FIGDIR_1 <- file.path(OUTDIR_1, "figures_GOBP&HALL"); dir.create(FIGDIR_1, showWarnings = FALSE)

library(ggplot2)
library(ggrepel)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(scales)

#Reading in dataset
fRNA_data <- read.csv("RNA-seq Prototype Data.csv", header = TRUE)

#Setting healthy state as baseline and fibrosis as experimental group 
FIBROSIS_IS_NEGATIVE_LOG2FC <- TRUE

flip <- if (FIBROSIS_IS_NEGATIVE_LOG2FC) -1 else 1
fRNA_data$stat_fib   <- flip * fRNA_data$stat
fRNA_data$log2FC_fib <- flip * fRNA_data$log2FoldChange
message(sprintf("    %d genes total; %d with a valid adjusted p-value.",
                nrow(fRNA_data), sum(!is.na(fRNA_data$padj))))

fRNA_data <- fRNA_data[!duplicated(fRNA_data$gene_id), ]

#Filter data
#Removing all grids with NA
fRNA_data <- fRNA_data %>% 
  dplyr::filter(!is.na(baseMean),
         !is.na(log2FoldChange), 
         !is.na(lfcSE), 
         !is.na(stat), 
         !is.na(pvalue), 
         !is.na(padj), 
         !is.na(gene_id), 
         !is.na(stat_fib),
         !is.na(log2FC_fib))



#Filtering significant genes
MA_data <- fRNA_data %>% 
  dplyr::mutate(significant = ifelse(padj < 0.05 & abs(log2FC_fib) > 1.25, "Yes", "No")) 
    
#MA plot
p1 <- ggplot(MA_data, aes(x = baseMean, y = log2FC_fib)) + 
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) +
  labs(title = "MA Plot", 
       x = "Average Expression (log10 scale) (A)", 
       y = "Log2 Fold Change (M)") + 
  theme_minimal() +
  scale_x_continuous(trans = "log10", labels = scales::label_comma())

print(p1)
 
ggsave(file.path(FIGDIR_1, "canvas_MA_plot.png"), p1,
       width = 14, height = 7, dpi = 300)



#Label genes for p2
genes_to_label_p2 <- MA_data %>% 
  dplyr::filter(significant == "Yes") %>% 
  dplyr::arrange(desc(abs(log2FC_fib))) %>% 
  dplyr::slice(1:20)
print(genes_to_label_p2)

#MA plot labelled
p2 <- ggplot(MA_data, aes(x = baseMean, y = log2FC_fib)) + 
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) +
  labs(title = "MA Plot", 
       x = "Average Expression (log10 scale) (A)", 
       y = "Log2 Fold Change (M)") + 
  ggrepel::geom_text_repel(data = genes_to_label_p2,
                           aes(label = gene_id),
                           size = 3, color = "black", 
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50") +
  theme_bw() +
  scale_x_continuous(trans = "log10", labels = scales::label_comma())

print(p2)

ggsave(file.path(FIGDIR_1, "labelled_MA_plot.png"), p2,
       width = 14, height = 7, dpi = 300)



#Updated MA Plots 
#Adding Significance Categories
fRNA_data$significance <- "Not Significant"
fRNA_data$significance[fRNA_data$padj < 0.05 & abs(fRNA_data$log2FC_fib) > 1.25] <- "Significant (padj < 0.05, |FC| > 1.25)"
fRNA_data$significance[fRNA_data$padj < 0.05 & abs(fRNA_data$log2FC_fib) <= 1.25] <- "Significant (padj < 0.05, |FC| ≤ 1.25)"
table(fRNA_data$significance)

#Specifying Upregulation or Downregulation
fRNA_data <- fRNA_data %>% 
  dplyr::mutate(significance_2 = case_when(
  padj < 0.05 & log2FC_fib > 1.25 ~ 'Significantly upregulated (padj < 0.05, FC > 1.25)',
  padj < 0.05 & log2FC_fib < -1.25 ~ 'Significantly downregulated (padj < 0.05, FC < -1.25)',
  .default = 'Not significant'
))



#Label genes for p3
genes_to_label_p3 <- fRNA_data %>% 
  dplyr::filter(significance == "Significant (padj < 0.05, |FC| > 1.25)") %>% 
  dplyr::arrange(desc(abs(log2FC_fib))) %>% 
  dplyr::slice(1:20)
print(genes_to_label_p3)

#MA plot sectioned based on padj
p3 <- ggplot(fRNA_data, aes(x = baseMean, y = log2FC_fib, color = significance)) + 
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c("Not Significant" = "grey",
               "Significant (padj < 0.05, |FC| ≤ 1.25)" = "lightblue",
               "Significant (padj < 0.05, |FC| > 1.25)" = "darkblue")
    ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_hline(yintercept = c(-1,1), linetype = "dotted", color = "blue", linewidth = 0.5) +
  labs(title = "MA Plot", 
       x = "Average Expression (log10 scale) (A)", 
       y = "Log2 Fold Change (M)", 
       color = "Differential Expression") + 
  ggrepel::geom_text_repel(data = genes_to_label_p3,
                           aes(label = gene_id),
                           size = 3, color = "black", 
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.title = element_text(size = 11), 
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_continuous(trans = "log10", labels = scales::label_comma())

print(p3)

ggsave(file.path(FIGDIR_1, "padj_MA_plot.png"), p3,
       width = 14, height = 7, dpi = 300)


 
#Label genes for p4
ugenes_to_label_p4 <- fRNA_data %>% 
  dplyr::filter(significance_2 == "Significantly upregulated (padj < 0.05, FC > 1.25)") %>% 
  dplyr::arrange(desc(abs(log2FC_fib))) %>% 
  dplyr::slice(1:10)
print(ugenes_to_label_p4)

dgenes_to_label_p4 <- fRNA_data %>% 
  dplyr::filter(significance_2 == "Significantly downregulated (padj < 0.05, FC < -1.25)") %>% 
  dplyr::arrange(desc(abs(log2FC_fib))) %>% 
  dplyr::slice(1:10)
print(dgenes_to_label_p4)

#Use regionReport for summary of results & plots

#MA plot sectioned based on upregulation or downregulation
p4 <- ggplot(fRNA_data, aes(x = baseMean, y = log2FC_fib, color = significance_2)) + 
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_hline(yintercept = c(-1,1), linetype = "dotted", color = "blue", linewidth = 0.5) +
  labs(title = "MA Plot", 
       x = "Average Expression (log10 scale) (A)", 
       y = "Log2 Fold Change (M)", 
       color = "Differential Expression") + 
  ggrepel::geom_text_repel(data = ugenes_to_label_p4,
                           aes(label = gene_id),
                           size = 3, color = "black", 
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50") +
  ggrepel::geom_text_repel(data = dgenes_to_label_p4,
                           aes(label = gene_id),
                           size = 3, color = "black", 
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50") +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.title = element_text(size = 11), 
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_continuous(trans = "log10", labels = scales::label_comma())

print(p4)

ggsave(file.path(FIGDIR_1, "UPDOWN_MA_plot.png"), p4,
       width = 14, height = 7, dpi = 300)



#Volcano plot
#Significant p-values
fRNA_data <- fRNA_data %>% 
  dplyr::mutate(significant_p = pvalue < 0.05)

#Label genes for p5
genes_to_label_p5 <- fRNA_data %>% 
  dplyr::filter(significant_p == "TRUE") %>% 
  dplyr::arrange(desc(abs(log2FC_fib) & abs(-log10(pvalue)))) %>% 
  dplyr::slice(1:20)
print(genes_to_label_p5)


p5 <- ggplot(fRNA_data, aes(x = log2FC_fib, y = -log10(pvalue))) +
  geom_point(aes(color = significant_p)) + 
  scale_color_manual(values = c("TRUE" = "steelblue1", "FALSE" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-13, 13), ylim = c(0,250)) +
  ggrepel::geom_text_repel(data = genes_to_label_p5,
                           aes(label = gene_id),
                           size = 3, color = "black", 
                           max.overlaps = Inf, box.padding = 0.4,
                           point.padding = 0.3, segment.color = "grey50") +
  theme_bw() +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~p~value),
    title = "Volcano Plot",
    color = "Significant (p < 0.05)"
  )
  
print(p5)
  
ggsave(file.path(FIGDIR_1, "Volcano_plot.png"), p5,
       width = 10, height = 14, dpi = 300)



#GSEA 
#Running GSEA using GO:BP (revised)
#Using clusterProfiler

library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)
library(enrichplot)
 

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


#Setting background genes 
msigdbr_species()

mm_ontology_sets <- msigdbr(
  species = "Mus musculus", 
  category = "C5", 
  subcollection = "GO:BP")

head(mm_ontology_sets)


#Setting pre-ranked gene list
#Always check if there are duplicated genes prior
lfc_vector <- RNA_data$stat_fib
names(lfc_vector) <- RNA_data$gene_id

lfc_vector <- na.omit(lfc_vector)

lfc_vector <- sort(lfc_vector, decreasing = TRUE)

head(lfc_vector)


#Mark for reproducible results 
set.seed(2026)


#Running GSEA
BP_gsea_results <- GSEA(
  geneList = lfc_vector, 
  minGSSize = 15,
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    mm_ontology_sets,
    gs_name,
    gene_symbol
  )
)

head(BP_gsea_results)

#Filtering out any row with NAs
BP_gsea_results@result <- BP_gsea_results@result[complete.cases(BP_gsea_results@result), ]


#Converting GSEA results to a data frame 
BP_gsea_df <- data.frame(BP_gsea_results)

BP_gsea_df <- na.omit(BP_gsea_df)

write.csv(BP_gsea_df, file.path(OUTDIR_1, "gsea_GOBP&HALL.csv"), row.names = FALSE)



#Visualizing results 
#Top 3 gene sets with most positive NES 
BP_gsea_df %>%
  dplyr::slice_max(NES, n = 3)


#Plotting most positive NES gene set 
most_positive_nes_plot <- enrichplot::gseaplot(
  BP_gsea_results,
  geneSetID = "GOBP_KERATINIZATION",
  title = "Keratinization",
  color.line = "cadetblue2"
)

most_positive_nes_plot

#To save
ggplot2::ggsave(file.path(FIGDIR_1, "MM_GOBP_gsea_enrich_positive_plot.png"), most_positive_nes_plot, width = 18, height = 15, dpi = 300)



#Top 3 gene sets with most negative NES
BP_gsea_df %>%
  dplyr::slice_min(NES, n = 3)


#Plotting most negative NES gene set 
most_negative_nes_plot <- enrichplot::gseaplot(
  BP_gsea_results,
  geneSetID = "GOBP_REGULATION_OF_ANIMAL_ORGAN_MORPHOGENESIS",
  title = "Regulation of animal organ morphogenesis",
  color.line = "cadetblue2"
)

most_negative_nes_plot

#To save
ggplot2::ggsave(file.path(FIGDIR_1, "MM_GOBP_gsea_enrich_negative_plot.png"), most_negative_nes_plot, width = 18, height = 15, dpi = 300)



#Creating more complex plots

#cnetplot 

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)


#Top 10 significantly enriched pathways (only plots 8 pathways due to fc_threshold) 
#Labelled all genes from each pathway

#Creating genelist to colour code genes based on their log2FC
genelist <- RNA_data$log2FC_fib
names(genelist) <- RNA_data$gene_id

genelist <- na.omit(genelist)
  
p6 <- cnetplot(BP_gsea_results,
         showCategory = 10, 
         foldChange = genelist,
         fc_threshold = 2,
         categorySizeBy = ~p.adjust, 
         color_category = 'skyblue'
         #, color_item = 'steelblue' (only run if foldChange does not exist)
         ) +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick2",
    midpoint = 0
  ) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

print(p6)

ggsave(file.path(FIGDIR_1, "Top_8_sig_enrichpathways.png"), p6, bg = "white", width = 20, height = 15, dpi = 300)



#Top 3 positively and negatively enriched pathways
#Included function for labeling only shared leading edge genes among pathways
#Genelist applied

top_up <- BP_gsea_df %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(p.adjust) %>%
  head(3) %>%
  dplyr::pull(ID)

top_down <- BP_gsea_df %>%
  dplyr::filter(NES < 0) %>%
  dplyr::arrange(p.adjust) %>%
  head(3) %>%
  dplyr::pull(ID)

target_pathways <- c(top_up, top_down)

reinserted_BP_gsea_results <- BP_gsea_results
reinserted_BP_gsea_results@result <- BP_gsea_results@result %>%
  dplyr::filter(ID %in% target_pathways)

p7 <- cnetplot(reinserted_BP_gsea_results, 
         showCategory = length(target_pathways),
         foldChange = genelist,
         fc_threshold = 2,
         categorySizeBy = ~p.adjust,
         color_category = 'skyblue'
         #, node_label = "share" (to label only share genes between pathways)
         ) +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick2",
    midpoint = 0
  ) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

print(p7)

ggsave(file.path(FIGDIR_1, "Top_3_pos_neg_enrichpathways.png"), p7, bg = "white", width = 20, height = 15, dpi = 300)


  
#Use networkD3 and RCytoscape for larger gene sets
library(clusterProfiler)
library(networkD3)
library(jsonlite)
library(dplyr)
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
activated_pathways <- BP_gsea_df %>% 
  dplyr::filter(NES > 0) %>% 
  dplyr::pull(Description)

suppressed_pathways <- BP_gsea_df %>% 
  dplyr::filter(NES < 0) %>% 
  dplyr::pull(Description)
  
#Creating Node and Links dataframe  
edge_list <- BP_gsea_df %>%
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



#Leading Edge Analysis (clusterProfiler does it automatically and gives it in the "core enrichment" column, but we want to seperate the genes now)
#To isolate individual genes into rows and create matrix (tot. 1963 genes)
# gene_list_per_pathway <- strsplit(BP_gsea_df$core_enrichment, "/")
# names(gene_list_per_pathway) <- BP_gsea_df$Description
# 
# all_shared_genes <- table(unlist(gene_list_per_pathway))
# shared_genes <- names(all_shared_genes[all_shared_genes > 1])
# 
# pathway_gene_matrix <- sapply(gene_list_per_pathway, function(pathway_genes) {
#   shared_genes %in% pathway_genes})
# 
# rownames(pathway_gene_matrix) <- shared_genes 
#  
# head(pathway_gene_matrix)



library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)

#Heatmap of top 50 leading edge genes shared among the top 20 significantly enriched pathways
heatmap_1 <- heatplot(BP_gsea_results, 
                      foldChange = genelist, 
                      showCategory = 20, 
                      showTop = 50)

plot(heatmap_1)

ggsave(file.path(FIGDIR_1, "Top_50_genes_shared_between_top_20_pathways.png"), heatmap_1, width = 20, height = 15, dpi = 300)



#Heatmap of top 100 leading edge genes shared among the top 20 significantly enriched pathways
heatmap_2 <- heatplot(BP_gsea_results, 
                      foldChange = genelist, 
                      showCategory = 20, 
                      showTop = 100)

plot(heatmap_2)

ggsave(file.path(FIGDIR_1, "Top_100_genes_shared_between_top_20_pathways.png"), heatmap_2, width = 20, height = 15, dpi = 300)



#Heatmap of leading edge genes from collagen metabolism 
select_path <- c("GOBP_COLLAGEN_METABOLIC_PROCESS")

heatmap_3 <- heatplot(BP_gsea_results,
                      foldChange = genelist,
                      showCategory = select_path)

plot(heatmap_3)

ggsave(file.path(FIGDIR_1, "Collagen_metabolism_genes.png"), heatmap_3, width = 15, height = 5, dpi = 300)



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



#Assessing shared leading edge genes
#Still learning to interpret it
#upsetplot(BP_gsea_results)



# ##Alternative GSEA using pre-downloaded dataset
# library(clusterProfiler)
# library(fgsea)
#   
# HALL_gmt <- read.gmt("h.all.v7.0.mouse_symbols.gmt")
#  
# HALL_gsea_results <- GSEA(
# geneList     = lfc_vector,  
# TERM2GENE    = HALL_gmt,
# pvalueCutoff = 0.05,
# pAdjustMethod = "BH",
# minGSSize = 15,
# maxGSSize = 500)
#    
# HALL_df <- data.frame(HALL_gsea_results)
# 
# 
# ##Plotting heatmaps
# heatmap_2 <- heatplot(HALL_gsea_results, 
#                       foldChange = genelist, 
#                       showCategory = 20, 
#                       showTop = 100)
#  
# plot(heatmap_2)



#Below this line is code that contains unexpected results and must be troubleshooted!
#Can use as template code and tweak as necessary
--------------------------------------------------------------------------------------------------------------



#GSEA 
#Running GSEA using GO:BP 
#Using fgsea
#Warning - results were unexpected (i.e., fibrotic state was set as a baseline reference instead of healthy state)

library(HDO.db)
library(data.table)
library(clusterProfiler)
library(tidyverse)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(ggplot2)


#Running GSEA using GO:BP
#Preparing background genes (BP)
gene_sets_df <- msigdbr(species = "Mus musculus", category = "C5", subcollection = "GO:BP")
head(gene_sets_df)

gene_sets <- gene_sets_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

cat(paste("Loaded", length(gene_sets), "gene sets\n"))


#Rank genes
#Creating a revised dataset
new_df <- fRNA_data %>% 
  select(gene_id, pvalue, padj, log2FoldChange)

fdf <- read.csv("Filtered RNA-seq Dataset.csv")
head(fdf)

rankings <- sign(fdf$log2FoldChange)*(-log10(fdf$pvalue))
names(rankings) <- fdf$gene_id
head(rankings)

rankings <- sort(rankings, decreasing = TRUE)
plot(rankings)


#Remember to filter out non-finite rankings
#max_ranking <- max(rankings[is.finite(rankings)])
#min_ranking <- min(rankings[is.finite(rankings)])
#rankings <- replace(rankings, rankings > max_ranking, max_ranking*10)
#rankings <- replace(rankings, rankings < min_ranking, min_ranking*10)
#rankings <- sort(rankings, decreasing = TRUE)


#Checking the gene rankings with ggplot
ggplot(data.frame(gene_ID = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_ID, ranks)) +
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#Running GSEA
GSEAres <- fgsea(pathways = gene_sets, 
                 stats = rankings,
                 scoreType = 'std',
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1)

head(GSEAres)



#Visualizing results 
#Top 6 enriched pathways (by p-value o/ p-adjusted value)
#Filtering out NAs (due to skewed data)

fGSEAres <- GSEAres %>% 
  filter(!is.na(pathway),
         !is.na(pval), 
         !is.na(padj), 
         !is.na(log2err), 
         !is.na(ES), 
         !is.na(NES), 
         !is.na(size), 
         !is.na(leadingEdge))

head(fGSEAres[order(pval), ])

sum(fGSEAres[, padj < 0.05])
sum(fGSEAres[, pval < 0.05])


#Plotting pathways
number_of_top_pathways_up <- 10
number_of_top_pathways_down <- 10

topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(gene_sets[topPathways], stats = rankings, fgseaRes = fGSEAres, gseaParam = 0.5)
#To export:
#pdf(file = paste0(filename, '_gsea_top20pathways.pdf'), width = 20, height = 15)
#dev.off()


#Selecting only independent pathways, removing redundancies/similar pathways (if required)
# #Warning was issued
# collapsedPathways <- collapsePathways(fGSEAres[order(pval)][pval < 0.05], gene_sets, rankings)
# mainPathways <- fGSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
# 
# plotGseaTable(gene_sets[mainPathways], rankings, fGSEAres, gseaParam = 0.5)
# #pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
# #dev.off()


#Plotting ES of most significantly enriched pathway (+) (ordered via padj)
plotEnrichment(gene_sets[[head(fGSEAres[order(padj), ], 1)$pathway]],
               rankings) +
  labs(title = head(fGSEAres[order(padj), ], 1)$pathway)


#Cleaning up the plot
p8 <- plotEnrichment(gene_sets[['GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS']],
                     rankings) +
  labs(title = 'Biological pathway: Regulation of cellular response to growth factor stimulus') +
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) + 
  scale_y_continuous('Enrichemnt score (ES)') +
  geom_line(aes(x = rank, y = ES), col = 'cadetblue3', linewidth = 1)

print(p8)


#Plotting ES of most significantly enriched pathway (-) (ordered via padj)
p9 <- plotEnrichment(gene_sets[['GOBP_NEGATIVE_REGULATION_OF_RELEASE_OF_CYTOCHROME_C_FROM_MITOCHONDRIA']],
                     rankings) +
  labs(title = 'Biological pathway: Negative regulation of cytochrome C release from mitochondria') +
  theme_classic() +
  scale_x_continuous('Rank', breaks = seq(0, 32000, 5000)) + 
  scale_y_continuous('Enrichemnt score (ES)') +
  geom_line(aes(x = rank, y = ES), col = 'plum3', linewidth = 1)

print(p9)


#Extracting leading edge genes
fGSEA_le <- fGSEAres %>% 
  select(pathway, leadingEdge)

head(fGSEA_le)


# #Saving results ()
# out_path <- "GSEA Analyses/"
# name_of_comparison <- 'fibroticvshealthy'
# background_genes <- 'gobp'
# filename <- paste0(out_path, 'GSEA/', name_of_comparison, '_', background_genes)
# saveRDS(fGSEAres, file = paste0(filename, '_gsea_results.RDS'))
# data.table::fwrite(fGSEAres, file = paste0(filename, '_gsea_results.tsv'), sep = "\t", sep2 = c("", " ", ""))


#Heatmap of diferentially expressed genes and gene sets
#Warning - troubleshooting heatmap, very messy
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)


#filter out significant gene sets
fGSEA_le <- fGSEAres %>% 
  filter(pval < 0.01)

head(fGSEA_le)


#Unnest leading edge and create a binary indicator
sig_pathway_df <- fGSEA_le %>% 
  select(pathway, leadingEdge) %>% 
  unnest(leadingEdge) %>% 
  mutate(present = 1) %>% 
  distinct()

head(sig_pathway_df)


#Convert data to a matrix
sig_pathway_mat <- sig_pathway_df %>% 
  pivot_wider(names_from = pathway, values_from = present, values_fill = 0) %>% 
  column_to_rownames("leadingEdge") %>% 
  as.matrix()

head(sig_pathway_mat)


#Basic heatmap
Basic_heatmap <- Heatmap(
  sig_pathway_mat, 
  col = c("white", "darkblue"), 
  name = "In leading edge",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  column_title = "Genes Expressed in Signicant Pathways",
  row_names_gp = gpar(fontsize = 3),
  )

draw(Basic_heatmap)





#Running GSEA using GO:MF
#Using fgsea
#Warning - results were unexpected (i.e., fibrotic state was set as a baseline reference instead of healthy state)

#Preparing background genes (MF)
MF_gene_sets_df <- msigdbr(species = "Mus musculus", category = "C5", subcollection = "GO:MF")
head(MF_gene_sets_df)

MF_gene_sets <- MF_gene_sets_df %>% 
  split(x = .$gene_symbol, f = .$gs_name)

cat(paste("Loaded", length(MF_gene_sets), "MF_gene sets\n"))


#Rank genes
#Creating a revised dataset
new_df <- fRNA_data %>% 
  select(gene_id, pvalue, padj, log2FoldChange)

fdf <- read.csv("Filtered RNA-seq Dataset.csv")
head(fdf)

rankings <- sign(fdf$log2FoldChange)*(-log10(fdf$pvalue))
names(rankings) <- fdf$gene_id
head(rankings)

rankings <- sort(rankings, decreasing = TRUE)
plot(rankings)


#Remember to filter out non-finite rankings
#max_ranking <- max(rankings[is.finite(rankings)])
#min_ranking <- min(rankings[is.finite(rankings)])
#rankings <- replace(rankings, rankings > max_ranking, max_ranking*10)
#rankings <- replace(rankings, rankings < min_ranking, min_ranking*10)
#rankings <- sort(rankings, decreasing = TRUE)


#Checking the gene rankings with ggplot
ggplot(data.frame(gene_ID = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_ID, ranks)) +
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#Running GSEA
MF_GSEAres <- fgsea(pathways = MF_gene_sets, 
                 stats = rankings,
                 scoreType = 'std',
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1)

head(MF_GSEAres)



#Visualizing results 
#Top 6 enriched pathways (by p-value o/ p-adjusted value)
#Filtering out NAs (due to skewed data)

MF_fGSEAres <- MF_GSEAres %>% 
  filter(!is.na(pathway),
         !is.na(pval), 
         !is.na(padj), 
         !is.na(log2err), 
         !is.na(ES), 
         !is.na(NES), 
         !is.na(size), 
         !is.na(leadingEdge))

head(MF_fGSEAres[order(pval), ])

sum(MF_fGSEAres[, padj < 0.05])
sum(MF_fGSEAres[, pval < 0.05])


#Plotting pathways
number_of_top_pathways_up <- 10
number_of_top_pathways_down <- 10

MF_topPathwaysUp <- MF_GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
MF_topPathwaysDown <- MF_GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
MF_topPathways <- c(MF_topPathwaysUp, rev(MF_topPathwaysDown))

plotGseaTable(MF_gene_sets[MF_topPathways], stats = rankings, fgseaRes = MF_fGSEAres, gseaParam = 0.5)



#Filter out significant Gene Sets 
ref_MF_fGSEAres <- MF_GSEAres %>% 
  filter(pval < 0.05)


#NES plot of significant gene sets
p10 <- ggplot(data = ref_MF_fGSEAres) +
  geom_bar(aes(x = pathway, y = NES, fill = pval), stat = 'identity') + 
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue") + 
  xlab("Molecular Function") +
  ylab("Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) + 
  scale_y_continuous(expand = c(0,0.25)) + 
  scale_x_discrete(expand = c(0,1)) +
  theme_bw()

print(p10)
