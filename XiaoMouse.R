#Paper written by Xiao et al. can be found with this link: https://www.nature.com/articles/s42255-018-0008-5

#Just a note about these mouse gene sets: I downloaded them separately and put them ALL in the same folder, along with this RStudio script.
#They're human gene sets that have been converted to mouse
#Then I set my WD to the same folder containing this script (In the options above, click session, set WD, To Source File Location)
#Below, make sure to adjust DE_File (the name of the dataset), GMT_Dir (the folder containing the mouse genesets), and an output for saving figures.
DE_FILE  <- "2024-05-16_differential_analysis_F20vN.txt"   
GMT_DIR  <- "genesets_mouse"                               
OUTDIR   <- "R_outputs_mouse"                              

#These are the signifciance values you can set
PADJ_CUT <- 0.05     
LFC_CUT  <- 1.25
#And if you want to switch whether the fibrotic case is the reference or the compared, you can adjust this line.
FIBROSIS_IS_NEGATIVE_LOG2FC <- TRUE


#Downloading stuff that you need
options(repos = c(CRAN = "https://cloud.r-project.org"))   # needed under source()

cran_pkgs <- c("ggplot2", "ggrepel", "patchwork")
bioc_pkgs <- c("fgsea")

loads     <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE))
need_cran <- cran_pkgs[!vapply(cran_pkgs, loads, logical(1))]
need_bioc <- bioc_pkgs[!vapply(bioc_pkgs, loads, logical(1))]
if (length(need_cran)) {
  message("Installing CRAN packages: ", paste(need_cran, collapse = ", "))
  install.packages(need_cran)
}
if (length(need_bioc)) {
  ## BiocManager is only needed if a Bioconductor package is actually missing.
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  message("Installing Bioconductor packages: ", paste(need_bioc, collapse = ", "))
  ## force = TRUE so a package's own missing dependencies get pulled in too.
  BiocManager::install(need_bioc, update = FALSE, ask = FALSE, force = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(fgsea)
})

dir.create(OUTDIR, showWarnings = FALSE)
FIGDIR <- file.path(OUTDIR, "figures"); dir.create(FIGDIR, showWarnings = FALSE)

## Fail early and clearly if inputs are missing.
if (!file.exists(DE_FILE))
  stop("Can't find DE file '", DE_FILE, "' in the working directory:\n  ", getwd(),
       "\n-> set the working directory (Session menu) or edit DE_FILE.")
gmt_files <- c(KEGG     = file.path(GMT_DIR, "c2.cp.kegg.v7.0.mouse_symbols.gmt"),
               REACTOME = file.path(GMT_DIR, "c2.cp.reactome.v7.0.mouse_symbols.gmt"),
               HALLMARK = file.path(GMT_DIR, "h.all.v7.0.mouse_symbols.gmt"))
missing_gmt <- gmt_files[!file.exists(gmt_files)]
if (length(missing_gmt))
  stop("Missing mouse gene-set file(s):\n  ", paste(missing_gmt, collapse = "\n  "),
       "\n-> keep the 'genesets_mouse' folder next to this script.")

#load the DE table
de <- read.csv(DE_FILE, stringsAsFactors = FALSE)
names(de)[names(de) == "gene_id"] <- "gene"

#I wasn't sure how to flip the fibrosis to be the negative/positive reference, so ChatGPT wrote this section for me
flip <- if (FIBROSIS_IS_NEGATIVE_LOG2FC) -1 else 1
de$stat_fib   <- flip * de$stat
de$log2FC_fib <- flip * de$log2FoldChange
message(sprintf("    %d genes total; %d with a valid adjusted p-value.",
                nrow(de), sum(!is.na(de$padj))))

tested   <- de[!is.na(de$padj) & !is.na(de$stat_fib), ]
tested_u <- tested[order(abs(tested$stat_fib), decreasing = TRUE), ]
tested_u <- tested_u[!duplicated(tested_u$gene), ]
ranks    <- sort(setNames(tested_u$stat_fib, tested_u$gene), decreasing = TRUE)


#reading the genesets from the folder
gene_sets <- lapply(gmt_files, fgsea::gmtPathways)   # $KEGG $REACTOME $HALLMARK
XiaoSet <- function(x)
  gsub("_", " ", sub("^(KEGG|REACTOME|HALLMARK)_", "", x))

#running GSEA
run_fgsea <- function(pathways, tag) {
  set.seed(42)
  res <- fgsea(pathways = pathways, stats = ranks, eps = 0,
               minSize = 15, maxSize = 500)
  res <- res[order(res$NES, decreasing = TRUE), ]
  out <- as.data.frame(res)
  out$leadingEdge <- vapply(out$leadingEdge, paste, collapse = ";", FUN.VALUE = "")
  write.csv(out, file.path(OUTDIR, sprintf("gsea_%s.csv", tag)), row.names = FALSE)
  message(sprintf("    %-9s: %d sets tested, %d significant (padj<0.05).",
                  tag, nrow(res), sum(res$padj < 0.05, na.rm = TRUE)))
  res
}
gsea <- list(KEGG     = run_fgsea(gene_sets$KEGG,     "KEGG"),
             REACTOME = run_fgsea(gene_sets$REACTOME, "REACTOME"),
             HALLMARK = run_fgsea(gene_sets$HALLMARK, "HALLMARK"))


#Drawing figure, similar to Fig. 1a
gsea_barplot <- function(res, tag, n = 14) {
  d <- as.data.frame(res); d <- d[!is.na(d$NES), ]
  keep <- d[d$padj < 0.25, ]
  if (nrow(keep) < 3) keep <- d
  keep <- head(keep[order(abs(keep$NES), decreasing = TRUE), ], n)
  keep$dir <- ifelse(keep$NES > 0, "Up in fibrosis", "Down in fibrosis")
  lev <- XiaoSet(keep$pathway[order(keep$NES)])
  keep$lab <- factor(XiaoSet(keep$pathway), levels = lev)
  ggplot(keep, aes(NES, lab, fill = dir)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Up in fibrosis" = "#c0392b",
                                 "Down in fibrosis" = "#2c6fbb")) +
    labs(x = "Normalized Enrichment Score",
         y = NULL, fill = NULL,
         title = sprintf("GSEA: %s pathways", tag)) +
    theme_classic(base_size = 11)
}
ggsave(file.path(FIGDIR, "gsea_KEGG_bar.png"),
       gsea_barplot(gsea$KEGG, "KEGG"),         width = 9, height = 6, dpi = 200)
ggsave(file.path(FIGDIR, "gsea_HALLMARK_bar.png"),
       gsea_barplot(gsea$HALLMARK, "Hallmark"), width = 9, height = 6, dpi = 200)

#Enrichment score curves.
key_paths <- c("KEGG_FATTY_ACID_METABOLISM", "KEGG_PPAR_SIGNALING_PATHWAY",
               "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "KEGG_RIBOSOME")
es_plot <- function(nm)
  fgsea::plotEnrichment(gene_sets$KEGG[[nm]], ranks) +
    labs(title = XiaoSet(nm)) +
    theme(plot.title = element_text(size = 9, face = "bold"))
es_grid <- (es_plot(key_paths[1]) | es_plot(key_paths[2])) /
  (es_plot(key_paths[3]) | es_plot(key_paths[4])) +
  plot_annotation(
    title    = "GSEA running enrichment score, key metabolic pathways"
    # subtitle = "leading edge genes marked prior to (positive enrichment) or following (negative enrichment) peak in ES magnitude"
    )
ggsave(file.path(FIGDIR, "gsea_enrichment_key.png"), es_grid,
       width = 10, height = 7, dpi = 200)

#Fig 1b and 4h are the over-representation analysis, so this is that
universe <- unique(tested$gene)
sig      <- tested[!is.na(tested$padj) & tested$padj < PADJ_CUT, ]
fib_up   <- unique(sig$gene[sig$log2FC_fib >=  LFC_CUT])   
fib_down <- unique(sig$gene[sig$log2FC_fib <= -LFC_CUT])   

ora <- function(query, sets, universe, min_set = 5) {
  query <- intersect(query, universe)      
  N <- length(universe); n <- length(query)
  res <- lapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], universe); K <- length(g)
    if (K < min_set) return(NULL)
    k <- length(intersect(query, g)); if (k < 2) return(NULL)
    data.frame(Term = nm, overlap = k, set_size = K, query_size = n, universe = N,
               pct_set = round(100 * k / K, 1),
               pval = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
               genes = paste(sort(intersect(query, g)), collapse = ";"),
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)               
  if (!is.null(res)) {
    res <- res[order(res$pval), ]
    res$FDR <- p.adjust(res$pval, method = "BH")
  }
  res
}
ora_res <- list()
for (coll in c("KEGG", "REACTOME"))
  for (dir in c("UP", "DOWN")) {
    q <- if (dir == "UP") fib_up else fib_down
    r <- ora(q, gene_sets[[coll]], universe)
    ora_res[[paste0(coll, "_", dir)]] <- r
    write.csv(r, file.path(OUTDIR, sprintf("ora_%s_fibrosis%s.csv", coll, dir)),
              row.names = FALSE)
  }

#Overrepresentation analysis bar graph
CAP <- 8   
ora_bardata <- function(direction, n = 6) {
  do.call(rbind, lapply(c("KEGG", "REACTOME"), function(coll) {
    r <- ora_res[[paste0(coll, "_", direction)]]
    if (is.null(r)) return(NULL)
    r <- r[r$pval < 0.05, ]; if (nrow(r) == 0) return(NULL)
    r <- head(r, n)                       # already sorted by p
    data.frame(collection = coll,
               term = paste0(XiaoSet(r$Term),
                             ifelse(direction == "DOWN", " ", "")),
               nlp  = pmin(-log10(r$pval), CAP),
               sig  = r$FDR < 0.05,
               direction = if (direction == "UP") "Up in fibrosis" else "Down in fibrosis",
               stringsAsFactors = FALSE)
  }))
}
ora_df <- rbind(ora_bardata("UP"), ora_bardata("DOWN"))
ora_df$direction <- factor(ora_df$direction,
                           levels = c("Up in fibrosis", "Down in fibrosis"))
ora_df$term <- factor(ora_df$term, levels = ora_df$term[order(ora_df$nlp)])

p_ora <- ggplot(ora_df, aes(nlp, term, fill = collection, alpha = sig)) +
  geom_col(width = 0.72, colour = "grey20") +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  facet_wrap(~ direction, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("KEGG" = "cornflowerblue", "REACTOME" = "darkorange1")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45),
                     breaks = c("TRUE", "FALSE"),
                     labels = c("FDR<0.05", "nominal p only"), name = NULL) +
  labs(x = sprintf("-log10(p)  (capped at %d)", CAP), y = NULL, fill = NULL,
       title = "Over-representation analysis") +
  theme_classic(base_size = 11)
ggsave(file.path(FIGDIR, "ora_bar.png"), p_ora, width = 9.5, height = 8, dpi = 200)
