# Load libraries
library(dplyr)
library(readr)
library(tibble)


# List all files ending in .txt
files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)

# Read the first file to start
counts_merged <- read_tsv(files[1], col_names = c("gene_id", "gene_name", tools::file_path_sans_ext(basename(files[1]))))

# Loop through the rest and merge
for (f in files[-1]) {
  sample_name <- tools::file_path_sans_ext(basename(f))
  temp <- read_tsv(f, col_names = c("gene_id", "gene_name", sample_name))
  counts_merged <- full_join(counts_merged, temp, by = c("gene_id", "gene_name"))
}

# Check result
head(counts_merged)

library(DESeq2)

library(dplyr)
library(tibble)

count_matrix <- counts_merged %>%
  dplyr::select(-gene_name) %>%
  column_to_rownames(var = "gene_id")

coldata <- data.frame(
  sample = colnames(count_matrix),
  condition = c("Senescent", "Senescent", "Senescent", "Non-senscent", "Non-senscent", "Non-senscent", "Conditioned_media", "Conditioned_media")
)
rownames(coldata) <- coldata$sample

# 4. Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)

# sanity checks (optional but helpful)
colnames(colData(dds))              # should include "condition"
levels(colData(dds)$condition)      # should include "Senescent" and "Non-senscent"
table(colData(dds)$condition)       # counts per group

# correct contrast: Mutant vs Wild
# Minimal fix using your existing level names
res <- results(dds, contrast = c("condition", "Senescent", "Non-senscent"))
head(res)


# Convert results to data frame and keep gene_id
res_df <- as.data.frame(res_SvN) %>%
  mutate(gene_id = rownames(res_SvN))

# Volcano plot with gene names for significant points
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1),
                  aes(label = gene_name),
                  size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "Volcano Plot Analysis Highlights Differential Gene Expression Landscape of Senescent Versus Non-Senescent Cells",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value")

print(ggplot)


# install.packages(c("BiocManager"))
# BiocManager::install(c("clusterProfiler","org.Hs.eg.db","enrichplot","ReactomePA"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(dplyr)
library(readr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
# Install if needed
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)

#Split into up/down
up_genes <- sig %>%
  filter(log2FoldChange > 0) %>%
  rownames()   # rownames of res_SvN usually hold Ensembl IDs

down_genes <- sig %>%
  filter(log2FoldChange < 0) %>%
  rownames()

print (up_genes)
print (down_genes)

#Convert Ensembl IDs → Entrez IDs & Symbols
conv <- bitr(c(up_genes, down_genes),
             fromType = "ENSEMBL",
             toType = c("ENTREZID", "SYMBOL"),
             OrgDb = org.Hs.eg.db)


#Install if you haven't already
if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ReactomePA")
}

# Load it
library(ReactomePA)
library(stringr)
library(DESeq2)
library(readr)
library(dplyr)

##### ===== Exploratory enrichment with relaxed thresholds ===== #####


if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (p in c("clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot")) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE)
}
for (p in c("dplyr","tibble","stringr","ggplot2")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

## --- Parameters (relaxed) ---
padj_cutoff     <- 0.10   # relaxed from 0.05
pvalue_cutoff   <- 0.05   # fallback if padj yields too few genes
lfc_cutoff      <- 0      # include all effect sizes for exploration
min_genes_needed <- 10    # minimum to run enrichment
ont_terms       <- "BP"
show_n          <- 20

## --- Prep ---
res_df <- res_SvN %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL_raw") %>%
  mutate(ENSEMBL = sub("\\.\\d+$", "", ENSEMBL_raw))

universe_df <- res_df %>%
  filter(!is.na(padj)) %>%
  distinct(ENSEMBL, .keep_all=TRUE)

## --- Try adjusted P first, then raw P if needed ---
sig_padj <- res_df %>%
  filter(!is.na(padj), padj < padj_cutoff, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cutoff)

sig_use <- sig_padj
used_metric <- "padj"

if (nrow(sig_use) < min_genes_needed) {
  message("padj yielded ", nrow(sig_use), " genes; falling back to raw p-value < ", pvalue_cutoff)
  sig_p <- res_df %>%
    filter(!is.na(pvalue), pvalue < pvalue_cutoff, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cutoff)
  sig_use <- sig_p
  used_metric <- "pvalue"
}

message("Using ", used_metric, " with ", nrow(sig_use), " genes.")

## --- Split up/down ---
up_df   <- sig_use %>% filter(log2FoldChange > 0)
down_df <- sig_use %>% filter(log2FoldChange < 0)

## --- Map Ensembl → ENTREZ/SYMBOL ---
conv_univ <- bitr(universe_df$ENSEMBL, fromType="ENSEMBL",
                  toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

universe_mapped <- universe_df %>%
  left_join(conv_univ, by=c("ENSEMBL"="ENSEMBL")) %>%
  filter(!is.na(ENTREZID)) %>% arrange(padj) %>% distinct(ENSEMBL, .keep_all=TRUE)

map_and_dedup <- function(df) {
  df %>% left_join(conv_univ, by=c("ENSEMBL"="ENSEMBL")) %>%
    filter(!is.na(ENTREZID)) %>% arrange(padj) %>% distinct(ENSEMBL, .keep_all=TRUE)
}
up_mapped   <- map_and_dedup(up_df)
down_mapped <- map_and_dedup(down_df)

genes_up    <- up_mapped$ENTREZID
genes_down  <- down_mapped$ENTREZID
universe_bg <- universe_mapped$ENTREZID

cat("Counts:\n up_mapped:", length(genes_up),
    "\n down_mapped:", length(genes_down),
    "\n universe:", length(universe_bg), "\n")

## --- Enrichment (guard against empties) ---
run_GO <- function(genes) if (length(genes) >= min_genes_needed) enrichGO(
  gene=genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont_terms,
  universe=universe_bg, pAdjustMethod="BH", qvalueCutoff=0.25, readable=TRUE) else NULL  # relaxed q

run_KEGG <- function(genes) if (length(genes) >= min_genes_needed) enrichKEGG(
  gene=genes, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.25) else NULL

run_REACT <- function(genes) if (length(genes) >= min_genes_needed) enrichPathway(
  gene=genes, organism="human", pAdjustMethod="BH", qvalueCutoff=0.25, readable=TRUE) else NULL

ego_up   <- run_GO(genes_up);    ego_down  <- run_GO(genes_down)
ek_up    <- run_KEGG(genes_up);  ek_down   <- run_KEGG(genes_down)
er_up    <- run_REACT(genes_up); er_down   <- run_REACT(genes_down)

## --- Print top terms (if any) ---
show_top <- function(x, label) if (!is.null(x) && nrow(as.data.frame(x))) {
  cat("\nTop ", label, ":\n", sep="")
  print(head(as.data.frame(x)[, c("ID","Description","GeneRatio","p.adjust")], 15))
} else cat("\nNo significant ", label, " (given current relaxed settings).\n", sep="")

show_top(ego_up,  "GO BP (UP)")
show_top(ego_down,"GO BP (DOWN)")
show_top(ek_up,   "KEGG (UP)")
show_top(ek_down, "KEGG (DOWN)")
show_top(er_up,   "Reactome (UP)")
show_top(er_down, "Reactome (DOWN)")

## --- Plots only if there are results ---
if (!is.null(ego_up) && nrow(as.data.frame(ego_up)))  print(dotplot(ego_up,  showCategory=show_n, title="GO BP — UP"))
if (!is.null(ego_down) && nrow(as.data.frame(ego_down))) print(dotplot(ego_down,showCategory=show_n, title="GO BP — DOWN"))
if (!is.null(ek_up) && nrow(as.data.frame(ek_up)))    print(dotplot(ek_up,   showCategory=show_n, title="KEGG — UP"))
if (!is.null(ek_down) && nrow(as.data.frame(ek_down)))print(dotplot(ek_down, showCategory=show_n, title="KEGG — DOWN"))
if (!is.null(er_up) && nrow(as.data.frame(er_up)))    print(dotplot(er_up,   showCategory=show_n, title="Reactome — UP"))
if (!is.null(er_down) && nrow(as.data.frame(er_down)))print(dotplot(er_down, showCategory=show_n, title="Reactome — DOWN"))
##### ============================================================ #####

dotplot(er_down, showCategory = 10, title = "Downregulated Reactome Pathways In Senescent Cells") +
  theme(axis.text.y = element_text(size = 10))
library(enrichplot)  # if not already loaded
library(ggplot2)

p <- dotplot(
  er_down,
  showCategory = 10,
  title = "Downregulated Reactome Pathways In Senescent Cells"
) + theme(axis.text.y = element_text(size = 10))

print(p)  
ggsave("reactome_up_barplot_circRNA_finder.pdf", p_bar_up,
       width = 7, height = 6, device = cairo_pdf)

dotplot(up, showCategory = 10, title = "Upregulated Reactome Pathways In Senescent Cells") +
  theme(axis.text.y = element_text(size = 10))
library(enrichplot)  # if not already loaded
library(ggplot2)

q <- dotplot(
  er_up,
  showCategory = 10,
  title = "Upregulated Reactome Pathways In Senescent Cells"
) + theme(axis.text.y = element_text(size = 10))

print(p)  




















