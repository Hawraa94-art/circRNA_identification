install.packages("readr")  # Only needed once
library(readr)             # Load the package each session
# Load tool1, sample 1
tool1_sample1 <- read_tsv("SRR27937869_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample1) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample1$tool <- "Tool1"
tool1_sample1$sample <- "Melanocytes were transduced with mutant BRAFV600E and treated with CM from keratinocytes expressing shCTL1"
tool1_sample1$circ_id <- paste0(tool1_sample1$chrom, ":", tool1_sample1$start, "|", tool1_sample1$end)

# Load tool1, sample 2
tool1_sample2 <- read_tsv("SRR27937870_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample2) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample2$tool <- "Tool1"
tool1_sample2$sample <- "Melanocytes were transduced with mutant BRAFV600E and treated with CM from keratinocytes expressing shCTL"
tool1_sample2$circ_id <- paste0(tool1_sample2$chrom, ":", tool1_sample2$start, "|", tool1_sample2$end)


# Load tool1, sample 3
tool1_sample3 <- read_tsv("SRR27937871_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample3) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample3$tool <- "Tool1"
tool1_sample3$sample <- "Melanocytes were transduced with mutant BRAFV600E and treated with CM from keratinocytes expressing shCTL"
tool1_sample3$circ_id <- paste0(tool1_sample3$chrom, ":", tool1_sample3$start, "|", tool1_sample3$end)

# Load tool1, sample 4
tool1_sample4 <- read_tsv("SRR27937874_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample4) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample4$tool <- "Tool1"
tool1_sample4$sample <- "Melanocytes were transduced with wild type BRAF and treated with CM from keratinocytes expressing shDsg1"
tool1_sample4$circ_id <- paste0(tool1_sample4$chrom, ":", tool1_sample4$start, "|", tool1_sample4$end)

# Load tool1, sample 5
tool1_sample5 <- read_tsv("SRR27937875_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample5) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample5$tool <- "Tool1"
tool1_sample5$sample <- "Melanocytes were transduced with wild type BRAF and treated with CM from keratinocytes expressing shDsg1"
tool1_sample5$circ_id <- paste0(tool1_sample5$chrom, ":", tool1_sample5$start, "|", tool1_sample5$end)

# Load tool1, sample 6
tool1_sample6 <- read_tsv("SRR27937876_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample6) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample6$tool <- "Tool1"
tool1_sample6$sample <- "Melanocytes were transduced with wild type BRAF and treated with CM from keratinocytes expressing shDsg1"
tool1_sample6$circ_id <- paste0(tool1_sample6$chrom, ":", tool1_sample6$start, "|", tool1_sample6$end)

# Load tool1, sample 7
tool1_sample7 <- read_tsv("SRR27937879_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample7) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample7$tool <- "Tool1"
tool1_sample7$sample <- "Conditioned media from keratinocytes expressing shDsg1"
tool1_sample7$circ_id <- paste0(tool1_sample7$chrom, ":", tool1_sample7$start, "|", tool1_sample7$end)

# Load tool1, sample 8
tool1_sample8 <- read_tsv("SRR27937880_pass_1filteredJunctions.bed", col_names = FALSE)
colnames(tool1_sample8) <- c("chrom", "start", "end", "name", "score", "strand")
tool1_sample8$tool <- "Tool1"
tool1_sample8$sample <- "Conditioned media from keratinocytes expressing shDsg1"
tool1_sample8$circ_id <- paste0(tool1_sample8$chrom, ":", tool1_sample8$start, "|", tool1_sample8$end)

#Step 2 count and compare detected circRNAs
# Count unique circRNAs for Tool1, Sample1
length(unique(tool1_sample1$circ_id))

# Count unique circRNAs for Tool1, Sample2
length(unique(tool1_sample2$circ_id))

# Count unique circRNAs for Tool1, Sample3
length(unique(tool1_sample3$circ_id))

# Count unique circRNAs for Tool1, Sample4
length(unique(tool1_sample4$circ_id))

# Count unique circRNAs for Tool1, Sample5
length(unique(tool1_sample5$circ_id))

# Count unique circRNAs for Tool1, Sample6
length(unique(tool1_sample6$circ_id))

# Count unique circRNAs for Tool1, Sample7
length(unique(tool1_sample7$circ_id))

# Count unique circRNAs for Tool1, Sample8
length(unique(tool1_sample8$circ_id))

#Combine all ids

all_ids <- unique(c(tool1_sample1$circ_id, tool1_sample2$circ_id, tool1_sample3$circ_id, tool1_sample4$circ_id, tool1_sample5$circ_id, tool1_sample6$circ_id, tool1_sample7$circ_id, tool1_sample8$circ_id))
length(all_ids)  # Total non-redundant circRNAs detected across tools

#Step 3 overlap analysis
# Example: List of circRNA data frames for 8 samples (from Tool1)
sample_names <- list(
  Senescent_cells1 = tool1_sample1,
  Senescent_cells2 = tool1_sample2,
  Senescent_cells3 = tool1_sample3,
  Non_Senescent1 = tool1_sample4,
  Non_Senescent2 = tool1_sample5,
  Non_Senescent3 = tool1_sample6,
  Conditioned_media1 = tool1_sample7,
  Conditioned_media2 = tool1_sample8
)

print(sample_names)

## 2) Extract the circRNA ID vector from each sample
id_sets <- lapply(sample_names, function(d) unique(as.character(d$circ_id)))

## 3) Build the pairwise shared-count matrix
n <- length(id_sets)
shared_matrix <- matrix(0, n, n, dimnames = list(names(id_sets), names(id_sets)))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    shared_matrix[i, j] <- length(intersect(id_sets[[i]], id_sets[[j]]))
  }
}

## 4) (optional) as data frame for viewing
shared_df <- as.data.frame(shared_matrix)
print(shared_df)



#UpSet Plot
# Install required package
install.packages("UpSetR")
library(UpSetR)

# Build presence-absence data frame
library(dplyr)

# Create list of unique circ_ids for each sample
circ_lists <- list(
  Senescent_cells1 = unique(tool1_sample1$circ_id),
  Senescent_cells2 = unique(tool1_sample2$circ_id),
  Senescent_cells3 = unique(tool1_sample3$circ_id),
  Non_Senescent1 = unique(tool1_sample4$circ_id),
  Non_Senescent2 = unique(tool1_sample5$circ_id),
  Non_Senescent3 = unique(tool1_sample6$circ_id),
  Conditioned_media = unique(tool1_sample7$circ_id),
  Conditioned_media = unique(tool1_sample8$circ_id)
)
print(circ_lists)
install.packages("writexl")
library(writexl)

# Convert each vector to a data frame
circ_lists_df <- lapply(circ_lists, function(x) data.frame(circ_id = x))

# Export: one sheet per sample
write_xlsx(circ_lists_df, path = "~/Desktop/circ_lists.xlsx")

# Check names
names(circ_lists)

# Fix duplicate names (example: add "_rep1", "_rep2")
names(circ_lists) <- make.unique(names(circ_lists))

# Convert vectors -> data frames
circ_lists_df <- lapply(circ_lists, function(x) data.frame(circ_id = x))

# Export to Excel
write_xlsx(circ_lists_df, path = "~/Desktop/circ_lists.xlsx")


# Make a binary matrix
all_circs <- unique(unlist(circ_lists))

print (all_circs)

# If all_circs might be a list, flatten + de-dup:
all_circs_vec <- unique(as.character(unlist(all_circs)))

# TXT: one circRNA per line
writeLines(all_circs_vec, "~/Desktop/all_circs.txt")



binary_matrix <- data.frame(circ_id = all_circs)
print(binary_matrix)

for (sample in names(circ_lists)) {
  binary_matrix[[sample]] <- binary_matrix$circ_id %in% circ_lists[[sample]]
}

# Convert to standard data.frame (not tibble) and drop circ_id column
binary_matrix_df <- as.data.frame(binary_matrix[ , -1])
rownames(binary_matrix_df) <- binary_matrix$circ_id  # Optional: set rownames
print(binary_matrix_df)

# Force all columns to be 0/1
binary_matrix_df[] <- lapply(binary_matrix_df, function(x) as.numeric(as.logical(x)))

# Recheck structure
str(binary_matrix_df)

install.packages("ComplexUpset")

library(ComplexUpset)
library(ggplot2)

# Ensure binary_matrix_df is in correct format: rows = circRNAs, columns = tools, values = 0/1
# Column names should be the tool names


# Override default bar color
update_geom_defaults("bar", list(fill = "DarkBlue"))

upset(
  binary_matrix_df,
  intersect = colnames(binary_matrix_df),
  name = "Detected by CircRNA_Finder Tool",
  width_ratio = 0.2,
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3)
    )
  ),
  matrix = intersection_matrix(
    geom = geom_point(shape = 21, size = 3, fill = "#d62728", color = "#d62728")
  )
) + 
  ggtitle("Shared and Unique circRNAs Across BRAFV600E Conditions as Detected by CircRNA_Finder Tool") +
  scale_fill_manual(values = c(
    'Senescent_cells' = '#1f77b4',
    'Non_Senescent' = '#2ca02c',
    'Wild_BRAFV600E' = '#d62728'
  ))

# Override default bar color
update_geom_defaults("bar", list(fill = "Dark_Blue"))

print (upset)
library(ComplexUpset)
library(ggplot2)


# Filter to only show intersections of at least 10 circRNAs

library(UpSetR)

# Only show top 30 intersections
library(ComplexUpset)

# Override default bar color
update_geom_defaults("bar", list(fill = "DarkBlue"))

upset(
  binary_matrix_df,
  intersect = colnames(binary_matrix_df),
  min_size = 100,   # show only intersections with >= 100 circRNAs
  name = "Detected by CircRNA_Finder Tool",
  matrix = intersection_matrix(
    geom = geom_point(shape = 21, size = 3, fill = "#d62728", color = "#d62728")
  )
) + 
  ggtitle("Top 30 Intersections of circRNAs (≥100)")



library(ComplexUpset)

# Map sample columns to conditions (adjust patterns/names as needed)
cond_map <- list(
  Senescent_cells   = grep("^Senescent_cells",  colnames(binary_matrix_df), value = TRUE),
  Non_Senescent     = grep("^Non_Senescent",    colnames(binary_matrix_df), value = TRUE),
  Conditioned_media  = grep("^Conditioned_media",colnames(binary_matrix_df), value = TRUE)
)

# If matrix has 0/1 or counts, treat >0 as "detected"
to_logic <- function(x) { (as.matrix(x) > 0) }

# Union across replicates (present in ANY replicate of the condition)
collapsed_any <- sapply(cond_map, function(cols) rowSums(to_logic(binary_matrix_df[, cols, drop=FALSE])) >= 1)
collapsed_any <- as.data.frame(collapsed_any)

# (Optional) Intersection across replicates (present in ALL replicates)
collapsed_all <- sapply(cond_map, function(cols) rowSums(to_logic(binary_matrix_df[, cols, drop=FALSE])) == length(cols))
collapsed_all <- as.data.frame(collapsed_all)


upset(
  collapsed_any,                               # or collapsed_all
  intersect   = colnames(collapsed_any),
  min_size    = 100,                            # keep your ≥100 filter
  name        = "Detected by CircRNA_Finder Tool",
  matrix      = intersection_matrix(
    geom = geom_point(shape = 21, size = 3,
                      fill = "#d62728", color = "#d62728")
  )
) +
  ggtitle("Top Intersections of circRNAs (≥100) across three conditions")




upset(
  collapsed_any,                               # or collapsed_all
  intersect   = colnames(collapsed_any),
  min_size    = 100,
  name        = "Detected by CircRNA_Finder Tool",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 6)   # numbers above vertical bars
    )
  ),
  matrix      = intersection_matrix(
    geom = geom_point(shape = 21, size = 4,   # bigger dots
                      fill = "#d62728", color = "#d62728")
  )
) +
  ggtitle("Top Intersections of CircRNAs (≥100) Across Three Conditions") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 14),    # condition labels (set names)
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x  = element_blank(),   # hide cluttered index labels
    axis.ticks.x = element_blank(),
    # this controls the numbers on the horizontal "set size" bars:
    strip.text.y = element_text(size = 12)
  )




#Assess detection consistency between samples

jaccard_matrix <- matrix(0, nrow = length(sample_names), ncol = length(sample_names),
                         dimnames = list(sample_names, sample_names))
print(jaccard_matrix)
for (i in 1:length(sample_names)) {
  for (j in 1:length(sample_names)) {
    inter <- length(intersect(circ_lists[[i]], circ_lists[[j]]))
    uni <- length(union(circ_lists[[i]], circ_lists[[j]]))
    jaccard_matrix[i, j] <- inter / uni
  }
}

print(round(jaccard_matrix, 3))

# install.packages("pheatmap") # if not installed
library(pheatmap)

pheatmap(jaccard_matrix,
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Jaccard Similarity Coefficients Between
         circRNA Profiles in Melanocytes Transduced 
         Senescent or Wild-Type BRAFV600E, or Treated with Conditioned Media",
         color = colorRampPalette(c("white", "orange", "red"))(100))

#check total counts per sample
sapply(circ_lists, length)

#Check specific intersections:
intersect(circ_lists$Sample1, circ_lists$Sample2)

#Cluster the Jaccard matrix
install.packages("dendextend")
library(dendextend)

dist_mat <- as.dist(1 - jaccard_matrix)
hclust_res <- hclust(dist_mat)
plot(hclust_res, main = "Dendrogram of Shared circRNAs in Melanocytes Transduced 
     Senescent or Wild-Type, or Treated with Conditioned Media" )
# Assign colors to groups
groups <- c("Senescent", "Non-Senescent", "Senescent", "Senescent", "Conditioned_media", "Conditioned_media", "Non_Senescent", "Non_Senescent")
names(groups) <- colnames(jaccard_matrix)  # Same order as the dendrogram
dend <- as.dendrogram(hclust_res)
group_colors <- c("Mutant_BRAFVF00E1" = "royalblue", "Mutant" = "firebrick", "CM" = "forestgreen")
labels_colors(dend) <- group_colors[groups[order.dendrogram(dend)]]
------
groups <- c("Mutant_BRAFVF00E1", "Mutant_BRAFVF00E2", "Mutant_BRAFVF00E3", "Wild_BRAFVF00E1", "Wild_BRAFVF00E1", "Wild_BRAFVF00E1", "Conditioned_media1", "Conditioned_media2")
names(groups) <- colnames(jaccard_matrix)  # Same order as the dendrogram

# Distance and clustering
dist_mat <- as.dist(1 - jaccard_matrix)
hclust_res <- hclust(dist_mat)

# Convert to dendrogram object
dend <- as.dendrogram(hclust_res)

# Assign colors to groups
group_colors <- c("Mutant_BRAFVF00E1", "Mutant_BRAFVF00E2", "Mutant_BRAFVF00E3" = "royalblue", "Wild_BRAFVF00E1", "Wild_BRAFVF00E1", "Wild_BRAFVF00E1" = "firebrick", "Conditioned_media1", "Conditioned_media2" = "forestgreen")
labels_colors(dend) <- group_colors[groups[order.dendrogram(dend)]]

# Plot with colored labels
plot(dend, main = "Clustering of Samples by Jaccard Similarity")
legend("topright", legend = names(group_colors), fill = group_colors)


sample_metrics$Percent_Unique <- round(
  100 * sample_metrics$Unique_circRNAs / sample_metrics$Total_circRNAs,
  2
)

tool1_samples <- list(
  Mutant_BRAFV600E1 = tool1_sample1,
  Mutant_BRAFV600E2 = tool1_sample2,
  Mutant_BRAFV600E3 = tool1_sample3,
  Wild_BRAFV600E1 = tool1_sample4,
  Wild_BRAFV600E2 = tool1_sample5,
  Wild_BRAFV600E3 = tool1_sample6,
  Conditioned_media1 = tool1_sample7,
  Conditioned_media2 = tool1_sample8
)

print(tool1_samples)

sample_metrics <- data.frame(
  Sample = names(tool1_samples),
  Total_circRNAs = sapply(tool1_samples, function(df) length(unique(df$circ_id))),
  stringsAsFactors = FALSE
)

str(tool1_samples)




#Total and Unique circRNAs per sample
# Total circRNAs per sample
sample_metrics <- data.frame(
  Sample = names(tool1_samples),
  Total_circRNAs = sapply(tool1_samples, function(df) length(unique(df$circ_id))),
  stringsAsFactors = FALSE
)

# Calculate how many circRNAs are unique to each sample
all_ids <- unique(c(tool1_sample1$circ_id, tool1_sample2$circ_id, tool1_sample3$circ_id, tool1_sample4$circ_id, tool1_sample5$circ_id, tool1_sample6$circ_id, tool1_sample7$circ_id, tool1_sample8$circ_id))
length(all_ids)  # Total non-redundant circRNAs detected across tools
print(all_ids)

all_idss <- lapply(tool1_samples, function(df) unique(df$circ_id))
all_combined <- unlist(all_idss)
print (all_idss)
print(all_combined)


# Count frequency of each circRNA across all samples
circ_counts <- table(all_combined)

# Add unique circRNA counts per sample
sample_metrics$Unique_circRNAs <- sapply(all_idss, function(ids) {
  sum(circ_counts[ids] == 1)
})

#Consistency Score (mean Jaccard with other samples)
# Replace diagonal with NA to exclude self-comparison
jaccard_matrix[diag(nrow(jaccard_matrix)) == TRUE] <- NA

# Add consistency as the average pairwise Jaccard for each sample
sample_metrics$Mean_Jaccard <- rowMeans(jaccard_matrix, na.rm = TRUE)

print(sample_metrics)

#Barplot of Mean Jaccard Similarity 
library(ggplot2)

sample_metrics$Sample <- factor(sample_metrics$Sample,
                                levels = sample_metrics$Sample[order(sample_metrics$Mean_Jaccard)])

ggplot(sample_metrics, aes(x = Sample, y = Mean_Jaccard)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Jaccard Similarity per Sample (Sorted)",
       y = "Mean Jaccard Similarity", x = "Sample") +
  ylim(0, 1)

