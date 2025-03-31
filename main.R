# Project: Tumor-associated microbiota suppresses anti-tumor immunity by driving cDC1 dysfunction in lung adenocarcinoma
# Description: Analysis of single cell RNA sequencing of CD8 isolated from GF and SPF tumor bearing mice
# Author: QIANG DONG
# Last Updated: 2025-03-24

# Setup ------------------------------------------------------------------------
# Clear workspace and set working directory
rm(list=ls())
setwd(getwd())

# Create necessary directories
directories <- c("Figures", "Tables", "rdsData", "RData", "ProjecTILs")
for(dir in directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Function to download files from Dropbox
download_from_dropbox <- function(url, destfile) {
  # Convert to direct download URL by adding dl=1
  direct_url <- gsub("dl=0", "dl=1", url)
  
  cat("Downloading from:", direct_url, "\n")
  
  # Download the file based on OS
  if(Sys.info()["sysname"] == "Windows") {
    result <- try(download.file(url = direct_url, destfile = destfile, mode = "wb"))
  } else {
    result <- try(download.file(url = direct_url, destfile = destfile))
  }
  
  # Check if download was successful
  if (!inherits(result, "try-error") && file.exists(destfile) && file.size(destfile) > 1000) {
    cat("Successfully downloaded:", destfile, "\n")
    return(TRUE)
  } else {
    cat("Failed to download file. Please download manually from:", url, "\n")
    cat("And save to:", destfile, "\n")
    return(FALSE)
  }
}

# Download KP_sc_data.RData to RData folder
kp_data_url <- "https://www.dropbox.com/scl/fi/n89uzggv42zl8tpbhoujx/KP_sc_data.RData?rlkey=qq2m1ec2enqtza0z2asg9eni4&st=5hdx88r2&dl=0"
kp_data_file <- file.path("RData", "KP_sc_data.RData")
kp_download_success <- download_from_dropbox(kp_data_url, kp_data_file)

# Download ref_TILAtlas_mouse_v1.rds to ProjecTILs folder
ref_data_url <- "https://www.dropbox.com/scl/fi/3la672wfyzkj2fd3en31k/ref_TILAtlas_mouse_v1.rds?rlkey=p3da4v6he9bigfxb17hguvit7&st=pzob2xbg&dl=0"
ref_data_file <- file.path("ProjecTILs", "ref_TILAtlas_mouse_v1.rds")
ref_download_success <- download_from_dropbox(ref_data_url, ref_data_file)

# Create a user-writeable library path
user_lib <- file.path(getwd(), "R_libs")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)

# Set the library path
.libPaths(c(user_lib, .libPaths()))

# First, separate CRAN and Bioconductor packages
cran_packages <- c(
  "Seurat", # v5.2.1
  "harmony", # v1.2.3
  "tidyverse", # v2.0.0
  "patchwork", # v1.3.0
  "ggpubr", # v0.6.0
  "reshape2", # v1.4.4
  "ggplot2", # v3.5.1
  "RColorBrewer", # v1.1-3
  "dplyr", # v1.1.4
  "Matrix", # v1.6-4
  "cowplot", # v1.1.3
  "rstatix", # v0.7.2
  "writexl", # v1.5.1
  "conflicted" # Add the conflicted package
)

bioc_packages <- c(
  "limma", # v3.62.1
  "ComplexHeatmap", # v2.22.0
  "EnhancedVolcano" # v1.13.2
)

# Check and install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}

# Set non-interactive mode for BiocManager
options(BiocManager.check_repositories = FALSE)

# Install and load CRAN packages
for(pkg in cran_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = user_lib)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Install and load Bioconductor packages
for(pkg in bioc_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib = user_lib, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Install and load ProjecTILs specifically
if(!requireNamespace("ProjecTILs", quietly = TRUE)) {
  BiocManager::install("ProjecTILs", lib = user_lib, update = FALSE, ask = FALSE)
}

# Try loading ProjecTILs with error handling
tryCatch({
  suppressPackageStartupMessages(library(ProjecTILs))
}, error = function(e) {
  cat("Error loading ProjecTILs:", e$message, "\n")
  cat("You may need to manually install ProjecTILs or provide the path to its location.\n")
})

# Resolve common conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("setdiff", "base")
conflict_prefer("union", "base")
conflict_prefer("unname", "base")  # Add this one
conflict_prefer("expand", "tidyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")

# Load downloaded data --------------------------------------------------------
# Check if KP_sc_data.RData exists and try to load it
if(file.exists(kp_data_file) && file.size(kp_data_file) > 1000) {
  cat("Attempting to load:", kp_data_file, "\n")
  kp_data_loaded <- try(load(kp_data_file))
  if(inherits(kp_data_loaded, "try-error")) {
    cat("Error loading KP_sc_data.RData. The file may be corrupted.\n")
    cat("Please download the file again from:", kp_data_url, "\n")
  } else {
    cat("Successfully loaded KP_sc_data.RData\n")
  }
} else {
  cat("KP_sc_data.RData not found or too small.\n")
  cat("Please download manually from:", kp_data_url, "\n")
  cat("And save to:", kp_data_file, "\n")
}

# Check if ref_TILAtlas_mouse_v1.rds file exists
if(file.exists(ref_data_file) && file.size(ref_data_file) > 1000) {
  cat("ref_TILAtlas_mouse_v1.rds is ready to use in the ProjecTILs folder\n")
} else {
  cat("ref_TILAtlas_mouse_v1.rds not found or too small.\n")
  cat("Please download manually from:", ref_data_url, "\n")
  cat("And save to:", ref_data_file, "\n")
}

# ================================================================
# PART 1: Data Integration
# ================================================================

# Load RData containing original data
load("RData/KP_sc_data.RData")

# Define function to process lane data
processLane <- function(lane_data, lane_name) {
  lane_data$sample <- lane_name
  lane_data$cell_barcode <- colnames(lane_data)
  lane_data$orig.ident <- paste(lane_data$sample, lane_data$cell_barcode, sep="")
  lane_data <- RenameCells(lane_data, new.names = lane_data$orig.ident)
  return(lane_data)
}

# Process individual lanes
lane_1 <- processLane(lsLungdata$`Jin1-CC_1`, 'lane_1')
lane_2 <- processLane(lsLungdata$`Jin1-CC_2`, 'lane_2')

# Prepare hash table
hash$orig.ident <- paste(hash$sample, hash$cell_barcode, sep="")

# Merge Seurat objects
merged_seurat <- merge(lane_1, lane_2)

# Add hash information to merged seurat object
merged_seurat <- AddMetaData(merged_seurat, metadata = hash$num_features, col.name = "hash_num_features")
merged_seurat <- AddMetaData(merged_seurat, metadata = hash$feature_call, col.name = "hash_feature_call")
merged_seurat <- AddMetaData(merged_seurat, metadata = hash$num_umis, col.name = "hash_num_umis")

# Save merged data
saveRDS(merged_seurat, "rdsData/KP_merged.rds")

# ================================================================
# PART 2: Data Processing
# ================================================================

# Read in merged data
sc_data <- readRDS("rdsData/KP_merged.rds")

# Define experimental groups based on HTO information
sc_data$group <- dplyr::case_when(
  sc_data$hash_feature_call %in% c("HTO_A0301", "HTO_A0302", "HTO_A0303", "HTO_A0304") ~ "GF",
  sc_data$hash_feature_call %in% c("HTO_A0305", "HTO_A0307", "HTO_A0308") ~ "SPF",
  sc_data$hash_feature_call == "HTO_A0306" ~ "Healthy",
  TRUE ~ "Other"
)

# Add mitochondrial percentage
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "mt-")

# Filter cells based on QC metrics
sc_data_filtered <- subset(
  sc_data,
  subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
    nCount_RNA > 1000 & nCount_RNA < 10000 &
    percent.mt < 15
)

# Report filtering results
cat("Cells before filtering:", dim(sc_data)[2], "\n")
cat("Cells after filtering:", dim(sc_data_filtered)[2], "\n")

# Standard Seurat workflow
sc_data_filtered <- sc_data_filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Run Harmony for batch correction
sc_data_filtered <- sc_data_filtered %>%
  RunHarmony(group.by.vars = "sample", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 1.0)

# ================================================================
# PART 3: Cluster Annotation
# ================================================================

# Combine the layers
sc_data_filtered <- JoinLayers(sc_data_filtered)

# Load cell type annotations
celltypist_annotations <- read.csv("RData/celltypist_predicted_labels_Adult_Mouse_Gut.csv")
cell_type_vector <- setNames(celltypist_annotations[,2], celltypist_annotations[,1])
sc_data_filtered$celltypist_annotation <- cell_type_vector[colnames(sc_data_filtered)]

# Clean up metadata 
sc_data_filtered$seurat_clusters_original <- sc_data_filtered$seurat_clusters
columns_to_remove <- c("hash_num_features", "hash_num_umis", "percent.mt", 
                       "RNA_snn_res.0.8", "seurat_clusters", "cell_barcode")
columns_to_keep <- base::setdiff(colnames(sc_data_filtered@meta.data), columns_to_remove)
sc_data_filtered@meta.data <- sc_data_filtered@meta.data[, columns_to_keep]
colnames(sc_data_filtered@meta.data)[colnames(sc_data_filtered@meta.data) == "hash_feature_call"] <- "mouse"

# ================================================================
# PART 4: Analysis of T and NK Cells
# ================================================================

# Extract T and NK cells based on cluster IDs
T_NK_cells <- subset(sc_data_filtered, 
                     seurat_clusters_original %in% c(3, 4, 6, 7, 8, 10, 13, 15, 16, 17, 18))

# Function to process cell subsets
process_subset <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, k.param = 20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, min.dist = 0.5)
  return(seurat_obj)
}

# Process T and NK cells
T_NK_cells <- process_subset(T_NK_cells)

# Filter out specific clusters (cluster 14)
T_NK_cells_filtered <- subset(T_NK_cells, subset = seurat_clusters != 14)
T_NK_cells_filtered <- process_subset(T_NK_cells_filtered)

# ================================================================
# Figure_S1A: T_and_NK_cells_UMAP
# ================================================================

# Extract marker gene expression
CD3e_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd3e", ]
CD3d_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd3d", ]
CD3g_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd3g", ]
CD8a_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd8a", ]
CD8b1_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd8b1", ]
CD4_exp <- GetAssayData(T_NK_cells_filtered, slot = "data")["Cd4", ]

# Initialize cell type annotation column
T_NK_cells_filtered$T_NK_cells <- "Unidentified"  # Default value for all cells

# First pass of annotation based on clusters
T_NK_cells_filtered$T_NK_cells[T_NK_cells_filtered$seurat_clusters %in% c(4,5,7,12)] <- "CD4 T cells"
T_NK_cells_filtered$T_NK_cells[T_NK_cells_filtered$seurat_clusters %in% c(0,3,8,9)] <- "CD8 T cells"
T_NK_cells_filtered$T_NK_cells[T_NK_cells_filtered$seurat_clusters == 2] <- "NK cells"
T_NK_cells_filtered$T_NK_cells[T_NK_cells_filtered$seurat_clusters == 11] <- "NKT cells"
T_NK_cells_filtered$T_NK_cells[T_NK_cells_filtered$seurat_clusters == 10] <- "gamma delta T cells"

# For unidentified clusters (1,6,13), use marker expression to classify
unidentified_idx <- which(T_NK_cells_filtered$seurat_clusters %in% c(1,6,13))

# CD4 T cells: CD4 > 0.2 & CD8a = 0 & CD8b1 = 0
cd4_cells <- unidentified_idx[
  CD4_exp[unidentified_idx] > 0.2 & 
    CD8a_exp[unidentified_idx] == 0 & 
    CD8b1_exp[unidentified_idx] == 0
]
T_NK_cells_filtered$T_NK_cells[cd4_cells] <- "CD4 T cells"

# CD8 T cells: CD4 = 0 & (CD8a > 0.2 | CD8b1 > 0.2)
cd8_cells <- unidentified_idx[
  CD4_exp[unidentified_idx] == 0 & 
    (CD8a_exp[unidentified_idx] > 0.2 | CD8b1_exp[unidentified_idx] > 0.2)
]
T_NK_cells_filtered$T_NK_cells[cd8_cells] <- "CD8 T cells"

# Custom colors for cell types
custom_colors <- c(
  "CD8 T cells" = "#377EB8",
  "CD4 T cells" = "#FF7F00",
  "NK cells" = "#4DAF4A",
  "NKT cells" = "#984EA3",
  "gamma delta T cells" = "#DB7093",
  "Unidentified" = "grey"
)

# Reorder factor levels so "Unidentified" is plotted last
T_NK_cells_filtered$T_NK_cells <- factor(T_NK_cells_filtered$T_NK_cells, 
                                         levels = c("CD4 T cells", "CD8 T cells", "NK cells", 
                                                    "NKT cells", "gamma delta T cells", "Unidentified"))

# Create Figure S1A: UMAP of T/NK cells with reordered plotting
T_and_NK_cells_UMAP <- DimPlot(T_NK_cells_filtered, 
                               group.by = 'T_NK_cells', 
                               cols = custom_colors, 
                               order = levels(T_NK_cells_filtered$T_NK_cells)) + 
  ggtitle("T and NK cells (Filtered)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_fixed() + xlim(-13, 15) + ylim(-13, 15)

ggsave("Figures/Figure_S1A.pdf", 
       plot = T_and_NK_cells_UMAP, width = 6, height = 4.5)

# ================================================================
# Figure_1A: T_and_NK_cells_DensityPlot & EnrichmentPlot
# ================================================================

# Create Density Plot by group
data <- data.frame(
  UMAP_1 = T_NK_cells_filtered@reductions$umap@cell.embeddings[,1],
  UMAP_2 = T_NK_cells_filtered@reductions$umap@cell.embeddings[,2],
  z = T_NK_cells_filtered$group
)
data$z <- factor(data$z, levels = c("SPF", "GF", "Healthy"))

T_and_NK_cells_DensityPlot <- ggplot(data, aes(x=UMAP_1, y=UMAP_2)) + 
  coord_fixed() + xlim(-13, 15) + ylim(-13, 15) +
  geom_density_2d_filled(contour_var = "density") + 
  facet_wrap(~z, nrow=1, ncol=3) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold")
  ) +
  scale_fill_viridis_d(n=10, option = "viridis", alpha = 1)

ggsave("Figures/Figure_1A_Densityplot.pdf", 
       plot = T_and_NK_cells_DensityPlot, width = 12, height = 4)

# ProjecTILs Analysis for enrichment
ref_file <- "ProjecTILs/ref_TILAtlas_mouse_v1.rds"
ref_name <- "ref_TILAtlas_mouse_v1"

ref <- readRDS(ref_file)
projected_data <- Run.ProjecTILs(T_NK_cells_filtered, ref, filter.cells = TRUE)
T_NK_cells_filtered[[paste0("ProjecTILs_", ref_name)]] <- projected_data$query$functional.cluster

# Create plots for target clusters
target_clusters <- c("CD8_EffectorMemory", "CD8_NaiveLike")
plot_list_projectils <- list()

for(cluster in target_clusters) {
  cells_in_cluster <- rownames(T_NK_cells_filtered@meta.data)[T_NK_cells_filtered[[paste0("ProjecTILs_", ref_name)]] == cluster]
  
  p <- DimPlot(T_NK_cells_filtered, 
               reduction = "umap",
               cells.highlight = cells_in_cluster,
               cols.highlight = "red",
               cols = "grey",
               pt.size = 0.5) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    ggtitle(cluster) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    coord_fixed() + xlim(-13, 15) + ylim(-13, 15)
  
  plot_list_projectils[[cluster]] <- p
}

# Module Score Analysis
target_gene_file <- "RData/genelist/Mm_S_epi_Tc17_v_Tc1.csv"
gene_data <- read.csv(target_gene_file, header = FALSE, stringsAsFactors = FALSE)
list_name <- gene_data[1,1]
genes <- gene_data[-1,1]
available_genes <- genes[genes %in% rownames(T_NK_cells_filtered)]

# Add module score
T_NK_cells_filtered <- AddModuleScore(
  object = T_NK_cells_filtered,
  features = list(available_genes),
  name = list_name
)

# Create module score plot
module_score_plot <- FeaturePlot(
  T_NK_cells_filtered,
  features = paste0(list_name, "1"),
  reduction = "umap",
  min.cutoff = "q40",
  max.cutoff = "q90",
  cols = c("grey", "red"),
  order = TRUE
) + 
  ggtitle(paste0(list_name, "\n(", length(available_genes), "/", length(genes), " genes)")) + 
  coord_fixed() + xlim(-13, 15) + ylim(-13, 15) +
  theme(panel.grid = element_blank())

# Correlation Analysis function
analyze_specific_correlation <- function(seurat_obj, group_value, feature_type, feature_name) {
  # Get cells in the specified group
  group_cells <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data$group == group_value & 
      !is.na(seurat_obj@meta.data$group)
  ]
  
  # Get feature-positive cells
  if (feature_type == "ProjecTILs") {
    feature_column <- "ProjecTILs_ref_TILAtlas_mouse_v1"
    feature_cells <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[feature_column]] == feature_name & 
        !is.na(seurat_obj@meta.data[[feature_column]])
    ]
  } else {
    feature_column <- paste0(feature_name, "1")
    threshold <- quantile(seurat_obj@meta.data[[feature_column]], 0.75, na.rm = TRUE)
    feature_cells <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[feature_column]] > threshold & 
        !is.na(seurat_obj@meta.data[[feature_column]])
    ]
  }
  
  # Calculate overlap and statistics
  overlap_cells <- base::intersect(group_cells, feature_cells)
  total_cells <- ncol(seurat_obj)
  n_group <- length(group_cells)
  n_feature <- length(feature_cells)
  n_overlap <- length(overlap_cells)
  
  # Calculate expected overlap by chance
  expected_overlap <- (n_feature / total_cells) * n_group
  
  # Calculate p-value and enrichment ratio
  p_value <- phyper(n_overlap - 1, n_feature, total_cells - n_feature, n_group, lower.tail = FALSE)
  enrichment_ratio <- n_overlap / expected_overlap
  
  return(list(
    group = group_value,
    feature = feature_name,
    feature_type = feature_type,
    p_value = p_value,
    enrichment_ratio = enrichment_ratio,
    n_group = n_group,
    n_feature = n_feature,
    n_overlap = n_overlap,
    total_cells = total_cells
  ))
}

# Define comparisons
comparisons <- list(
  list(group = "SPF", feature_type = "Module", feature_name = "Mm_S_epi_Tc17_v_Tc1"),
  list(group = "GF", feature_type = "ProjecTILs", feature_name = "CD8_EffectorMemory"),
  list(group = "Healthy", feature_type = "ProjecTILs", feature_name = "CD8_NaiveLike")
)

# Run analyses
results <- list()
for (comparison in comparisons) {
  result <- analyze_specific_correlation(
    T_NK_cells_filtered, 
    comparison$group, 
    comparison$feature_type, 
    comparison$feature_name
  )
  
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Create summary table for visualization
summary_table <- data.frame(
  Group = sapply(results, function(x) x$group),
  Feature = sapply(results, function(x) x$feature),
  Feature_Type = sapply(results, function(x) x$feature_type),
  P_Value = sapply(results, function(x) x$p_value),
  Enrichment_Ratio = sapply(results, function(x) x$enrichment_ratio),
  Group_Cells = sapply(results, function(x) x$n_group),
  Feature_Cells = sapply(results, function(x) x$n_feature),
  Overlap_Cells = sapply(results, function(x) x$n_overlap),
  stringsAsFactors = FALSE
)

# Create label and format p-values for display
summary_table$label <- paste0(summary_table$Group, " + ", summary_table$Feature)
summary_table$p_display <- ifelse(
  summary_table$P_Value < 0.001, "p < 0.001",
  ifelse(summary_table$P_Value < 0.01, sprintf("p = %.3f", summary_table$P_Value),
         sprintf("p = %.1e", summary_table$P_Value))
)

# Export results to CSV
write.csv(summary_table, "Tables/Figure_1A_enrichment_analysis.csv", row.names = FALSE)

# Enrichment plot
enrichment_plot <- ggplot(summary_table, aes(x = reorder(label, Enrichment_Ratio), y = Enrichment_Ratio, fill = Feature_Type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = p_display), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3.5, fontface = "bold") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Feature Enrichment Analysis",
    x = "",
    y = "Enrichment Ratio",
    fill = "Feature Type"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Combine all plots
plot_list <- list()

# Add ProjecTILs plots
for (cluster in names(plot_list_projectils)) {
  plot_list[[paste0("ProjecTILs_", cluster)]] <- plot_list_projectils[[cluster]]
}

# Add Module Score plot
plot_list[["ModuleScore"]] <- module_score_plot

# Add enrichment plot
plot_list[["Enrichment"]] <- enrichment_plot

# Combine plots into a single figure
combined_plot <- wrap_plots(plot_list, ncol = 2)

# Add annotations to combined figure
combined_plot_with_annotation <- combined_plot + 
  plot_annotation(
    title = 'Figure 1A: T and NK Cell Analysis',
    subtitle = 'ProjecTILs classification, Module Score analysis, and Enrichment results',
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save individual enrichment plot
ggsave("Figures/Figure_1A_enrichment_plot.pdf", plot = enrichment_plot, width = 8, height = 6)

# Save combined figure
ggsave("Figures/Figure_1A_combined_plot.pdf", 
       plot = combined_plot_with_annotation, 
       width = 12, height = 10, 
       units = "in")

# ================================================================
# Figure_1C: CD8_T_cells_UMAP
# ================================================================

# Create logical vectors for selection
CD3_positive <- CD3e_exp > 0.2 | CD3d_exp > 0.2 | CD3g_exp > 0.2 
CD8_positive <- CD8a_exp > 0.2 | CD8b1_exp > 0.2
CD4_negative <- CD4_exp == 0

# Combine conditions for CD8+ T cells
CD8_cells <- CD8_positive & CD4_negative & CD3_positive
CD8_T_cells <- subset(T_NK_cells_filtered, cells = names(CD8_cells)[CD8_cells])

# Keep only SPF and GF groups
CD8_T_cells <- subset(CD8_T_cells, subset = group %in% c("SPF", "GF"))

# Process CD8 T cells
CD8_T_cells <- NormalizeData(CD8_T_cells)
CD8_T_cells <- FindVariableFeatures(CD8_T_cells, selection.method = "vst", nfeatures = 2000)
CD8_T_cells <- ScaleData(CD8_T_cells, features = rownames(CD8_T_cells))
CD8_T_cells <- RunPCA(CD8_T_cells, features = VariableFeatures(CD8_T_cells))
CD8_T_cells <- FindNeighbors(CD8_T_cells, dims = 1:30)
CD8_T_cells <- FindClusters(CD8_T_cells, resolution = 0.5)
CD8_T_cells <- RunUMAP(CD8_T_cells, 
                       dims = 1:30, 
                       n.neighbors = 60, 
                       min.dist = 0.5, 
                       spread = 0.5, 
                       metric = "cosine")

# Re-assign cluster IDs for better interpretation
new_numbers <- c("0" = "1", "1" = "3", "2" = "2", "3" = "4")
CD8_T_cells$seurat_clusters <- plyr::mapvalues(CD8_T_cells$seurat_clusters,
                                               from = names(new_numbers),
                                               to = new_numbers)
Idents(CD8_T_cells) <- "seurat_clusters"

# Create Figure 1C: CD8 T cell UMAP
CD8Umap <- DimPlot(CD8_T_cells, reduction = "umap", label = TRUE) + 
  coord_fixed() + xlim(-4.5, 5) + ylim(-4.5,5)
ggsave("Figures/Figure_1C.pdf", plot = CD8Umap, width = 6, height = 4)

# ================================================================
# Figure_1D: CD8_Cluster_distribution
# ================================================================

# Calculate percentages of each cluster in SPF vs GF
metadata <- CD8_T_cells@meta.data
cluster_proportions <- metadata %>%
  dplyr::filter(group %in% c("GF", "SPF")) %>%
  group_by(mouse, group, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(mouse, group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Function to create statistical comparison plots
create_cluster_plot <- function(data, cluster_num) {
  cluster_data <- data %>% dplyr::filter(seurat_clusters %in% cluster_num)
  stat_test <- cluster_data %>%
    t_test(percentage ~ group, var.equal = TRUE, paired = FALSE) %>%
    add_significance("p") %>%
    mutate(y.position = max(cluster_data$percentage) * 1.2,
           p = format(p, digits = 3))
  
  ggplot(cluster_data, aes(x = group, y = percentage, color = group)) +
    geom_point(size = 3.5, position = position_jitter(width = 0.2, seed = 123)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.4, linewidth = 0.7, color = "black") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.6, linewidth = 0.7, color = "black") +
    scale_color_manual(values = c("GF" = "#F71480", "SPF" = "black")) +
    stat_pvalue_manual(stat_test, label = "p = {p}", bracket.size = 0.8, tip.length = 0.02) +
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.line = element_line(linewidth = 0.8)) +
    labs(title = paste("Cluster", cluster_num),
         y = ifelse(cluster_num == 0, "Percentage of Cells (%)", ""),
         x = "") +
    scale_x_discrete(limits = c("SPF", "GF")) +
    scale_y_continuous(limits = c(0, max(cluster_data$percentage) * 1.3),
                       expand = expansion(mult = c(0, 0)))
}

# Create and save cluster comparison plots
clusters <- sort(unique(CD8_T_cells$seurat_clusters))
plot_list <- lapply(clusters, function(x) create_cluster_plot(cluster_proportions, x))
combined_plot <- patchwork::wrap_plots(plot_list, nrow = 1)
ggsave("Figures/Figure_1D.pdf", combined_plot, width = 9, height = 3.5, dpi = 300)

# ================================================================
# Figure_1E: CD8_Cluster_Feature_genes_Balloonplot
# ================================================================

# Create dotplot of key marker genes
markers <- c('Lef1','Sell', 'Cd44', 'Pdcd1', 'Tcf7', 'Tox','Ccr5', 'Gzmk', 'Cxcr5','Slamf6','Havcr2', 'Cxcr6', 'Id2', 'Lgals3')
colour <- ifelse(markers %in% c('Sell', 'Cd44', 'Pdcd1', 'Tcf7'), "red", "black")

# First, create a new factor with the exact order you want
new_cluster_levels <- c("1", "2", "3", "4")
CD8_T_cells <- SetIdent(CD8_T_cells, value = "seurat_clusters")
CD8_T_cells@active.ident <- factor(CD8_T_cells@active.ident, levels = new_cluster_levels)

# Create dotplot with custom settings
Balloon_plot <- DotPlot(CD8_T_cells, features = markers, scale = TRUE, scale.by = 'size') +
  scale_colour_gradientn(colors = rev(brewer.pal(30, 'RdYlBu')), limits = c(-2.6, 1.6)) +
  scale_size_continuous(range = c(1, 6)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x.bottom = element_text(color = colour))

# Remove the existing point layer
Balloon_plot$layers[[1]] <- NULL

# Add back custom layers - first filled dots, then outlines
Balloon_plot <- Balloon_plot + 
  geom_point(data = Balloon_plot$data, 
             aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
  geom_point(data = Balloon_plot$data, 
             aes(x = features.plot, y = id, size = pct.exp),
             shape = 1, color = "black", stroke = 0.5)

# Save the plot
ggsave("Figures/Figure_1E.pdf", plot = Balloon_plot, width = 7.5, height = 3.3)

# ================================================================
# Figure_S1D: CD8_Cluster_Feature_genes_Heatmap
# ================================================================

# Find markers for CD8 T cell clusters
CD8_T_cells$seurat_clusters <- factor(CD8_T_cells$seurat_clusters, levels = c("1", "2", "3", "4"))
CD8_markers <- FindAllMarkers(CD8_T_cells, 
                              only.pos = TRUE,  
                              min.pct = 0.25,   
                              logfc.threshold = 0.25)

# Process marker results
CD8_markers$cluster <- factor(CD8_markers$cluster, levels = c("1", "2", "3", "4"))
CD8_markers <- CD8_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::rename(
    Log2_Fold_Change = avg_log2FC,
    Percent_In_Cluster = pct.1,
    Percent_Out_Cluster = pct.2,
    Adjusted_P_value = p_val_adj
  )

# Save the complete marker details to CSV
write.csv(CD8_markers, "Tables/CD8_markers_complete.csv", row.names = FALSE)

# Function to create heatmaps
DoAllHeatmaps <- function(markers, dt, save_name='1', n=30, selected_genes=NULL) {
  # Ensure clusters are properly ordered
  cluster_order <- c("1", "2", "3", "4")
  dt$seurat_clusters <- factor(dt$seurat_clusters, levels = cluster_order)
  markers$cluster <- factor(markers$cluster, levels = cluster_order)
  
  # Get markers in the right order
  ordered_markers <- vector("character")
  for(clust in cluster_order) {
    cluster_markers <- markers %>%
      filter(cluster == clust) %>%
      arrange(desc(Log2_Fold_Change)) %>%
      slice_head(n = n)
    ordered_markers <- c(ordered_markers, cluster_markers$gene)
  }
  
  # Scale data
  dt <- ScaleData(dt, features = ordered_markers)
  
  # Create data frame for selected gene labels
  gene_positions <- data.frame(
    gene = ordered_markers,
    position = seq(length(ordered_markers), 1, length.out = length(ordered_markers))
  ) 
  
  if (!is.null(selected_genes)) {
    gene_positions <- gene_positions %>% filter(gene %in% selected_genes)
  }
  
  # Create plot
  plot <- DoHeatmap(dt, 
                    features = ordered_markers,
                    group.by = "seurat_clusters",
                    label = FALSE,
                    cells = WhichCells(dt, expression = seurat_clusters %in% cluster_order)) + 
    scale_fill_gradientn(colors = c("#4F0060", "black", "#FFF309")) +
    ggtitle(save_name) +
    geom_text(data = gene_positions,
              aes(x = Inf, y = position, label = gene),
              hjust = -0.1,
              size = 3) +
    theme(plot.margin = margin(r = 100, unit = "pt"))
  
  # Save the gene list used for the heatmap
  gene_list_df <- data.frame(
    Heatmap_Gene = ordered_markers,
    Cluster = rep(cluster_order, each = n)[1:length(ordered_markers)]
  )
  
  # Save the gene list to CSV
  write.csv(gene_list_df, paste0("Tables/", save_name, "_heatmap_gene_list.csv"), row.names = FALSE)
  
  # Also save just the selected genes if applicable
  if (!is.null(selected_genes)) {
    selected_gene_df <- data.frame(Selected_Genes = selected_genes)
    write.csv(selected_gene_df, paste0("Tables/", save_name, "_selected_genes.csv"), row.names = FALSE)
  }
  
  return(list(plot = plot, gene_list = gene_list_df))
}

# Create heatmap with key genes
genes_to_label <- c('Ccr7', 'Lef1', 'Ly6c2', 'Ifngas1', 'Il12rb2', 'Gzmk', 'Eomes', 'Ccr5', 'Maf', 'Fgl2', 'Nrp1', 'Itga1', 'Cxcr6', 'Lgals3')

# Create the heatmap and get the returned list
result <- DoAllHeatmaps(markers = CD8_markers, 
                        dt = CD8_T_cells, 
                        save_name = 'CD8', 
                        n = 30,
                        selected_genes = genes_to_label)

# Extract the plot from the result
Heatmap <- result$plot

# Save the heatmap
ggsave('Figures/Figure_S1D.pdf', Heatmap, width = 7, height = 12)


# ================================================================
# Figure_S1E: Marker_genes_check_in_Tex&Tpex
# ================================================================
# Balance sample sizes between SPF and GF groups
set.seed(42)
n_SPF <- sum(CD8_T_cells$group == "SPF")
GF_cells <- WhichCells(CD8_T_cells, expression = group == "GF")
GF_downsampled <- sample(GF_cells, size = n_SPF)
CD8_downsampled <- subset(CD8_T_cells, cells = c(GF_downsampled, WhichCells(CD8_T_cells, expression = group == "SPF")))

# Focus on Tpex (cluster 3) and Tex (cluster 4) clusters
umap_coords <- data.frame(
  UMAP_1 = CD8_downsampled[["umap"]]@cell.embeddings[,1],
  UMAP_2 = CD8_downsampled[["umap"]]@cell.embeddings[,2],
  cluster = CD8_downsampled$seurat_clusters,
  row.names = colnames(CD8_downsampled)
)

# Remove outliers using Mahalanobis distance
get_non_outliers <- function(cluster_num, threshold = 2) {
  coords <- umap_coords[umap_coords$cluster == cluster_num, c("UMAP_1", "UMAP_2")]
  center <- colMeans(coords)
  cov_matrix <- cov(coords)
  distances <- mahalanobis(coords, center, cov_matrix)
  return(rownames(coords)[distances < mean(distances) + threshold * sd(distances)])
}

cells_to_keep <- c(get_non_outliers("3"), get_non_outliers("4"))
CD8_filtered <- subset(CD8_downsampled, cells = cells_to_keep)
SPF_CD8_filtered <- subset(CD8_filtered, group == "SPF")
GF_CD8_filtered <- subset(CD8_filtered, group == "GF")

# Common plotting parameters
umap_limits <- list(
  x = c(-1.25, 3.75),
  y = c(-0.25, 4.75)
)

# Create a function that applies all your settings for the dimension plots
apply_dimplot_settings <- function(p, title) {
  p + 
    scale_color_manual(values = c("3" = "#F8756B", "4" = "#B59E00")) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6)
    ) +
    coord_fixed() +
    xlim(umap_limits$x) +
    ylim(umap_limits$y) +
    ggtitle(title)
}

# Create dimension plots
P1 <- apply_dimplot_settings(DimPlot(CD8_filtered, group.by = 'seurat_clusters'), 'Tpex vs Tex')
P2 <- apply_dimplot_settings(DimPlot(SPF_CD8_filtered, group.by = 'seurat_clusters'), 'SPF')
P3 <- apply_dimplot_settings(DimPlot(GF_CD8_filtered, group.by = 'seurat_clusters'), 'GF')

# Create feature plots for key marker genes
features <- c('Eomes', 'Ccr5', 'Gzmk', 'Cxcr6', 'Havcr2', 'Lgals3')

# Function to apply consistent settings to feature plots
apply_featureplot_settings <- function(p, gene_name) {
  p +
    scale_color_gradientn(
      colors = c("lightgrey", "#FFFFA0", "#FBF719", "#F85A3E", "#dc0000"),
      name = NULL # Remove legend title
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "none", # Remove legend for cleaner layout
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 6)
    ) +
    coord_fixed() +
    xlim(umap_limits$x) +
    ylim(umap_limits$y) +
    ggtitle(gene_name)
}

# Create feature plots
feature_plots <- list()
for(gene in features) {
  if(gene %in% rownames(CD8_filtered)) {
    p <- FeaturePlot(CD8_filtered, features = gene, order = TRUE, pt.size = 0.8)
    feature_plots[[gene]] <- apply_featureplot_settings(p, gene)
  }
}

# Arrange plots in a 3x3 grid using patchwork
library(patchwork)

Figure_S1E <- (
  P1 + P2 + P3 +
    feature_plots[["Eomes"]] + feature_plots[["Ccr5"]] + feature_plots[["Gzmk"]] +
    feature_plots[["Cxcr6"]] + feature_plots[["Havcr2"]] + feature_plots[["Lgals3"]]
) + 
  plot_layout(ncol = 3)

# Save the combined visualization with square aspect ratio
ggsave("Figures/Figure_S1E.pdf", plot = Figure_S1E, width = 8, height = 8, dpi = 600)

# ================================================================
# Save final processed objects
# ================================================================

saveRDS(T_NK_cells_filtered, "rdsData/T_NK_cells_final.rds")
saveRDS(CD8_T_cells, "rdsData/CD8_T_cells_final.rds")

cat("Analysis pipeline completed successfully!\n")