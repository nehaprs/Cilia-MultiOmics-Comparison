# ===== Setup =====
install.packages("clustermole")  
BiocManager::install("GSVA")
BiocManager::install("singscore")
install.packages("XML")
library(dplyr)
library(tidyr)
library(stringr)
library(clustermole)
library(singscore)
library(readxl)
library(writexl)
# ===== Inputs =====
setwd("~/BINF/yushi scrnaseq/all six/threshold0/split/pax3")
markers_df <- read_excel("filtmarkers_resolution_1.5.xlsx")

# Parameters
species <- "mm"          # mouse
top_n_per_cluster <- 50  # genes per cluster to use
min_fc <- 0              # filter weak markers if desired

# Harmonize column names
if (!"avg_log2FC" %in% names(markers_df) && "avg_logFC" %in% names(markers_df)) {
  markers_df <- markers_df %>% rename(avg_log2FC = avg_logFC)
}

# ===== Prepare marker lists per cluster =====
marker_lists <-
  markers_df %>%
  filter(avg_log2FC > min_fc) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n_per_cluster) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

# ===== Annotate each cluster via CellMarker + PanglaoDB overlap =====
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

annot_one <- function(genes_input, species = "mm", top_per_source = 3) {
  # flatten any list/factor to unique character
  genes_chr <- unique(stats::na.omit(as.character(unlist(genes_input, use.names = FALSE))))
  if (length(genes_chr) == 0L) return(tibble())
  
  res <- clustermole::clustermole_overlaps(genes = genes_chr, species = species)
  if (is.null(res) || nrow(res) == 0L) return(tibble())
  res <- as_tibble(res)
  
  # standardize names
  names(res) <- tolower(names(res))
  if (!"source" %in% names(res)) {
    res$source <- dplyr::coalesce(res$db, res$database, "unknown")
  }
  if (!"cell_type" %in% names(res)) {
    res$cell_type <- dplyr::coalesce(res$celltype, res$label)
  }
  if (!"overlap" %in% names(res)) {
    res$overlap <- dplyr::coalesce(res$k, res$hits, NA_real_)
  }
  if (!"p_value" %in% names(res)) {
    res$p_value <- dplyr::coalesce(res$`p.val`, res$pval, res$p, NA_real_)
  }
  # reference size if provided under various names
  if (!"n_ref" %in% names(res)) {
    res$n_ref <- dplyr::coalesce(res$ref_size, res$reference_size, res$n_marker, res$markers, NA_real_)
  }
  
  # derive precision/recall if missing
  if (!"precision" %in% names(res)) {
    res$precision <- as.numeric(res$overlap) / length(genes_chr)
  }
  if (!"recall" %in% names(res)) {
    res$recall <- ifelse(is.na(res$n_ref), NA_real_, pmin(1, as.numeric(res$overlap) / as.numeric(res$n_ref)))
  }
  
  res %>%
    mutate(
      source = toupper(as.character(source)),
      score  = as.numeric(overlap)
    ) %>%
    arrange(desc(score), desc(precision), p_value) %>%
    select(cell_type, score, precision, recall, p_value, source) %>%
    group_by(source) %>%
    slice_head(n = top_per_source) %>%
    ungroup()
}

# build annotations
annot_tbl <- marker_lists %>%
  mutate(hits = purrr::map(genes, ~annot_one(.x, species)))

# pick best label per cluster
best_labels <- annot_tbl %>%
  tidyr::unnest(hits) %>%
  group_by(cluster) %>%
  arrange(desc(score), desc(precision), p_value) %>%
  slice_head(n = 1) %>%
  ungroup()


# ===== Export helpful artifacts =====
# 1) Table of best labels
writexl::write_xlsx(best_labels, "cluster_labels_clustermole.xlsx")

# 2) Full ranked hits per cluster
full_hits <-
  annot_tbl %>%
  unnest(hits) %>%
  arrange(cluster, desc(score), desc(precision), p_value)

write_xlsx(full_hits, "cluster_label_candidates_full.xlsx")

