#################
#find the relevant clusters and marker genes
##################

library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(stringr)
library(purrr)

setwd("~/BINF/others/cilia/results")

#load metadata
meta <- read_excel("BINF/others/cilia/rna-paper/supplementary tables/table1.xlsx")
#sum(rowSums(is.na(st24mcc)))

stage = "st22"

st = meta[meta$Stages == stage,]
unique(st$CellType)
stmcc = st[st$CellType == "Multiciliated",]
#remove NA rows
stmcc = stmcc[rowSums(is.na(stmcc)) != ncol(stmcc), ]

#rename barcode by removing prefix st24_

##choose########
stmcc = mutate(stmcc, CellBarcode = sub("^st22_", "", CellBarcode))

#load cluster matric 
library(readr)
clusters <- read_csv("~/BINF/others/cilia/rna-paper/GSE158088_RAW_ALL/GSM4790542_scCapSt22_count.v3.out/scCapSt22_count.v3/outs/analysis/clustering/graphclust/clusters.csv")
#remove -1 in barcode name, which is added by cellranger
clusters = mutate(clusters, Barcode = sub("-1$","", Barcode))

mcc_clusters = inner_join(stmcc, clusters, by = c("CellBarcode" = "Barcode"))
write_xlsx(mcc_clusters, "mcc_clusters22.xlsx")


#diffexp
markers <- read_excel("mcc_markers.xlsx")
markers = markers[markers$`Cluster 6 Adjusted p value` < 0.05,]
markersup = markers[markers$`Cluster 6 Log2 fold change` > 0,]
write_xlsx(markersup,"mccUpGenes_st22.xlsx")

#############
#compare proteome and rnaseq data
rnaseq = markersup

prot <- read_excel("~/BINF/others/cilia/prot.xlsx")

#convert XELEAV....mg into XELEAV...m
prot = prot %>%
  mutate(
    Genes = str_split(Genes, "\\s*;\\s*") |>
      map(~ str_replace(.x, "(?i)^(XELAEV_.*m)g$", "\\1")) |>
      map_chr(~ paste(.x, collapse = ";"))
  )


#separate gene names in rnaseq data
rnaseq$Genes.rna = sub("^.*\\|\\s*", "", markersup$`Gene Name`)

#clean proteome data and join

#normalize rnaseq key
rna_key = rnaseq %>% mutate(geneKey = str_to_lower(str_trim(Genes.rna)))

#expand proteome Genes into tokens
prot_tok = prot %>% mutate(token = str_split(Genes, "\\s*;\\s*")) %>%
  unnest(token) %>%
  mutate(token = str_to_lower(token))



#inner join
joined_all_matches = inner_join(prot_tok, rna_key, 
                                by = c("token" = "geneKey"))

#collapse back to one row per proteome record
joined_per_prot = joined_all_matches %>%
  group_by(across(names(prot))) %>%
  summarise(
    matched_tokens = paste(unique(token), collapse = ";"),
    matched_Gene.rnaseq = paste(unique(Genes.rna), collapse = ";"),
    .groups = "drop"
    
  )

write_xlsx(joined_per_prot, "st22_overlap.xlsx")

#for overlap in any stage: copy-pasted, now find uniquw
overlap_in_any_stage <- read_excel("~/BINF/others/cilia/results/overlap_in_any_stage.xlsx")
overlap_in_any_stage <- unique(overlap_in_any_stage)
write_xlsx(overlap_in_any_stage,"overlap_in_any_stage_unique.xlsx")
