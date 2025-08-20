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
sum(rowSums(is.na(st24mcc)))



st24 = meta[meta$Stages == "st24",]
unique(st24$CellType)
st24mcc = st24[st24$CellType == "Multiciliated",]
#remove NA rows
st24mcc = st24mcc[rowSums(is.na(st24mcc)) != ncol(st24mcc), ]

#rename barcode by removing prefix st24_


st24mcc = mutate(st24mcc, CellBarcode = sub("^st24_", "", CellBarcode))

#load custer matric 
library(readr)
clusters <- read_csv("~/BINF/others/cilia/rna-paper/GSE158088_RAW/GSM4790543_scCapSt24_count.out/scCapSt24_count/outs/analysis/clustering/graphclust/clusters.csv")
#remove -1 in barcode name, which is added by cellranger
clusters = mutate(clusters, Barcode = sub("-1$","", Barcode))

mcc_clusters = inner_join(st24mcc, clusters, by = c("CellBarcode" = "Barcode"))
write_xlsx(mcc_clusters, "mcc_clusters24.xlsx")


#diffexp
markers <- read_excel("mcc_markers.xlsx")
markers = markers[markers$`Cluster 9 Adjusted p value` < 0.05,]
markersup = markers[markers$`Cluster 9 Log2 fold change` > 0,]
write_xlsx(markersup,"mccUpGenes_st24.xlsx")

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

write_xlsx(joined_per_prot, "st24_overlap.xlsx")
