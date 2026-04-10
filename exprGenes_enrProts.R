#comparing all genes in the filtered bc matrix with enriched proteins

library(Seurat)
library(Matrix)
library(writexl)


library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(stringr)
library(purrr)
# Read 10x .h5 file
mat0 = Read10X_h5("~/BINF/others/cilia/rna-paper/GSE158088_RAW_ALL/GSM4790538_scCapSt13_count.out/scCapSt13_count/outs/filtered_gene_bc_matrices_h5.h5")
if (is.list(mat0)) {
  mat0 <- mat0[["Gene Expression"]]
}
mcc27 <- read_excel("results/s27/mcc_clusters27.xlsx")
mcc27$CellBarcode <- paste0(mcc27$CellBarcode, "-1")

##
table2 <- table1 %>%
  filter(str_starts(CellBarcode, "st27_"),
         CellType == "Multiciliated") %>%
  mutate(new_barcode = paste0(str_remove(CellBarcode, "^st27_"), "-1"))
keep_cells <- intersect(colnames(mat0), table2$new_barcode)
mat <- mat0[, keep_cells]

library(stringr)



n_cells_expr <- Matrix::rowSums(mat > 0)

# Keep genes expressed in ≥1 cell
genes_expr <- n_cells_expr[n_cells_expr >= 1]

# Create output table
exp_genes <- data.frame(
  gene = names(genes_expr),
  n_cells = as.numeric(genes_expr)
)


# Save
write_xlsx(exp_genes, "genes_1st27xlsx")


##############
prot <- read_excel("~/BINF/others/cilia/prot.xlsx")
prot = prot %>%
  mutate(
    Genes = str_split(Genes, "\\s*;\\s*") |>
      map(~ str_replace(.x, "(?i)^(XELAEV_.*m)g$", "\\1")) |>
      map_chr(~ paste(.x, collapse = ";"))
  )

#expand proteome Genes into tokens
prot_tok = prot %>% mutate(token = str_split(Genes, "\\s*;\\s*")) %>%
  unnest(token) %>%
  mutate(token = str_to_lower(token))
prot_tok$token = tolower(prot_tok$token)


exp_genes$rna_key = tolower(sub("^.*\\|", "", exp_genes$gene))
joined_all_matches = inner_join(prot_tok, exp_genes, 
                                by = c("token" = "rna_key"))


setwd("~/BINF/others/cilia")
rnalist =  as.data.frame(exp_genes$rna_key)
protlist = as.data.frame(prot_tok$token)
colnames(rnalist) = "name"
colnames(protlist) = "name"
match = inner_join(protlist,rnalist)
match = unique(match)
no_match = anti_join(protlist,rnalist)
no_match = unique(no_match)
write_xlsx(match,"match27.xlsx")
write_xlsx(no_match,"no_match27.xlsx")
rnalist[grepl("sstr5\\.s", rnalist$name), ]

#####################
match27 <- read_excel("match27.xlsx")
match24 <- read_excel("match24.xlsx")

match2724 = rbind(match27, match24)
match2724 = unique(match2724)
write_xlsx(match2724,"match2724.xlsx")

match20 <- read_excel("match20.xlsx")
match272420 = rbind(match2724, match20)
match272420 = unique(match272420) 
aa = anti_join(match20, match2724)
match22 <- read_excel("match22.xlsx")

match2724 <- bind_rows(match27, match24) %>%
  distinct(name, .keep_all = TRUE)

match272422 <- bind_rows(match2724, match22) %>%
  distinct(name, .keep_all = TRUE)


match27242220 <- bind_rows(match272422, match20) %>%
  distinct(name, .keep_all = TRUE)

match2724222018 <- bind_rows(match27242220, match18) %>%
  distinct(name, .keep_all = TRUE)

match272422201816 <- bind_rows(match2724222018, match16) %>%
  distinct(name, .keep_all = TRUE)

match_all <- bind_rows(match272422201816, match13) %>%
  distinct(name, .keep_all = TRUE)

write_xlsx(match_all,"rnaexp_protenr.xlsx")

no_match_any = anti_join(protlist, match_all)
no_match_any = unique(no_match_any)
write_xlsx(no_match_any,"no_match_any.xlsx")
match_all = inner_join(match_all, protlist)
match_all = unique(match_all)

##############
matches = joined_all_matches$token

matches = unique(matches)
print(matches)

write_xlsx(prot_tok,"prot_tok.xlsx")
write_xlsx(exp_genes,"exp_genes.xlsx")


protlist = as.data.frame(prot_tok$token)
colnames(protlist) = "name"