library(Seurat)
library(Matrix)
library(writexl)


library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(stringr)
library(purrr)

setwd("~/BINF/others/cilia/expgenes-enrprots")
mcc_markers13 <- read_excel("~/BINF/others/cilia/results/s13/mcc_markers13.xlsx")
exp_genes = mcc_markers13
exp_genes$rna_key = tolower(sub("^.*\\|", "", exp_genes$`Gene Name`))

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

joined_all_matches = inner_join(prot_tok, exp_genes, 
                                by = c("token" = "rna_key"))

match13 = joined_all_matches$token
match13 = as.data.frame(unique(match13))
write_xlsx(match13,"match13.xlsx")



match <- bind_rows(
  match27 %>% rename(gene = 1),
  match24 %>% rename(gene = 1),
  match22 %>% rename(gene = 1),
  match20 %>% rename(gene = 1),
  match18 %>% rename(gene = 1),
  match16 %>% rename(gene = 1),
  match13 %>% rename(gene = 1)
) %>%
  distinct(gene)

write_xlsx(match,"match_all.xlsx")
