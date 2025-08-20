######################
#gather useful  information from supplementary table 1: cell metadata
######################
library(readxl)
library(dplyr)
library(writexl)

setwd("~/BINF/others/cilia/results")

#load metadata
meta <- read_excel("BINF/others/cilia/rna-paper/supplementary tables/table1.xlsx")


st24 = meta[meta$Stages == "st24",]
unique(st24$CellType)

st24mcc = st24[st24$CellType == "Multiciliated",]
#remove NA rows
st24mcc = st24mcc[rowSums(is.na(st24mcc)) != ncol(st24mcc), ]
#282 mcc in st24, of which 261 aren't NA total 3107 in st24
sum(rowSums(is.na(meta))) #4186 total NA rows

st27 = meta[meta$Stages == "st27",]
st27mcc = st27[st27$CellType == "Multiciliated",]
st27mcc = st27mcc[rowSums(is.na(st27mcc)) != ncol(st27mcc), ]
#st27: 3371 total, 293 mcc, 264 non-NA mccs

##Per Stage Louvain represented and their Phenograph

perstage_louvain_vs_phenograph <- meta %>%
  filter(Stages == "st24", CellType == "Multiciliated") %>%
  count(`Per Stage Louvain`, Phenograph, name = "n") %>%
  arrange(desc(n))

louv24 = st24 %>% count(`Per Stage Louvain`, CellType, Phenograph, name = "n") %>%
  arrange(desc(n))

write_xlsx(louv24, "lovian24.xlsx")



louv27 = st27 %>% count(`Per Stage Louvain`, CellType, Phenograph, name = "n") %>%
  arrange(desc(n))

write_xlsx(louv27, "lovian27.xlsx")

