library(readxl)
pathway<- read_excel("pathway.xlsx",col_names = TRUE)
CPSS<- read_excel("CPSS.xlsx",col_names = TRUE)
library(dplyr)
genes <- CPSS$Gene
pathway2 <- lapply(pathway, function(x) {x[x %in% genes]})
pathway2
pathway2_csv <- pathway2 %>%stack() %>%
  rename(Gene = values,Pathway = ind) %>%
  select(Pathway, Gene)
write.csv(pathway2_csv,"pathway2.csv",row.names = FALSE)
library(tibble)
VST <- read_excel("VST.xlsx",col_names = TRUE)
expr <- VST %>%distinct(GeneSymbol,.keep_all = TRUE)
expr_matrix <- expr %>%column_to_rownames("GeneSymbol") %>%as.matrix()
######################################
library(GSVA)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
gene_sets <- pathway2
# pathway genes overlap with VST
gene_sets_filtered <- lapply(gene_sets,
  function(x){intersect(x, rownames(expr_matrix))})
# check overlap number genes
gene_number <- sapply(gene_sets_filtered,length)
gene_number
# only keep pathway at least 5 genes
gene_sets_filtered <- gene_sets_filtered[gene_number >= 5]
#### GSVA
gsva_param <- gsvaParam(exprData = expr_matrix,
                        geneSets = gene_sets_filtered,kcdf = "Gaussian")
gsva_scores <- gsva(gsva_param)
write.csv(gsva_scores,"Peroxisome_GSVA_scores.csv")
