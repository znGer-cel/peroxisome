library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
###
go_terms <- c(
  #Peroxisome biogenesis & protein import
  "GO:0016558" = "protein_import_into_peroxisome_matrix",
  "GO:0007031" = "peroxisome_organization",
  #ROS detoxification
  "GO:0042744" = "hydrogen_peroxide_catabolic_process",
  "GO:0098869" = "cellular_oxidant_detoxification",
  #fatty acid β-oxidation
  "GO:0006635" = "fatty_acid_beta_oxidation",
  #fatty acidα-oxidation
  "GO:0001561" = "fatty_acid_alpha_oxidation",
  #acyl-CoA metabolism
  "GO:0035337" = "fatty_acyl_CoA_metabolic_process",
  #ether lipid metabolism
  "GO:0046485" = "ether_lipid_metabolic_process",
  #bile acid metabolism
  "GO:0008206" = "bile_acid_metabolic_process",
  #fatty acid transport 
  "GO:0015908" = "fatty_acid_transport",
  #lipid transport
  "GO:0006869" = "lipid_transport",
  #isoprenoid metabolism
  "GO:0006720" = "isoprenoid_metabolic_process",
  #amino acid & polyamine metabolism
  "GO:0009063" = "amino_acid_catabolic_process",
  "GO:0006595" = "polyamine_metabolic_process")
go_genes <- AnnotationDbi::select(org.Mm.eg.db,
  keys = names(go_terms),keytype = "GOALL",columns = "SYMBOL")
go_genes <- go_genes %>%filter(!is.na(SYMBOL)) %>%distinct(GOALL, SYMBOL)
go_genes$Pathway <- go_terms[go_genes$GOALL]
gene_table <- go_genes %>%group_by(Pathway) %>%mutate(row = row_number()) %>%
  dplyr::select(-GOALL) %>%
  pivot_wider(names_from = Pathway,values_from = SYMBOL) %>%
  dplyr::select(-row)
gene_table
write.csv(gene_table,"CPSS_Function_GO_gene.csv",row.names = FALSE)
##GO:0051216 cartilage development
library(GO.db)
go_id <- "GO:0051216"
offspring <- GOBPOFFSPRING[[go_id]]
go_all <- unique(c(go_id, offspring))
length(go_all)
go_all
genes <- AnnotationDbi::select(org.Mm.eg.db,keys = go_all,keytype = "GO",columns = c("GO", "SYMBOL"))
#Remove duplicates and NA values.
genes <- genes[!is.na(genes$SYMBOL), ]
genes <- unique(genes)
gene_symbols <- sort(unique(genes$SYMBOL))
length(gene_symbols)
go_name <- Term(GOTERM[[go_id]])
col_name <- gsub(" ", "_", go_name)
gene_list <- data.frame(gene_symbols,stringsAsFactors = FALSE)
colnames(gene_list) <- col_name
write.csv(gene_list,"GO0051216_with_offspring_Mouse_Genes.csv",row.names = FALSE,quote = FALSE)
