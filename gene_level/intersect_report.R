library(zellkonverter) ## for read in h5 files
library(Seurat)
library(biomaRt)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringi)

## currently for mouse
# mart <- useEnsembl(
#     biomart = "ensembl",
#     dataset = "mmusculus_gene_ensembl",
#     mirror = "useast"
# )

make_upsetr_table <- function(sheet, name, species) {
    long <- sheet %>%
        pivot_longer(everything(), names_to = "comparison", values_to = "genes")
    
    long2plot <- na.omit(dcast(long, genes ~ comparison, fun.aggregate = length))
    dataset <- paste(species, "_gene_ensembl", sep = "", collapse = "")
    filters <- ifelse(species == "mmusculus", "mgi_symbol", "hgnc_symbol")
    mart <- useEnsembl(
        biomart = "ensembl",
        dataset = dataset,
        mirror = "useast"
    )
    annoteTable <- biomaRt::getBM(
        attributes = c(filters, "description"),
        filters = filters,
        values = long2plot$genes,
        mart = mart
    )
    ## report all genes including those without a description
    saveTable <- merge(annoteTable, long2plot, by.x = filters, by.y = "genes", all.y=T)
    ## numeric the columns except first two
    ncolumn <- ncol(saveTable)
    ## sort rows by sum of the numeric columns
    if (ncolumn>4){
    saveTable[, 3:ncolumn] <- apply(saveTable[, 3:ncolumn], 2, as.numeric)
    saveTable <- saveTable[order(-apply(saveTable[, 3:ncolumn], 1, sum)), ]
    }
    write.csv(saveTable, paste(c(name, "upsetr.csv"), collapse = "_"), row.names = F, quote = TRUE)
}
## name could include path
mk_upsetr_table_from_list <- function(gene_list, name, species) {
    gene_df <- stri_list2matrix(gene_list)
    make_upsetr_table(as.data.frame(gene_df), name, species)
}