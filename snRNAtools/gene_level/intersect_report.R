library(zellkonverter) ## for read in h5 files
library(Seurat)
library(biomaRt)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringi)


make_upsetr_table <- function(sheet, name, species) {
    long <- sheet %>%
        pivot_longer(everything(), names_to = "comparison", values_to = "genes")
    ## deal with zero length
    long2plot <- na.omit(dcast(long, genes ~ comparison, fun.aggregate = length)) ## this is what I need
    mart <- useEnsembl(
            biomart = "ensembl",
            ## mmusculus
            ## hsapiens
            dataset = paste(c(species, "_gene_ensembl"), sep = "", collapse = ""),
            mirror = "useast"
    )
    gene_symbol <- ifelse(species == "hsapiens", "hgnc_symbol", "mgi_symbol")
    annoteTable <- biomaRt::getBM(
        attributes = c(gene_symbol, "description"),
        filters = gene_symbol,
        values = long2plot$genes,
        mart = mart
    )
    saveTable <- merge(annoteTable, long2plot, by.x = gene_symbol, by.y = "genes")
    ## numeric the columns except first two
    saveTable[, 3:ncol(saveTable)] <- apply(saveTable[, 3:ncol(saveTable)], 2, as.numeric)
    ## sort rows by sum of the numeric columns
    saveTable <- saveTable[order(-apply(saveTable[, 3:ncol(saveTable)], 1, sum)), ]
    write.csv(saveTable, paste(c(name, "upsetr.csv"), collapse = "_"), row.names = F, quote = TRUE)
}
## name could include path
mk_upsetr_table_from_list <- function(gene_list, name, species) {
    gene_df <- stri_list2matrix(gene_list)
    ## currently for mouse
    make_upsetr_table(as.data.frame(gene_df), name, species)
}