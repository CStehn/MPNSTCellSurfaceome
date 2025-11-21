# Helper functions for Stehn et al. 2026
# Christopher Stehn
# November 21, 2025

library(tidyverse)
library(httr)
library(jsonlite)
library(xml2)
library(XML)

# Below function is to take files containing expression data (MS1 iBAQ) that are from ProteomicsDB in a directory for a specific protein
# Need to change it, as currently it is only able to pull tables listed by their gene name, not the protein isoform
# This requires a manual renaming of each file after manually downloading from PDB's website
parsePDBfiles <- function(genes, dir) {
  # get all files containing gene names in genes variable
  # pull them all together into table

  # just in case there is a trailing '/' in dir, remove it quickly
  # make glob of gene names, then get files containing
  dir <- gsub("/$", "", dir)
  geneglob <- paste0(genes, collapse = "|")
  files <- list.files(dir, pattern = geneglob, full.names = TRUE)

  # use lapply to read all tables into their own data frame
  # use select over this to ensure we only get Tissue and Average Intensity (Normalized)
  tables <- lapply(files, function(x) dplyr::select(read.table(x, sep = ";", header = TRUE), "Tissue", "Average.Normalized.Intensity"))
  # use reduce to do a full join and keep all possible tissues but get values for every protein in a data frame
  pdbtable <- purrr::reduce(tables, dplyr::full_join, by = "Tissue")
  # before assigning column names to our genes, we need to make sure we actually have all the genes
  # some genes might not have been available from PDB; if this is the case, our column names will not work
  # also, list.files() lists files in alphanumeric order and our genes are not, so arrange them to be so
  filegenes <- gsub(".*\\s(.*).csv", "\\1", files)
  genesremaining <- genes[match(filegenes, genes)]

  colnames(pdbtable)[-1] <- genesremaining
  return(pdbtable)
}


# NOTE: FOR A GET REQUEST, OFTEN WILL NEED TO USE '' CHARACTERS, SO ENCASE URL IN ""
# in request, the resulting object will be an http response, encoding errors, site info, response times, and site content
# site content is usually a file, either html, or with APIs: json or xml
# displaying object on console will let you know what type of file it is
# Proteomics DB files are available in xml or json
# strongly recommend the json files
getPDBexpression <- function(protein) {
  # get json file for a protein and its expression in various human tissues from ProteomicsDB
  # will narrow the file to only get normalized_expression, tissue_source, protein ID
  url <- "proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams(PROTEINFILTER="
  protformat <- paste0("'", protein, "'")
  protSymbol <- protein
  try(protSymbol <- suppressMessages(clusterProfiler::bitr(protein, "UNIPROT", "SYMBOL", org.Hs.eg.db::org.Hs.eg.db)$SYMBOL), silent = TRUE)

  tailend <- ",MS_LEVEL=1,TISSUE_ID_SELECTION='',TISSUE_CATEGORY_SELECTION='tissue;fluid',SCOPE_SELECTION=1,CALCULATION_METHOD=0,GROUP_BY_TISSUE=1,EXP_ID=-1)/Results?$select=UNIQUE_IDENTIFIER,TISSUE_NAME,NORMALIZED_INTENSITY&$format=json"
  httpreq <- GET(paste0(url, protformat, tailend))
  # format of content is in raw characters, need to convert
  # format of resulting json file is a table, results, with lines stored under d always
  jsontbl <- jsonlite::fromJSON(rawToChar(httpreq$content))$d$results
  # use tryCatch to return final data frame if length has elements, return error message if not
  tryCatch(
    {
      exp <- dplyr::select(jsontbl, TISSUE_NAME, NORMALIZED_INTENSITY)
      message(paste0("Table returned for ", protSymbol, "."))
      return(exp)
    },
    error = function(e) {
      message(paste0("Request for ", protSymbol, " table unsuccessful."))
      return(NULL)
    }
  )
}

plotHumanRNA <- function(longRNA, prot) {
  # simple wrapper to plot a protein from the filtered protein list
  # requires RNA data frame in long format with specified disease states
  humanRNAlong %>%
    filter(SYMBOL == prot) %>%
    ggplot(aes(x = sampleType, y = TPM, fill = sampleType)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    guides(fill = "none") +
    ylab("Expression (TPM)") +
    scale_x_discrete(labels = c("PN", "pNF", "ANNUBP", "MPNST")) +
    ggtitle(paste0(prot))
}
