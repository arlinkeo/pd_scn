# Cell-type markers
setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)
library(homologene)
library(ComplexHeatmap)
library(RColorBrewer)

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

# Function for gene conversion
probeInfo <- read.csv("../AHBA_Arlin/probe_info_2018-11-18.csv")
entrezId2Name <- function (x) {probeInfo$gene_symbol[match(x, probeInfo$entrez_id)]} #Input is vector
name2EntrezId <- function (x) {as.character(probeInfo$entrez_id[match(x, probeInfo$gene_symbol)])} #Input is vector

# brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))
sample_info <- readRDS("resources/sample_info.rds")

# Read markers
dir <- "C:/Users/dkeo/surfdrive/Neuroexpresso_markers/markerGenesNCBI/All"
filenames <- list.files(dir)
conversion_table <- sapply(filenames, function(f){
  m <- unlist(read.table(paste0(dir, "/", f)))
  # Convert mouse entrez IDs to human ortholog entrez IDs
  homologene(m, inTax = 10090, outTax = 9606)
}, simplify = FALSE)
sapply(conversion_table, nrow)
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, rownames(brainExpr$donor9861))# Filter for genes presentin AHBA
}) 

# duplicates
all_markers <- unlist(markerlist)
dup <- which(duplicated(all_markers))
all_markers[dup]

# Sample info per donor
ontology <- read.csv("../AHBA_Arlin/Ontology.csv")
sample_annot <- lapply(donorNames, function(d){
  read.csv(paste0("../AHBA_Arlin/sample_info_", d, "_2018-11-18.csv"))
})

#Heatmaps

celltype <- rep(names(markerlist), sapply(markerlist, length))#[-dup]
row_annot <- data.frame(celltype = celltype)
row_color <- rainbow(length(markerlist))
names(row_color) <- names(markerlist)
row_color <- row_color[celltype]

ahba_color <- ontology$color_hex_triplet
ahba_color <- paste0("#", ahba_color)
names(ahba_color) <- ontology$id
col_color <- list(ahba_color = ahba_color)

expr <- 
  lapply(donorNames, function(d){
    e <- brainExprNorm[[d]]
    samples <- sample_info[[d]]
    networks <- apply(samples[, -c(1:3)], 2, function(i){
      t <- as.matrix(e[all_markers, i==1])
      
      col_annot <- ontology[match(colnames(t), ontology$id), c("id")]
      
      Heatmap(celltype, name = "cell-type", width = unit(5, "mm")) +
      Heatmap(t, column_title = "Samples", row_title = "Markers",
              row_names_gp = gpar(fontsize = 3), split = celltype,
              column_names_side = c("top")#,
              # top_annotation = HeatmapAnnotation(col_annot, col = col_color)
              )
      
    })
     
  })

  