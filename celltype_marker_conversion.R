# Cell-type marker conversion
library(homologene)

# Read markers
dir <- "../../Neuroexpresso_markers/markerGenesNCBI/All"
filenames <- list.files(dir)
conversion_table <- sapply(filenames, function(f){
  m <- unlist(read.table(paste0(dir, "/", f)))
  homologene(m, inTax = 10090, outTax = 9606) # Convert mouse entrez IDs to human ortholog entrez IDs
}, simplify = FALSE)

# Convert to human gene entrez IDs and filter for genes present in AHBA
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, ahba.genes())# Filter for genes present in AHBA
})

# Table with number of markers in mouse and human
df_size <- data.frame(Mouse = sapply(conversion_table, nrow), Human = sapply(markerlist, length))
df_size <- cbind('Cell-type' = rownames(df_size), df_size)
write.table(df_size, file = "output/number_of_celltype_markers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Check for duplicate markers 
all_markers <- unlist(markerlist)
dup <- all_markers[duplicated(all_markers)]
conversion_table <- lapply(names(conversion_table), function(x) {
  cbind(x, conversion_table[[x]])
})
conversion_table <- Reduce(rbind, conversion_table)
conversion_table[conversion_table$`9606_ID` %in% dup, ]