# Cell-type marker conversion
library(homologene)

# Read markers
dir <- "../../Neuroexpresso_markers/markerGenesNCBI/All"
filenames <- list.files(dir)
conversion_table <- sapply(filenames, function(f){
  m <- unlist(read.table(paste0(dir, "/", f)))
  t <- homologene(m, inTax = 10090, outTax = 9606) # Convert mouse entrez IDs to human ortholog entrez IDs
}, simplify = FALSE)
names(conversion_table)[names(conversion_table) %in% c("Microglia_activation", "Microglia_deactivation")] <- c("Activated microglia", "Deactivated microglia")

# Convert to human gene entrez IDs and filter for genes present in AHBA
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, ahba.genes())# Filter for genes present in AHBA
})

# Check for duplicate markers 
all_markers <- unlist(markerlist)
dup <- all_markers[duplicated(all_markers)]
conversion_table <- sapply(names(conversion_table), function(x) {
  t <- conversion_table[[x]]
  t <- t[!(t$`9606_ID` %in% dup), ]
  cbind(x, t)
}, simplify = FALSE)
markerlist <- lapply(conversion_table, function(x) {
  l <- x[, "9606_ID"] # Human entrez IDs
  intersect(l, ahba.genes())# Filter for genes present in AHBA
})
saveRDS(markerlist, file = "output/markerlist.rds")

# Table with number of markers in mouse and human
df_size <- data.frame(Mouse = sapply(conversion_table, nrow), Human = sapply(markerlist, length))
df_size <- cbind('Cell-type' = rownames(df_size), df_size)
df_size <- do.call("cbind", split(df_size, c(rep(1, nrow(df_size)/2), rep(2, nrow(df_size)/2))))
write.table(df_size, file = "output/number_of_celltype_markers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

