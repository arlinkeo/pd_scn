# Cell-type enrichment of differentially expressed genes
library(reshape2)
library(plyr)

# Filter cell-types with at least 6 markers
markerlist6 <- markerlist[sapply(markerlist, length) >= 6]
cat(paste("Cell-types with at least 6 genes:\n", 
            tolower(gsub("_", " ", paste(names(markerlist6), collapse = ", ")))))

# Run cell-type enrichment with hyper.test
ct_enrichment <- lapply(degs, function(n){
  de_genes <- lapply(n, rownames)
  hyper.test.table(de_genes, markerlist6)
})
ct_enrichment <- simplify2array(ct_enrichment)
ct_enrichment <- alply(ct_enrichment, c(2,3), function(x){
  format(x[x < 0.05], digits = 3, scientific = TRUE)
})
names(ct_enrichment) <- apply(attributes(ct_enrichment)$split_labels, 1, paste, collapse = "_")
ct_enrichment
lapply(ct_enrichment, function(x){
  paste(paste0(names(x), " cells (P", " = ", x, ")"), collapse = ", ")
})
