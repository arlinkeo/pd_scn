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

# Check any overlap of DEGs and marker genes
ct_degs <- lapply(markerlist, function(l1){
  lapply(degs, function(l2){
    intersect(l1, rownames(l2$upregulated))
  })
})
ct <- sapply(ct_degs, function(m){
  !any(sapply(m, function(n){
    length(n)
  }) == 0)
})
ct_degs <- ct_degs[ct]

df <- summary_ttest[, "summary",unique(unlist(ct_degs)), c("Estimate", "BH")]
df <- asplit(df, 1)
df <- Reduce(cbind, df)
colnames(df) <- paste0(c("C.", "D."), colnames(df))
df[, c(1,3)] <- round(df[, c(1,3)], digits = 2)
df[, c(2,4)] <- format(df[, c(2,4)], digits = 3, scientific = TRUE)
ct <- unique(melt(ct_degs)[, -2])[,2]
df <- data.frame(Gene = entrezId2Name(rownames(df)), Marker = ct, df)
write.table(df, file = "output/celltype_degs.txt", row.names = FALSE, quote = FALSE, sep = "\t")  
  