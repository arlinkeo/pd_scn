# Cell-type enrichment of differentially expressed genes
library(reshape2)

# Hypergeometric test
hyper.test <- function(a, b, total){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

# Test for two lists of gene sets and correct P for cell-types tested
hyper.test.table <- function(l1, l2){ # two lists of gene sets
  pvalue <- sapply(names(l1), function(n){
    print(paste0(which(names(l1)==n), ": ", n))
    set <- l1[[n]]
    # Overlap with each module
    sapply(l2, function(mod_genes){
      hyper.test(mod_genes, set, length(ahba.genes()))
    })
  })
  apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
}

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