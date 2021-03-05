# GO enrichment

library("clusterProfiler")
library("org.Hs.eg.db")

bg <- ahba.genes()
go_enrichment <- lapply(ll_degs, function(g){
  enrichGO(gene = g, OrgDb = "org.Hs.eg.db", ont = "ALL", universe = bg)
})
  
go_tables <- lapply(go_enrichment, function(x){
  x@result
})  
sapply(go_tables, nrow)  

# Write resulting tables
lapply(names(go_tables), function(n){
  t <- go_tables[[n]]
  write.table(t, file = paste0("output/go_enrichment_", gsub("\\.", "_", gsub(" ", "_", n)), ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})
 
# Find overlap in GO terms
go_lists <- lapply(go_tables, function(t){
  t$Description
})
v <- venn(go_lists)
attributes(v)
