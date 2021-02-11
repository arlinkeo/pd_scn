# Functional enrichment with Reactome
library('ReactomePA')

# GSEA with Reactome PA
pathways <- lapply(ll_degs[3:4], function(g){
  reactome <- enrichPathway(g, universe = ahba.genes(), readable = TRUE)
  p <- reactome@result
  p <- p[p$p.adjust< 0.05, ]
  p
})

lapply(names(pathways), function(n){
  p <- pathways[[n]]
  p <- p[-c(1,3:5,7,8)]
  p$p.adjust <- format(p$p.adjust, digits = 3, scientific = TRUE)
  colnames(p) <- c("Pathway", "BH", "Gene count")
  write.table(p, file = paste0("output/reactomePA_", n, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
})
saveRDS(pathways, file = "output/reactomePA_pathways.rds")
# pathways <- readRDS("output/reactomePA_pathways.rds")

# # Read results from pathway enrichment
# pathways <- lapply(names(pathways), function(n){
#   p <- read.delim(file = paste0("output/reactomePA_", n, ".txt"))
# })

# overlap between C and D
intersect(pathways$`Network C.upregulated`$Description, pathways$`Network D.upregulated`$Description)