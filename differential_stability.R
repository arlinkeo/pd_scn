# Differential Stability
library(reshape2)
library(ggplot2)

pairs <- combn(donorNames, 2) # pairwise donor combinations

r <- apply(pairs, 2, function(x){
  d1 <- x[1]
  d2 <- x[2]
  e1 <- brainExpr[[d1]]
  e2 <- brainExpr[[d2]]
  c <- intersect(colnames(e1), colnames(e2))
  e1 <- e1[, c]
  e2 <- e2[, c]
  genes <- rownames(e1)
  sapply(genes, function(g){
    cor(t(e1[g, c]), t(e2[g, c]))
  })
})
ds <- data.frame(
  ds = apply(r, 1, mean)
)
ds$rank <- rank(ds$ds)
saveRDS(ds, "output/differential_stability.rds")

ds$label <- sapply(rownames(ds), function(x){
  if (any(unlist(ll_degs) %in% x)) entrezId2Name(x) else ""
})
ds$color <- ifelse(ds$label == "", "non-DEG", "DEG")
ds <- ds[ds$rank, ]
head(ds)
pdf("output/differential_stability_degs.pdf", 8, 3)
ggplot(ds, aes(rank, ds, color = color)) +
  geom_point(size = 0.25) +
  ylab("Differential Stability") +
  geom_vline(xintercept = floor(length(ahba.genes())*0.9), color = "gray") +
  # geom_text(aes(label=label),hjust=0, vjust=0) +
  theme_classic() +
  theme(legend.title = element_blank())
dev.off()

topdecile <- sort(ds$ds, decreasing = TRUE)[floor(length(ahba.genes())*0.1)]
sapply(ll_degs, function(g){
  sum <- sum(ds[g, "ds"] > topdecile)
  print(paste0(sum, "/", length(g), " (", round(sum/length(g)*100, digits = 1), "%)"))
}, simplify = FALSE)

# Write tables of DEGs with DS
lapply(names(degs), function(network){
  direction <- degs[[network]]
  lapply(names(direction), function(d){
    t <- degs[[network]][[d]]
    t <- t[, -c(2,6)]
    t <- data.frame(
      'Gene name'= entrezId2Name(rownames(t)), 
      'Entrez ID' = rownames(t), 
      t,
      'Differential stability' = ds[rownames(t), "ds"],
      check.names = FALSE)
    t[, c(3:5,8)] <- round(t[, c(3:5,8)], digits = 2)
    t[, c(6:7)] <- format(t[, c(6:7)], digits = 3, scientific = T)
    colnames(t) <- paste0(toupper(substr(colnames(t), 1, 1)), substring(colnames(t), 2))
    colnames(t)[c(3,6)] <- c("Fold-change", "P-value")
    # t[is.na(t)] <- ""
    write.table(t, file = paste0("output/degs_ds_", network,"_", d, ".txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  })
})
