# Disease enrichment
library(reshape2)
library(DescTools)

# Disease-associated genes
disease_table <- read.delim("../all_gene_disease_associations.tsv")
diseases <- unique(disease_table$diseaseName)
disease_genes <- sapply(diseases, function(x){
  disease_table$geneId[disease_table$diseaseName %in% x]
}, simplify = FALSE)
saveRDS(disease_genes, "output/disease_genes.rds")
# disease_genes <- readRDS("output/disease_genes.rds")

# Disease enrichment of DEGs
or <- odd.ratio.table(ll_degs[3:4], disease_genes)
or <- round(or, digits = 2)
pval <- hyper.test.table(ll_degs[3:4], disease_genes)
rows <- apply(pval, 1, function(x) any(x < 0.05))
pval <- -log10(pval)
t <- abind('OR' = or[rows, ], 'P-value' =  pval[rows, ], along = 3, use.anon.names = TRUE)
dimnames(t)[[2]] <- gsub(".upregulated", "", dimnames(t)[[2]])
dimnames(t)[[2]] <- gsub("_", " ", dimnames(t)[[2]])
rowOrder <- order(apply(t[,,"OR"], 1, max))
t <- t[rowOrder, ,]
df <- melt(t)
colnames(df) <- c("disease", "network", "measure", "value")
df$disease <- factor(df$disease, levels = rev(unique(df$disease)))
p <- ggplot(df, aes(x=disease, y=value, fill=network)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0), 
        legend.title = element_blank(),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  facet_grid(measure ~., scales = "free")
pdf("output/disease_enrichment.pdf", 11, 5)
p
dev.off()