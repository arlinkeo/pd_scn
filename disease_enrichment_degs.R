# Disease enrichment
# library('disgenet2r')
library(reshape2)
library(DescTools)
# 
# res <- gene2disease(as.numeric(ll_degs$Network_C.upregulated), database = "ALL", verbose = TRUE)
# plot(res, class = "Network", prop = 10)
# plot(res, class = "Heatmap", limit = 100)
# 
# 
# res <- gene2disease(c("NPPA", "SOSTDC1", "TYRP1") , database = "ALL", verbose = TRUE)
# plot(res, class = "Network", prop = 10)
# plot(res, class = "DiseaseClass", limit = 200)

# Disease-associated genes
disease_table <- read.delim("../all_gene_disease_associations.tsv")
diseases <- unique(disease_table$diseaseName)
disease_genes <- sapply(diseases, function(x){
  disease_table$geneId[disease_table$diseaseName %in% x]
}, simplify = FALSE)
saveRDS(disease_genes, "output/disease_genes.rds")
disease_genes <- readRDS("output/disease_genes.rds")

# Disease enrichment of DEGs
or <- odd.ratio.table(ll_degs[3:4], disease_genes)
or <- round(or, digits = 2)
pval <- hyper.test.table(ll_degs[3:4], disease_genes)
rows <- apply(pval, 1, function(x) any(x < 0.05))
pval <- -log10(pval)
t <- abind('OR' = or[rows, ], 'P-value' =  pval[rows, ], along = 3, use.anon.names = TRUE)
dimnames(t)[[2]] <- gsub(".upregulated", "", dimnames(t)[[2]])
dimnames(t)[[2]] <- gsub("_", " ", dimnames(t)[[2]])
rowOrder <- order(apply(t[,,"or"], 1, max))
t <- t[rowOrder, ,]
df <- melt(t)
colnames(df) <- c("disease", "network", "measure", "value")
df$disease <- factor(df$disease, levels = rev(unique(df$disease)))
p <- ggplot(df, aes(x=disease, y=value, fill=network)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0), 
        legend.title = element_blank(),
        axis.title = element_blank()) +
  facet_grid(measure ~., scales = "free")
pdf("output/disease_enrichment.pdf", 11, 5)
p
dev.off()

# Same heat colors for all heatmaps
# t_log <- t(t_log)
# q1 <- quantile(t_log, 0.95)
# col_fun <- colorRamp2(c(0, q1), c("#EEEEEE", "red"))
# Heatmap(t_log, name = "-log10\nBH-corrected\nP-value",
#         col = col_fun,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         column_names_side = c("top"), 
#         row_names_side = c("left"),
#         width = unit(ncol(t_log)*.8, "lines"), 
#         height = unit(nrow(t_log)*.8, "lines")
#         )
