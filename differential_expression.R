# Differential expression between C or D and remaining networks A, B, E, F, G, H, and I
library('plyr')
library('metafor')
library('venn')
library('ggplot2')
library('ggrepel')

# Differential expression between network and rest of the brain (t-test)
ttest <- lapply(donorNames, function(d){
  expr <- brainExpr[[d]]
  samples <- sample_info[[d]][, -c(1:3)]
  abefghi <- samples[, -which(colnames(samples) %in% c("Network_C", "Network_D"))] # Networks unaffected in PD
  abefghi <- apply(abefghi, 1, function(x) Reduce(bitwOr, x))
  cd <- samples[, c("Network_C", "Network_D")] # Networks affected in PD
  tab <- alply(cd, 2, function(s){
    a <- expr[, abefghi==1]# expression in network ABEFGHI
    b <- expr[, s==1] # expression in network C or D
    t.test.table(a, b) # in compared to out
  }, .dims = TRUE) # keep names
  simplify2array(tab) # 3D array: genes x measures x networks
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x networks x donors
saveRDS(ttest, file = "output/ttest.rds")
# ttest <- readRDS("output/ttest.rds")

# Meta-analysis of differential expression across donors
summary_ttest <- aaply(ttest, c(1,3), function(g){
  gene <- t(g)
  t <- escalc(measure = "MD",  # Get estimates, variance (needed for meta-analysis) and confidence intervals (region B vs. A)
              m1i = gene[, "meanB"], m2i = gene[, "meanA"], # estimate
              n1i = gene[, "sizeB"], n2i = gene[, "sizeA"], 
              sd1i = sqrt(gene[, "varB"]), sd2i = sqrt(gene[, "varA"]))
  t <- summary(t)[, -c(3,4)]
  rownames(t) <- rownames(gene)
  colnames(t) <- c("Estimate", "Var", "lower95", "upper95")
  summary <- rma(t$Estimate, t$Var, method = "DL", test = "t")
  # gene <- cbind(gene, weight  = weights(summary))
  t <- cbind(t, pvalue = gene[, "pvalue"], weight = weights(summary))
  t <- rbind(t, 'summary' = list(summary$beta, summary$se^2 , summary$ci.lb, summary$ci.ub,
                            summary$pval, sum(weights(summary))))
  as.matrix(t)
}) # 4D-array: genes x networks x donors x measures

summary_ttest <- aaply(summary_ttest, c(2,3), function(t){ # P-value corrected for genes
  b <- p.adjust(t[, "pvalue"], method = "BH")
  cbind(t, BH = b)
}) # 4D-array: networks x donors x genes x measures (+BH)
saveRDS(summary_ttest, file = "output/summary_ttest.rds")
# summary_ttest <- readRDS("output/summary_ttest.rds")

# Print number of diff. expr. genes per donor and summary
df <- aaply(summary_ttest, c(1,2), function(x){
  down <- x[which(x[, "Estimate"] < -1 & x[,"BH"] < 0.05),]
  up <- x[which(x[, "Estimate"] > 1 & x[,"BH"] < 0.05),]
  n1 <- nrow(down)
  n2 <- nrow(up)
  n1 <- ifelse(is.null(n1), 0, n1)
  n2 <- ifelse(is.null(n2), 0, n2)
  c(down = n1, up = n2)
})
df <- alply(df, 1, function(x)x)
df <- Reduce(cbind, df)
colnames(df) <- c("Downregulated in network C", "Upregulated in network C", "Downregulated in network D", "Upregulated in network D")
rownames(df)[7] <- "Summary"
df <- cbind(Donor = gsub("donor", "Donor ", rownames(df)), df)
write.table(df, file = paste0("output/number_of_degs_per_donor.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Get tables of DEGs
degs <- alply(summary_ttest[, "summary", , ], 1, function(x){
  down <- x[which(x[, "Estimate"] < -1 & x[,"BH"] < 0.05),]
  down <- down[order(down[,"BH"]),]
  up <- x[which(x[, "Estimate"] > 1 & x[,"BH"] < 0.05),]
  up <- up[order(up[,"BH"]),]
  list(downregulated = down, upregulated = up)
}, .dims = TRUE)
# saveRDS(degs, file = "output/degs.rds")

# Write tables of DEGs
lapply(names(degs), function(network){
  direction <- degs[[network]]
  lapply(names(direction), function(d){
    t <- degs[[network]][[d]]
    t <- t[, -c(2,6)]
    t[, c(1:3)] <- round(t[, c(1:3)], digits = 2)
    t[, c(4:5)] <- format(t[, c(4:5)], digits = 3, scientific = T)
    t <- cbind('gene'= entrezId2Name(rownames(t)), 'gene_id' = rownames(t), t)
    colnames(t) <- paste0(toupper(substr(colnames(t), 1, 1)), substring(colnames(t), 2))
    colnames(t)[c(2,3,6)] <- c("Entrez ID", "Fold-change", "P-value")
    write.table(t, file = paste0("output/degs_", network,"_", d, ".txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  })
})

# Print number of diff. expr. genes
df <- sapply(degs, function(l) sapply(l, nrow))
rownames(df) <- paste0(toupper(substr(rownames(df), 1, 1)), substring(rownames(df), 2))
write.table(df, file = "output/number_of_degs.txt", sep = "\t", quote = FALSE)

# Overlap of DEGs between network C and D
ll_degs <- lapply(degs, function(l) lapply(l, rownames))
ll_degs <- unlist(ll_degs, recursive = FALSE)[c(1,3,2,4)]
names(ll_degs) <- gsub("_", " ", names(ll_degs))
pdf("output/venn_degs.pdf", 6, 5)
venn <- venn(ll_degs, ellipse = TRUE, zcolor = "style", cexil = 1.2, cexsn = 1)
dev.off()
overlap <- c(attributes(venn)$intersections$`Network C.downregulated:Network D.downregulated`, 
             attributes(venn)$intersections$`Network C.upregulated:Network D.upregulated`)

# # Check presence PD variant-associated genes
# pdGenes <- list(hiImpact = c("SNCA", "LRRK2", "GBA", "VPS35", "PARK2", "PINK1", "PARK7", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6", "SYNJ1", 
#                              "EIF4G1", "DNAJC13", "CHCHD2", "C20orf30", "RIC3", "LRP10"), #TMEM230 is C20orf30
#                 jansen2017 = c("INPP5F", "TMEM175", "ASH1L", "MAPT", "RIT1", "C14orf83", "STK39", "GPNMB", "BST1", 
#                                "SIPA1L2", "DLG2", "NUCKS1", "GCH1", "MCCC1", "FAM47E", "BCKDK", "TMPRSS9", "UBOX5", 
#                                "CCDC62", "SYNJ1", "EIF4G1", "FBXO7", "C20orf30", "POLG", "VPS13C", "PLA2G6"),
#                 hla = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1"),
#                 'Chang et al. 2017' = read.table("../../pd_braak/chang2017_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1], 
#                 'Nalls et al. 2014' = read.table("../../pd_braak/nalls2014_riskgenes.txt", comment.char = "#", sep = "\n", row.names = NULL, stringsAsFactors = FALSE)[, 1]
# )
# pdGenesID <- lapply(pdGenes, name2EntrezId)
# pdGenesID <- lapply(pdGenesID, function(x) x[!is.na(x)])
# lapply(pdGenesID, function(pd){
#   lapply(degs, function(network){
#     lapply(network, function(t){
#       intersect(rownames(t), pd)
#     })
#   })
# })

xlim <- max(summary_ttest[, "summary", , "Estimate"])
ylim <- c(0, max(-log10(summary_ttest[, "summary", , "pvalue"])))

summary <- summary_ttest[, "summary", , c("Estimate", "pvalue", "BH")]
volcanoplots <- lapply(dimnames(summary)[[1]], function(name){
  df <- data.frame(summary[name, ,])
  df$info <- ifelse(abs(df$Estimate) > 1 & df$BH < 0.05, 1, 0)
  df$info <- ifelse(rownames(df) %in% overlap, 2, df$info)
  df$info <- as.factor(df$info)
  df$logp <- -log10(df$pvalue)
  df$label <- entrezId2Name(rownames(df))
  df$label[-tail(order(abs(df$Estimate)), 10)] <- "" # labels of top 10
  cbf <- if (name == "Network_C") c("#999999", "#0072B2", "#F0E442") else c("#999999", "#CC79A7", "#F0E442") # color blind friendly pallette
  
  ggplot(df, aes(Estimate, logp, colour = info, label = label)) +
    geom_point(size = .1) +
    geom_text_repel(force = 5, colour = "black", size = 2.5, nudge_y = 0.1, 
                    fontface = "italic", segment.size = .1) +
    scale_colour_manual(values = cbf) + 
    scale_x_continuous(limits = c(-xlim, xlim)) +
    scale_y_continuous(limits = ylim) +
    labs(x = expression('-log'[2]*' '*'fold-change'), y = expression('-log'[10]*' '*italic('P')*'-value')) +
    ggtitle(network_names[name]) +
    theme_classic() + 
    theme(legend.position = "none", plot.title = element_text(size = 10))
})
pdf("output/volcanoplots.pdf", 4, 3)
volcanoplots
dev.off()