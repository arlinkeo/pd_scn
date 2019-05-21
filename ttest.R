setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)
library(plyr)
library(metafor)
library(venn)
source("../pd_braak/PD/t.test.table.R")

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
sample_info <- readRDS("resources/sample_info.rds")

# Differential expression between network and rest of the brain (t-test)
ttest <- lapply(donorNames, function(d){
  expr <- brainExpr[[d]]
  samples <- sample_info[[d]][, -c(1:3)]
  abefghi <- samples[, -which(colnames(samples) %in% c("Network_C", "Network_D"))] # Networks unaffected in PD
  abefghi <- apply(abefghi, 1, function(x) Reduce(bitwOr, x))
  cd <- samples[, c("Network_C", "Network_D")] # Networks affected in PD
  tab <- alply(cd, 2, function(s){
    a <- expr[, s==1] # expression in network C or D
    b <- expr[, abefghi==1]# expression in network ABEFGHI
    t.test.table(a, b) # in compared to out
  }, .dims = TRUE) # keep names
  simplify2array(tab) # 3D array: genes x measures x networks
})
ttest <- simplify2array(ttest) # 4D array: genes x measures x networks x donors
saveRDS(ttest, file = "resources/ttest.rds")
ttest <- readRDS("resources/ttest.rds")

# Meta-analysis of differential expression across donors
summary_ttest <- aaply(ttest, c(1,3), function(g){ # For each Braak region pair and gene
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
saveRDS(summary_ttest, file = "resources/summary_ttest.rds")
summary_ttest <- readRDS("resources/summary_ttest.rds")

# # Print number of diff. expr. genes per donor and summary
# df <- aaply(summary_ttest, c(1,2), function(x){
#   down <- x[which(x[, "Estimate"] < -1 & x[,"BH"] < 0.05),]
#   up <- x[which(x[, "Estimate"] > 1 & x[,"BH"] < 0.05),]
#   n1 <- nrow(down)
#   n2 <- nrow(up)
#   n1 <- ifelse(is.null(n1), 0, n1)
#   n2 <- ifelse(is.null(n2), 0, n2)
#   c(down = n1, up = n2)
# })
# df <- alply(df, 1, function(x)x)
# df <- Reduce(cbind, df)
# colnames(df) <- c("Downregulated in network C", "Upregulated in network C", "Downregulated in network D", "Upregulated in network D")
# df <- cbind(Donor = gsub("donor", "Donor ", rownames(df)), df)

# Get tables of DEGs
degs <- alply(summary_ttest[, "summary", , ], 1, function(x){
  down <- x[which(x[, "Estimate"] < -1 & x[,"BH"] < 0.05),]
  down <- down[order(down[,"BH"]),]
  up <- x[which(x[, "Estimate"] > 1 & x[,"BH"] < 0.05),]
  up <- up[order(up[,"BH"]),]
  list(downregulated = down, upregulated = up)
}, .dims = TRUE)
saveRDS(degs, file = "resources/degs.rds")

# Print number of diff. expr. genes
df <- sapply(degs, function(l){
  sapply(l, nrow)
})
write.table(df, file = "number_of_degs.txt", sep = "\t", quote = FALSE)

# Overlap of DEGs between network C and D
ll_degs <- lapply(degs, function(l) lapply(l, rownames))
ll_degs <- unlist(ll_degs, recursive = FALSE)[c(1,3,2,4)]
names(ll_degs) <- gsub("_", " ", names(ll_degs))
pdf("venn_degs.pdf", 6, 5)
venn <- venn(ll_degs, ellipse = TRUE, zcolor = "style", cexil = 1.2, cexsn = 1)
dev.off()