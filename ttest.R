setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)
library(plyr)
library(metafor)
source("../pd_braak/PD/t.test.table.R")

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
sample_info <- readRDS("resources/sample_info.rds")

# check order of samples
sapply(donorNames, function(d){ 
  identical(colnames(brainExpr[[d]]), sample_info[[d]]$Sample_id)
})

# Differential expression between network and rest of the brain (t-test)
ttest <- lapply(donorNames, function(d){
  expr <- brainExpr[[d]]
  samples <- sample_info[[d]][, -c(1:3)]
  tab <- alply(samples, 2, function(s){
    expr_in <- expr[, s==1] # samples in one network
    expr_out <- expr[, s==0]# samples outside one network
    t.test.table(expr_out, expr_in) # in compared to out
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

# Print number of diff. expr. genes
df <- t(apply(summary_ttest, c(1,2), function(x){
  sum(x[, "BH"] < 0.05 & abs(x[, "Estimate"]) > 1)
}))
write.table(df, file = "number_of_degs.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Get tables of DEGs
degs <- alply(summary_ttest[, "summary", , ], 1, function(x){
  down <- x[which(x[, "Estimate"] < -1 & x[,"BH"] < 0.05),]
  down <- down[order(down[,"BH"]),]
  up <- x[which(x[, "Estimate"] > 1 & x[,"BH"] < 0.05),]
  up <- up[order(up[,"BH"]),]
  list(down = down, up = up)
}, .dims = TRUE)
saveRDS(degs, file = "resources/degs.rds")