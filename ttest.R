setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)
source("../pd_braak/PD/t.test.table.R")


donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")
samples_C <- readRDS("resources/samples_C.rds")
samples_D <- readRDS("resources/samples_D.rds")

# check order of samples
sapply(donorNames, function(d){
  data.frame(
  C = identical(colnames(brainExpr[[d]]), samples_C[[d]]$Sample_id),
  D = identical(colnames(brainExpr[[d]]), samples_D[[d]]$Sample_id)
  )
})

# T-test between network and whole brain/cortex?
ttest_C <- lapply(donorNames, function(d){
  expr <- brainExpr[[d]]
  samples <- samples_C[[d]]$Inside.y.n
  expr_in <- expr[, samples==1]
  expr_out <- expr[, samples==0]
  t.test.table(expr_in, expr_out)
})
ttest_C <- simplify2array(ttest_C)
saveRDS(ttest_C, file = "ttest_C.rds" )

ttest_D <- lapply(donorNames, function(d){
  expr <- brainExpr[[d]]
  samples <- samples_D[[d]]$Inside.y.n
  expr_in <- expr[, samples==1]
  expr_out <- expr[, samples==0]
  t.test.table(expr_in, expr_out)
})
ttest_D <- simplify2array(ttest_D)
saveRDS(ttest_D, file = "ttest_D.rds")

# Print number of diff. expr. genes
df <- data.frame(
  'Donors' = gsub("donor", "Donor ", donorNames),
  'Network_C' = apply(ttest_C, 3, function(x){sum(x[, "BH"] < 0.05 & abs(x[, "meanDiff"]) > 1)}),
  'Network_D' = apply(ttest_D, 3, function(x){sum(x[, "BH"] < 0.05 & abs(x[, "meanDiff"]) > 1)})
)
write.table(df, file = "number_of_degs.txt", row.names = FALSE, sep = "\t", quote = FALSE)
