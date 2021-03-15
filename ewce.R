# EWCE
# Based on tutorial: https://nathanskene.github.io/EWCE/articles/EWCE.html

library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
theme_set(theme_cowplot())
library(reshape2)

################################################################################
# Prepare Cell type data

# Load NeuroExpresso annotation data
load("../Neuroexpresso_expression/n_expressoSamples.rda")
n_expressoSamples <- n_expressoSamples
apply(n_expressoSamples[, c(5:11)], 2, function(x) data.frame(table(x)))
n_expressoSamples <- n_expressoSamples[, c("sampleName", "MajorType", "CellTypes")] # keep only expression info
colnames(n_expressoSamples) <- c("cell_id", "level1class", "level2class")
rownames(n_expressoSamples) <- n_expressoSamples$cell_id

# Load NeuroExpresso data (not single cell)
load("../Neuroexpresso_expression/n_expressoExpr.rda")
ag <- n_expressoExpr$Gene.Symbol[which(duplicated(n_expressoExpr$Gene.Symbol))] # 69 ambiguous, duplicate genes
sum(n_expressoExpr$Gene.Symbol %in% ag) #
n_expressoExpr <- n_expressoExpr[!(n_expressoExpr$Gene.Symbol %in% ag), ] # remove
rownames(n_expressoExpr) <- n_expressoExpr$Gene.Symbol
n_expressoExpr <- as.matrix(n_expressoExpr[, -c(1:6)])

# There are more samples in annotation data, select intersection
samples <- intersect(colnames(n_expressoExpr), rownames(n_expressoSamples))
n_expressoExpr <- n_expressoExpr[, samples]

# Check expression of genes in cell-type data
gene <- c("Nppa", "Sostdc1", "Tyrp1")
cellExpDist = data.frame(t(n_expressoExpr[gene,]),l2=n_expressoSamples[colnames(n_expressoExpr),]$level2class)
cellExpDist = melt(cellExpDist, value.name = "e")
ggplot(cellExpDist) + geom_boxplot(aes(x=l2,y=e)) + xlab("Cell type") + ylab("Expression") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(variable~., scales = "free")

# Copy function from https://github.com/NathanSkene/EWCE/blob/fd47e07dbe12b27217d06bf2e558923eb941f6c0/R/generate.celltype.data.r
calculate.specificity.for.level <- function(ctd_oneLevel) {
  normalised_meanExp <- t(t(ctd_oneLevel$mean_exp) * 
                            (1 / colSums(ctd_oneLevel$mean_exp)))
  ctd_oneLevel$specificity <- normalised_meanExp / 
    (apply(normalised_meanExp, 1, sum) + 0.000000000001)
  return(ctd_oneLevel)
}

# Generate specificity matrix with only informative genes (cell type markers)
markers <- lapply(conversion_table, function(t) t$'10090')
markers <- unlist(markers) # 647 markers in total
markers <- intersect(markers, rownames(n_expressoExpr)) # 629 markers? markers from only purified bulk cell-type data
exp_dropped <- n_expressoExpr[markers,]
celltypedata <- list(
  level1class =list(
    annot = n_expressoSamples$level1class
  ),
  level2class = list(
    annot = n_expressoSamples$level2class
  )
)
celltypedata <- lapply(names(celltypedata), function(n){
  ct <- celltypedata[[n]]
  ct$mean_expr <- apply(exp_dropped, 1, function(g){
    df <- data.frame(g, n_expressoSamples[, n])
    colnames(df)[2] <- n
    l <- split(df, df[, n])
    sapply(l, function(t)mean(t$g))
  })
  ct <- calculate.specificity.for.level(ct)
  ct$mean_expr <- t(ct$mean_expr)
  ct$specificity <- t(ct$specificity)
  ct
})

################################################################################
# EWCE with cell-types defined by Zeisel et al. (already processed by Skene et al.)

set.seed(1234)
library(EWCE)
library(homologene)

h2m <- homologene(ahba.genes(), inTax = 9606, outTax = 10090)
mouse.bg <- h2m[, "10090"]

# EWCE
ewce <- sapply(c(1:2), function(l){
  ewce <- lapply(degs, function(nw){
    lapply(nw, function(g){
      mouse.hits <- h2m[h2m$'9606_ID' %in% rownames(g), "10090"]
      full_results = bootstrap.enrichment.test(sct_data=celltypedata,hits=mouse.hits,bg=mouse.bg,
                                               reps=1000,annotLevel=l)
      full_results
    })
  })
  ewce <- unlist(ewce, recursive = FALSE)
  ewce <- sapply(names(ewce), function(x){
    data.frame(ewce[[x]]$results, list = x)
  }, simplify = FALSE)
  Reduce(rbind, ewce)
}, simplify = FALSE)

# plot results
lapply(c(1:2), function(l){
  width <- if(l==1) 4 else 8
  pdf(paste0("output/ewce_level", l, ".pdf"), width, 6)
  e <- ewce[[l]]
  e <- e[grep("upregulated", e$list), ]
  e$list <- gsub(".upregulated", "", e$list)
  e$list <- gsub("_", " ", e$list)
  p <- print(ewce.plot(total_res=e,mtc_method="BH")$plain)
  dev.off()
  p$data
})

lapply(ewce, function(t){
  # t$bh <- p.adjust(t$p, method = "BH")
  t[t$p<0.05, ]
})

# ##### Using Zeisel's single cell data #####
# 
# # drop genes with no mouse orthologs
# data("mouse_to_human_homologs")
# m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
# mouse.bg  = unique(m2h$MGI.symbol)
# 
# # EWCE
# ewce <- sapply(c(1:2), function(l){
#   ewce <- lapply(degs, function(nw){
#     lapply(nw, function(g){
#       g <- entrezId2Name(rownames(g))
#       print(unique(m2h[m2h$HGNC.symbol %in% g,"HGNC.symbol"]))
#       mouse.hits = unique(m2h[m2h$HGNC.symbol %in% g,"MGI.symbol"])
#       print(mouse.hits)
#       full_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,bg=mouse.bg,
#                                                reps=1000,annotLevel=l)
#       
#     })
#   })
#   ewce <- unlist(ewce, recursive = FALSE)
#   ewce <- sapply(names(ewce), function(x){
#     data.frame(ewce[[x]]$results, list = x)
#   }, simplify = FALSE)
#   Reduce(rbind, ewce)
# }, simplify = FALSE)