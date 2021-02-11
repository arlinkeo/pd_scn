# EWCE
# Based on tutorial: https://nathanskene.github.io/EWCE/articles/EWCE.html

library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
theme_set(theme_cowplot())

################################################################################
# Prepare SCT data

# Load single cell transcriptomic data from cortex
# download.file(url ="https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt",
#               destfile="expression_mRNA_17-Aug-2014.txt") 
path = "expression_mRNA_17-Aug-2014.txt"
cortex_mrna  = load.linnarsson.sct.data(path)
data(cortex_mrna)

# Check SCT data
gene <- "Snca"
cellExpDist = data.frame(e=cortex_mrna$exp[gene,],l1=cortex_mrna$annot[colnames(cortex_mrna$exp),]$level1class)
ggplot(cellExpDist) + geom_boxplot(aes(x=l1,y=e)) + xlab("Cell type") + ylab("Unique Molecule Count") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

# Drop genes to reduce compute time
nKeep = 1000
must_keep = c("Apoe","Gfap","Gapdh")
set.seed(123458)
keep_genes = c(must_keep,sample(rownames(cortex_mrna$exp),997))
cortex_mrna$exp = cortex_mrna$exp[keep_genes,]

# Transform data to normalize for differences in cell size
library(sctransform)
scT = sctransform::vst(cortex_mrna$exp, return_cell_attr = TRUE)
cortex_mrna$exp_scT = correct_counts(scT, cortex_mrna$exp) # umi_corrected
cortex_mrna$exp_scT_normed = Matrix::t(Matrix::t(cortex_mrna$exp_scT)*(1/Matrix::colSums(cortex_mrna$exp_scT)))

# Generate celltype data for just the cortex/hippocampus data
exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=cortex_mrna$exp_scT_normed,level2annot = cortex_mrna$annot$level2class)
annotLevels = list(level1class=cortex_mrna$annot$level1class,level2class=cortex_mrna$annot$level2class)
fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,groupName="kiCortexOnly")
print(fNames_CortexOnly)
fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
print(fNames_CortexOnly)
load(fNames_CortexOnly[1])

# # Download the hypothalamus data and unzip
# if(!file.exists("GSE74672_expressed_mols_with_classes.xlsx")){
#   download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz", destfile="GSE74672_expressed_mols_with_classes.xlsx.gz")
#   system("gunzip GSE74672_expressed_mols_with_classes.xlsx.gz")
# }
# 
# # Read in the hypothalamus data
# hypo_dat = read_excel("hypoth_moldata_classification08-Mar-2017.xlsx")# name of unzipped file
# 
# # Extract the expression data, gene symbols and annotation data
# exp = data.matrix(hypo_dat[12:dim(hypo_dat)[1],2:dim(hypo_dat)[2]])
# rownames(exp) = data.frame(hypo_dat[12:dim(hypo_dat)[1],1])[,1]
# level1class = data.frame(level1class=t(hypo_dat[1,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
# level2class = data.frame(leve2class=t(hypo_dat[2,2:dim(hypo_dat)[2]]),stringsAsFactors = FALSE)[,1]
# cell_id     = colnames(hypo_dat)[2:dim(hypo_dat)[2]]
# hypo_annot  = data.frame(cell_id=cell_id,level1class=level1class,level2class=level2class,stringsAsFactors = FALSE)
# 
# # Drop the glia and unclassified cells (which don't have level 2  annotations)
# hypo_annot  = hypo_annot[!is.na(hypo_annot$level2class) & !hypo_annot$level2class=="uc",]
# hypo_exp    = exp[,hypo_annot$cell_id]
# 
# # Make the celltype names more aesthetically pleasing
# hypo_annot$level2class=gsub(",",";",hypo_annot$level2class)
# hypo_annot$level1class[grep("Oxt;|^Avp",hypo_annot$level2class)] = "Oxytocin / Vasopressin Expressing Neurons"
# hypo_annot$level1class[grep("^Th;|^Dopamine",hypo_annot$level2class)] = "Hypothalamic Dopaminergic Neurons"
# hypo_annot$level1class[grepl("^Vglut2|^Trh|^Qrfp|^Hcrt|^Pmch|^Adcyap1|^Npvf|^Ghrh|^Hmit|^Nms|^Vip;|^Per2|Tnr$|^Gad-low;Gnrh",hypo_annot$level2class) & grepl("neurons",hypo_annot$level1class)] = "Hypothalamic Glutamatergic Neurons"
# hypo_annot$level1class[grepl("GABA|^Sst|^Crh|^Npy|^Pomc|^Galanin|^Otof|Pnoc$|^Calcr-high",hypo_annot$level2class) & grepl("^neurons$",hypo_annot$level1class)] = "Hypothalamic GABAergic Neurons"
# hypo_annot$level2class[hypo_annot$level2class!=""] = sprintf("Hypothalamic %s Neuron",hypo_annot$level2class[hypo_annot$level2class!=""])
# 
# # Fix bad MGI symbols
# hypo_exp_CORRECTED = fix.bad.mgi.symbols(hypo_exp)




# map AHBA genes with SCT genes
sctGenes <- toupper(rownames(cortex_mrna$exp))
gene.map <- data.frame(AHBA = entrezId2Name(ahba.genes()))
x=match(sctGenes, gene.map$AHBA)
length(x[!is.na(x)]) # 13,667 overlapping genes between two datasets
gene.map$SCT <- NA
gene.map$SCT[x[!is.na(x)]] <- rownames(cortex_mrna$exp)[!is.na(x)]
rownames(gene.map) <- ahba.genes()

################################################################################
# load("../Neuroexpresso_expression/n_expressoExpr.rda")
# load("../Neuroexpresso_expression/n_expressoSamples.rda")
# load("../Neuroexpresso_expression/TasicMouseExp.rda")
# load("../Neuroexpresso_expression/TasicMouseMeta.rda")
# 
# annot <- n_expressoSamples$CellTypes
# neuroExpr <- n_expressoExpr
# rownames(neuroExpr) <- neuroExpr$Gene.Symbol
# rownames(neuroExpr) <- n_expressoExpr$NCBIids
# 
# # Generate celltype data for just the cortex/hippocampus data
# neuroExpr = drop.uninformative.genes(exp=neuroExpr,level1annot = annot)
# annotLevels = list(level1class=cortex_mrna$annot$level1classs)
# fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,groupName="kiCortexOnly")
# print(fNames_CortexOnly)
# fNames_CortexOnly = filter.genes.without.1to1.homolog(fNames_CortexOnly)
# print(fNames_CortexOnly)
# load(fNames_CortexOnly[1])
# 


################################################################################
data("ctd")
set.seed(1234)

library(EWCE)

# drop genes with no mouse orthologs
data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
#mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
mouse.bg  = unique(m2h$MGI.symbol)

# EWCE
ewce <- sapply(c(1:2), function(l){
  ewce <- lapply(degs, function(nw){
    lapply(nw, function(g){
      g <- entrezId2Name(rownames(g))
      print(unique(m2h[m2h$HGNC.symbol %in% g,"HGNC.symbol"]))
      mouse.hits = unique(m2h[m2h$HGNC.symbol %in% g,"MGI.symbol"])
      print(mouse.hits)
      full_results = bootstrap.enrichment.test(sct_data=ctd,hits=mouse.hits,bg=mouse.bg,
                                               reps=1000,annotLevel=l)
      
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
  width <- if(l==1) 8 else 11
  pdf(paste0("output/ewce_level", l, ".pdf"), width, 9)
  print(ewce.plot(total_res=ewce[[l]],mtc_method="BH")$plain)
  dev.off()
})

