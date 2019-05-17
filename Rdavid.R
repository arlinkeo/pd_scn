setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)
library("RDAVIDWebService")
brainExpr <- readRDS("../AHBA_Arlin/gene_expr.RDS")

# Functional enrichment of DEGs in network C and network D
degs <- readRDS("resources/degs.rds")

#Functional enrichment of  genes correlated greater or smaller than 0
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
bg_list <- rownames(brainExpr$donor9861)
rm(brainExpr)
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold

# Enrichment of positively and negatively correlated progression genes
lapply(names(degs), function(n){
  lapply(names(degs[[n]]), function(l){
    name <- paste(n, l, sep = "_")
    genes <- rownames(degs[[n]][[l]])
    result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = name, listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    getFunctionalAnnotationChartFile(david, paste0("Functional_analyses/", name, "_goterms.txt"), threshold=t, count=2L)
    getClusterReportFile(david, paste0("Functional_analyses/", name, "_termclusters.txt"), type = c("Term"))
  })
})
Sys.time()


