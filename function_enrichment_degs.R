# Functional enrichment of DEGs in network C and D
library("RDAVIDWebService")

# RDavid
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
bg_list <- rownames(brainExpr$donor9861)
rm(brainExpr)
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
t <- 0.05 # EASE p-value threshold

lapply(names(degs), function(n){
  lapply(names(degs[[n]]), function(l){
    name <- paste(n, l, sep = "_")
    genes <- rownames(degs[[n]][[l]])
    result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = name, listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    getFunctionalAnnotationChartFile(david, paste0("output/Functional_analyses/", name, "_goterms.txt"), threshold=t, count=2L)
    getClusterReportFile(david, paste0("output/Functional_analyses/", name, "_termclusters.txt"), type = c("Term"))
  })
})

#Function to read Rdavid output
read.RdavidOutput <- function(fileName){
  if (file.exists(fileName)){
    terms <- read.csv(fileName, header = TRUE, sep = "\t", colClasses = "character")
  } else {
    NULL
  }
}

go <- sapply(names(degs), function(n){
  sapply(names(degs[[n]]), function(l){
    name <- paste(n, l, sep = "_")
    fName <- paste0("output/Functional_analyses/", name, "_goterms.txt")
    print(fName)
    t <- read.RdavidOutput(fName)
    rows <- t$Benjamini < 0.05
    t <- t[rows, c("Term", "Count", "Benjamini")]
    t[order(t$Benjamini), ]
  }, simplify = FALSE)
}, simplify = FALSE)
go <- unlist(go, recursive = FALSE)
go <- go[-which(sapply(go, nrow) == 0)]

lapply(names(go), function(n){
  t <- go[[n]]
  write.table(t, file = paste0("output/Functional_analyses/", n, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
})
