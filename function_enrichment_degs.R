# Functional enrichment of DEGs in network C and D
  library("RDAVIDWebService")

# RDavid
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl",
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
setTimeOut(david, 200000)
bg <- addList(david, ahba.genes(), idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg

output_dir <- "output/Functional_enrichment/"
mkdirs(output_dir)

lapply(names(degs), function(n){
  lapply(names(degs[[n]]), function(l){
    name <- paste(n, l, sep = "_")
    genes <- rownames(degs[[n]][[l]])
    result <- addList(david, genes, idType = "ENTREZ_GENE_ID", listName = name, listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    getFunctionalAnnotationChartFile(david, paste0(output_dir, name, "_goterms.txt"), threshold=0.05, count=2L)
    getClusterReportFile(david, paste0(output_dir, name, "_termclusters.txt"), type = c("Term"))
  })
})

go <- sapply(names(degs), function(n){
  sapply(names(degs[[n]]), function(l){
    name <- paste(n, l, sep = "_")
    fName <- paste0(output_dir, name, "_goterms.txt")
    print(fName)
    t <- read.csv(fName, header = TRUE, sep = "\t", colClasses = "character")
    rows <- t$Benjamini < 0.05
    t <- t[rows, c("Term", "Count", "Benjamini")]
    t[order(t$Benjamini), ]
  }, simplify = FALSE)
}, simplify = FALSE)
go <- unlist(go, recursive = FALSE)
go <- go[-which(sapply(go, nrow) == 0)]

lapply(names(go), function(n){
  t <- go[[n]]
  write.table(t, file = paste0(output_dir, "/BHcorrected_", n, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
})
