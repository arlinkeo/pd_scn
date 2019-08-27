# Thalamus analysis

# Thalamic samples within networks

thalamus_id <- ontology[ontology$name == "thalamus", "id"]

nws <- c("Network_C", "Network_D")

l <- sapply(nws, function(n){
  s <- lapply(donorNames, function(d){
    info <- networks[[n]][[d]]
    s <- info[info$Inside.y.n == 1, "Sample_id"]
    ontology_rows <- match(s, ontology$id)
    t <- ontology[ontology_rows, c("acronym", "name", "structure_id_path")]
    t[grep(thalamus_id, t$structure_id_path), ]
  })
  s <- melt(s)
  colnames(s)[4] <- "donor"
  s <- s[order(s$structure_id_path),]
  
  parent_structures <- lapply(s$structure_id_path, function(s){
    s <- unlist(strsplit(s, "/"))
    s <- s[-c(1:which(s == thalamus_id))]
    s <- ontology[ontology$id %in% s, "name"]
    paste(s, collapse = "/")
  })
  s$structure_name_path <- parent_structures
  s[, -c(2,3)]
}, simplify = FALSE)
