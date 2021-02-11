# Select samples within SCNs

# Function to read mapping of AHBA samples to network
read.network <- function(nw, d){
  l <- lapply(c(1:length(d)), function(i){
    col_idx <- seq(i*5-4, i*5-1)
    tab <- nw[, col_idx]
    colnames(tab) <- colnames(nw)[1:4]
    tab$Sample_id <- gsub("'", "", tab$Sample_id)
    na_rows <- which(is.na(tab$Slab_num))
    tab <- if (length(na_rows) == 0) tab else tab[-na_rows, ]
    tab
  }) 
  names(l) <- d
  l
}

# Read and merge info of all networks
networks <- paste0("Network_", LETTERS[1:9])
networks <- sapply(networks, function(n){
  file <- paste0("../Results_Oleh/", n, ".txt")
  f <- read.delim(file, skip = 1)
  read.network(f, donorNames)
}, simplify = FALSE)
sample_info <- lapply(donorNames, function(d){
  shared <- networks[[1]][[d]][, c(1:3)] # info shared across networks
  nw_info <- sapply(networks, function(n){
    n[[d]][, -c(1:3)]
  })
  cbind(shared, nw_info)
})
# saveRDS(sample_info, file = "output/sample_info.rds")

# Number of samples within network
df <- data.frame(
  t(sapply(sample_info, function(t){ apply(t[, -c(1:3)], 2, sum)}))
)
df <- rbind(df, 'Total' = apply(df, 2, sum))
df <- cbind('Donors' = gsub("donor", "Donor ", rownames(df)), df)
write.table(df, "output/number_of_samples.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Function to select all sample IDs of a brain structure
sample.ids <- function(roi){
  id <- ontology[ontology$name == roi, "id"]
  rows <- grep(id, ontology$structure_id_path)
  ontology$id[rows]
}
# Get column indices of roi per donor
sample.idx <- function(id){
  lapply(brainExpr, function(e){
    colnames <- colnames(e)
    ids <- intersect(id, colnames)
    cols <- as.numeric(colnames %in% ids)
    # names(cols) <- colnames[cols]
    cols # column indices
  })
}

# Percentage cortical samples
cc <- sample.idx(sample.ids("cerebral cortex"))
network_idx <- lapply(networks, function(x) lapply(x, function(y) y $Inside.y.n))
network_idx$abefghi <- lapply(donorNames, function(d) { # non PD-related networks
  v <- lapply(networks[-c(3,4)], function(x) {
    x[[d]]$Inside.y.n
  })
  Reduce(bitwOr, v)
})
network_idx$whole_brain <- lapply(network_idx$abefghi, function(x){ # all AHBA samples (whole brain)
  rep(1, length(x))
})
size <- sapply(network_idx, function(x){
  sapply(donorNames, function(d){
    n <- x[[d]] # samples in network
    c <- cc[[d]] # samples in cortex
    i <- bitwAnd(n, c) # both present in cortex and network
    # sum(i)
    sum(i)/sum(n)*100
  })
})
size

#####
# Create phenotypic vector for Fulcher method
network_phenotype <- t(sapply(networks, function(l){
  l <- lapply(l, function(d){
    d$Inside.y.n
  })
  l <- unlist(l, use.names = F)
}))
write.csv(network_phenotype, "output/network_phenotype.csv")
