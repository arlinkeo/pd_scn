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
saveRDS(sample_info, file = "output/sample_info.rds")

# Number of samples within network
df <- data.frame(
  t(sapply(sample_info, function(t){ apply(t[, -c(1:3)], 2, sum)}))
)
df <- rbind(df, 'Total' = apply(df, 2, sum))
df <- cbind('Donors' = gsub("donor", "Donor ", rownames(df)), df)
write.table(df, "output/number_of_samples.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# # Check overlap of samples with stress network
# scz_network <- lapply(donorNames, function(d){
#   read.table(file = paste0("C:/Users/dkeo/surfdrive/MandySCZpaper/", d, ".txt"), header = TRUE, sep ="\t")
# })
# sampleOverlap <- lapply(donorNames, function(d){
#   sapply(networks, function(n){
#     scz <- scz_network[[d]]$Inside.Y.N
#     scn <- n[[d]]$Inside.y.n
#     df <- data.frame(scz, scn, both = bitwAnd(scn, scz))
#     apply(df, 2, sum)
#   })
# })
