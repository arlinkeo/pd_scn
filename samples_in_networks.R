# Select samples within SCNs
setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn")
options(stringsAsFactors = FALSE)

donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

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

# Read and merge info of both networks C and D
samples_networkC <- read.delim("Results_Oleh/Network_C.txt", skip = 1)
samples_networkD <- read.delim("Results_Oleh/Network_D.txt", skip = 1)
samples_C <- read.network(samples_networkC, donorNames)
samples_D <- read.network(samples_networkD, donorNames)
sample_info <- lapply(donorNames, function(d){
  cbind(samples_C[[d]][, c(1:3)],
    'Network_C' = samples_C[[d]]$Inside.y.n, 'Network_D' = samples_D[[d]]$Inside.y.n)
})
saveRDS(sample_info, file = "resources/sample_info.rds")

# Number of samples within network
df <- data.frame(
  'Donors' = gsub("donor", "Donor ", donorNames),
  'All_samples' = sapply(samples_C, nrow),
  t(sapply(sample_info, function(t){ apply(t[, c(4,5)], 2, sum)}))
)
write.table(df, "number_of_samples.txt", row.names = FALSE, quote = FALSE, sep = "\t")