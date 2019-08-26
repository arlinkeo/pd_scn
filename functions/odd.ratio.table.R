# Odds ratios for two lists of gene sets
odd.ratio.table <- function(l1, l2){
  sapply(names(l1), function(n1){
    print(paste0(which(names(l1)==n1), ": ", n1))
    set <- l1[[n1]]
    # Overlap with each gene set
    sapply(l2, function(n2){
      odds.ratio(n2, set, length(ahba.genes())) # odds of l2
    })
  })
}
  