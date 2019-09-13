# Expression heatmap of cell-type markers

library(ComplexHeatmap)
library(circlize)
library(abind)

# Normalize expression across samples
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))

# Networks to plots
########## Prepare expression data ##########

# Expression of cell-type markers in networks
expr1 <- sapply(names(networks), function(n){
  l <- lapply(donorNames, function(d){
    brainExprNorm[[d]][unlist(markerlist), sample_info[[d]][, n]==1]
  })
}, simplify = FALSE)
names(expr1) <- gsub("_", " ", names(expr1))

# Mean across regions with same acronym within network
expr2 <- lapply(expr1, function(n){
  lapply(n, function(t){  
    ontology_rows <- match(colnames(t), ontology$id)
    colnames(t) <- ontology[ontology_rows, c("acronym")]
    t <- t[, order(ontology[ontology_rows, c("graph_order")])]
    sapply(unique(colnames(t)), function(rg){
      cols <- which(colnames(t) == rg)
      if(length(cols) > 1) apply(t[, cols], 1, mean) else t[, cols]
    })
  })
})

# Mean across cell-types
expr3 <- lapply(expr2, function(n){
  lapply(n, function(t){
    t <- t(sapply(markerlist, function(markers){
      e <- t[markers, ]
      if(length(markers) > 1) apply(e, 2, mean) else e
    }))
    rownames(t) <- gsub("_", " ", rownames(t))
    t
  })
})

# Mean within network
expr4 <- lapply(expr3, function(n){
  sapply(n, function(t){
    apply(t, 1, mean)
  })
})
expr4 <- simplify2array(expr4) # cell-types x donors x networks
expr4 <- abind(Average = apply(expr4, c(1,3), mean), expr4, along = 2) # Add mean across donors

########## Heatmaps ##########

# Column (sample) annotation
ahba_color <- paste0("#", ontology$color_hex_triplet)
names(ahba_color) <- ontology$acronym
col_color <- list(ahba_color = ahba_color)

# Same heat colors for all heatmaps
q1 <- max(sapply(brainExprNorm, function(x) quantile(abs(x), 0.9)))
col_fun <- colorRamp2(c(-q1, 0, q1), c("blue", "#EEEEEE", "red"))

# Simplified heatmap
t <- alply(expr4, 2, data.frame, .dims = TRUE)
names(t) <- gsub("donor", "Donor ", names(t))
split <- rep(names(t), sapply(t, ncol))
t <- t(Reduce(cbind, t))
rownames(t) <- gsub("\\.", " ", rownames(t))
hm <- Heatmap(t, name = 'Z-Score\nexpression', 
              row_split = split,
              col = col_fun,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 8),
              # row_names_side = c("left"),
              row_title_gp = gpar(fontsize = 10),
              row_title_rot = 0,
              column_names_side = c("top"),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 45,
              width = unit(ncol(t)*.8, "lines"), 
              height = unit(nrow(t)*.8, "lines")
)
pdf("output/heatmap_celltypes_simple.pdf", 7, 11.1)
hm
dev.off()

# Heatmap of thalamus cholinergic marker genes
heatmaps <- lapply(donorNames, function(d){
  hm <- lapply(names(expr2)[3:4], function(n){ # Only for network C & D
    t <- expr2[[n]][[d]][unique(unlist(ct_degs)), ]
    rownames(t) <- entrezId2Name(rownames(t))
    ht_opt(heatmap_row_names_gp = gpar(fontface = "italic"))
    ha <- HeatmapAnnotation(
      df = data.frame(ahba_color = colnames(t)), 
      col = col_color, 
      show_legend = FALSE,
      show_annotation_name = FALSE)
    Heatmap(t, name = 'Z-Score\nexpression',
            col = col_fun,
            column_title = n,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_names_side = c("top"), 
            row_names_side = c("left"),
            top_annotation = ha,
            width = unit(ncol(t)*.8, "lines"), 
            height = unit(nrow(t)*.8, "lines")
    )
  })
  ht_list <- Reduce('+', hm)
  draw(ht_list, column_title = gsub("donor", "Donor ", d), column_title_gp = gpar(fontface = "bold", fontsize = 18))
})

pdf_size <- apply(sapply(heatmaps, function(h){
  s <- sapply(h@ht_list, function(d){
    p <- d@matrix_param
    c(width = p$width, height = p$height)
  })
  c(width = sum(s["width",])/5+3, height = s["height",1]/5+2)
}), 1, max)

pdf("output/heatmap_celltypes_degs.pdf", pdf_size[1], pdf_size[2])
heatmaps
dev.off()

acronyms <- unique(unlist(lapply(heatmaps, function(d){
  unlist(lapply(d@ht_list, function(n){
    colnames(n@matrix)
  }))
})))
acronyms <- ontology[ontology$acronym %in% acronyms, c("acronym", "name")]
acronyms$name <- gsub(", left|, right", "", acronyms$name)
acronyms <- unique(acronyms)
acronyms <- acronyms[-(which(duplicated(acronyms$acronym))-1),] # remove duplicates with typo's
substr(colnames(acronyms), 1, 1) <- toupper(substr(colnames(acronyms), 1, 1))
substr(acronyms$Name, 1, 1) <- toupper(substr(acronyms$Name, 1, 1))
write.table(acronyms, file = "output/region_acronyms.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Thalamic samples within networks
thalamus_id <- ontology[ontology$name == "thalamus", "id"]
l <- sapply(c("Network_C", "Network_D"), function(n){
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

