# Expression heatmap of cell-type markers

library(ComplexHeatmap)
library(circlize)
library(abind)

# Normalize expression across samples
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))

# Networks to plots
nws <- networks #c("Network_C", "Network_D")

########## Prepare expression data ##########

# Expression of cell-type markers in networks
expr1 <- sapply(names(nws), function(n){
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
heatmaps <- lapply(dimnames(expr4)[[2]], function(d){
  t <- expr4[, d,]
  Heatmap(t, name = 'Z-Score\nexpression', 
          col = col_fun,
          column_title = gsub("donor", "Donor ", d),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          column_names_side = c("top"), 
          row_names_side = c("left"),
          width = unit(ncol(t)*.8, "lines"), 
          height = unit(nrow(t)*.8, "lines")
  )
})
heatmaps <- Reduce('+',heatmaps)
pdf("output/heatmap_celltypes_simple.pdf", 13.3, 6)
heatmaps
dev.off()

# Heatmap ThalamusCholin DEGs
tc_degs <- lapply(degs, function(n){
  lapply(n, function(l){
    intersect(rownames(l), markerlist6$ThalamusCholin)
  })
})
sapply(tc_degs, function(n){
  sapply(n, function(l){
    paste(entrezId2Name(l), collapse = ", ")
  })
})
tc_degs <- tc_degs$Network_D$upregulated # select markers upregulated in D

heatmaps <- sapply(c("Network C", "Network D"), function(n){
  heatmaps <- lapply(donorNames, function(d){
    t <- expr2[[n]][[d]][tc_degs, ]
    rownames(t) <- entrezId2Name(rownames(t))
    ht_opt(heatmap_row_names_gp = gpar(fontface = "italic"))
    ha <- HeatmapAnnotation(
      df = data.frame(ahba_color = colnames(t)), 
      col = col_color, 
      show_legend = FALSE,
      show_annotation_name = FALSE)
    Heatmap(t, name = 'Z-Score\nexpression',
            col = col_fun,
            column_title = gsub("donor", "Donor ", d),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_names_side = c("top"), 
            row_names_side = c("left"),
            top_annotation = ha,
            width = unit(ncol(t)*.8, "lines"), 
            height = unit(nrow(t)*.8, "lines")
    )
  })
  ht_list <- Reduce('+', heatmaps)
  draw(ht_list, column_title = n, column_title_gp = gpar(fontface = "bold", fontsize = 18))
}, simplify = FALSE) 

pdf_size <- apply(sapply(heatmaps, function(h){
  s <- sapply(h@ht_list, function(d){
    p <- d@matrix_param
    c(width = p$width, height = p$height)
  })
  c(width = sum(s["width",])/5+3, height = s["height",1]/5+2)
}), 1, max)

pdf("output/heatmap_celltypes_thalamuscholin_degs.pdf", pdf_size[1], pdf_size[2])
heatmaps
dev.off()

acronyms <- unique(unlist(lapply(heatmaps, function(n){
  unlist(lapply(n@ht_list, function(d){
    colnames(d@matrix)
  }))
})))
acronyms <- ontology[ontology$acronym %in% acronyms, c("acronym", "name")]
acronyms$name <- gsub(", left|, right", "", acronyms$name)
acronyms <- unique(acronyms)
acronyms <- acronyms[-(which(duplicated(acronyms$acronym))-1),] # remove duplicates with typo's
acronyms <- unname(apply(acronyms, 1, function(x)paste(x, collapse = ": ")))
acronyms <- paste(acronyms, collapse = "; ")
acronyms
