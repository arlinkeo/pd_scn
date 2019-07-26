# Expression heatmap of cell-type markers

library(ComplexHeatmap)
library(RColorBrewer)

# Normalize expression across samples
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))

# Networks to plots
nws <- networks #c("Network_C", "Network_D")

# Expression of cell-types in networks
expr1 <- lapply(donorNames, function(d){
  sapply(names(nws), function(n){
    # Mean expression of cell-types
    t <- t(sapply(markerlist, function(markers){
      e <- brainExprNorm[[d]][markers, sample_info[[d]][, n]==1]
      if(length(markers) > 1) apply(e, 2, mean) else e
    }))
    ontology_rows <- match(colnames(t), ontology$id)
    colnames(t) <- ontology[ontology_rows, c("acronym")]
    t <- t[, order(ontology[ontology_rows, c("graph_order")])]
    
    # Mean across groups of regions based on acronym
    t <- sapply(unique(colnames(t)), function(rg){
      cols <- which(colnames(t) == rg)
      if(length(cols) > 1) apply(t[, cols], 1, mean) else t[, cols]
    })
    rownames(t) <- gsub("_", " ", rownames(t))
    t
  }, simplify = FALSE)
})

# Expression averaged per celltype and network
expr2 <- lapply(names(expr1), function(d){
  t <- sapply(names(expr1[[d]]), function(n){
    e <- expr1[[d]][[n]]
    apply(e, 1, mean)
  })
  colnames(t) <- gsub("_", " ", colnames(t))
  t
})

# Expression averaged per celltype and network across donors
expr3 <- apply(simplify2array(expr2), c(1,2), mean)

#Heatmaps

# Column (sample) annotation
ahba_color <- paste0("#", ontology$color_hex_triplet)
names(ahba_color) <- ontology$acronym
col_color <- list(ahba_color = ahba_color)

# Heatmaps for each donor
heatmaps <- lapply(names(expr1), function(d){
  heatmaps <- lapply(names(expr1[[d]]), function(n){
    t <- expr1[[d]][[n]]
    ha <- HeatmapAnnotation(
      df = data.frame(ahba_color = colnames(t)), 
      col = col_color, 
      show_legend = FALSE,
      show_annotation_name = FALSE)
    Heatmap(t, name = 'Z-Score\nexpression',
            column_title = gsub("_", " ", n),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_names_side = c("top"), 
            row_names_side = c("left"),
            top_annotation = ha,
            width = unit(ncol(t)*0.8, "lines")
    )
  })
  heatmaps <- Reduce('+',heatmaps)
  heatmaps <- draw(heatmaps, column_title = gsub("donor", "Donor ", d), column_title_side = "top")
})
pdf("output/heatmap_celltypes_expanded.pdf", 50, 7)
heatmaps
dev.off()

heatmaps <- lapply(expr2, function(t){
  Heatmap(t, name = 'Z-Score\nexpression',
          column_title = gsub("donor", "Donor ", d),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          column_names_side = c("top"), 
          row_names_side = c("left")
  )
})
heatmaps <- Reduce('+',heatmaps)
heatmap_average <- Heatmap(expr3, name = 'Z-Score\nexpression',
          column_title = "Averaged",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          column_names_side = c("top"), 
          row_names_side = c("left")
)
heatmaps <- heatmap_average + heatmaps
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
tc_degs <- unique(unlist(tc_degs))

heatmaps <- lapply(donorNames, function(d){
  heatmaps <- lapply(names(expr1[[d]]), function(n){
    print(paste(d, "; ", n))
    
    # Expression of thalmuscholin DEGs
    t <-  brainExprNorm[[d]][tc_degs, sample_info[[d]][, n]==1]
    rownames(t) <- entrezId2Name(rownames(t))
    ontology_rows <- match(colnames(t), ontology$id)
    colnames(t) <- ontology[ontology_rows, c("acronym")]
    t <- t[, order(ontology[ontology_rows, c("graph_order")])]
   
    # Mean across groups of regions based on acronym
    t <- sapply(region_groups, function(rg){
      cols <- which(colnames(t) == rg)
      if(length(cols) > 1) apply(t[, cols], 1, mean) else t[, cols]
    })
    
    ha <- HeatmapAnnotation(
      df = data.frame(ahba_color = colnames(t)), 
      col = col_color, 
      show_legend = FALSE,
      show_annotation_name = FALSE)
    Heatmap(t, name = 'Z-Score\nexpression',
            column_title = gsub("_", " ", n),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_names_side = c("top"), 
            row_names_side = c("left"),
            top_annotation = ha,
            width = unit(ncol(t)*0.8, "lines")
    )
  })
  heatmaps <- Reduce('+',heatmaps)
  heatmaps <- draw(heatmaps, column_title = gsub("donor", "Donor ", d), column_title_side = "top")
})

pdf("output/heatmap_celltypes_thalamuscholin_degs.pdf", 60, 3.5)
heatmaps
dev.off()
