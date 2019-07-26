# Expression heatmap of cell-type markers

library(ComplexHeatmap)
library(RColorBrewer)

# Normalize expression across samples
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))

#Heatmaps

# Row (marker) annotations
# row_annot <- rep(names(markerlist), sapply(markerlist, length))#[-dup]
# row_annot <- data.frame(celltype = celltype)
# row_color <- rainbow(length(markerlist))
# names(row_color) <- names(markerlist)
# row_color <- row_color[row_annot]

# Column (sample) annotation
ahba_color <- paste0("#", ontology$color_hex_triplet)
names(ahba_color) <- ontology$acronym
col_color <- list(ahba_color = ahba_color)

# Networks to plots
nws <- networks #c("Network_C", "Network_D")

# Heatmaps for each donor
heatmaps <- lapply(donorNames, function(d){
  expr <- as.matrix(brainExprNorm[[d]])

  heatmaps <- lapply(names(nws), function(n){
    # Mean expression of cell-types
    t <- t(sapply(markerlist, function(markers){
      e <- expr[markers, sample_info[[d]][, n]==1]
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

# Heatmaps simplified: mean expression of cell-types within network
expr <- lapply(donorNames, function(d){
  expr <- as.matrix(brainExprNorm[[d]])
  t <- sapply(names(nws), function(n){
    sapply(markerlist, function(markers){
      e <- expr[markers, sample_info[[d]][, n]==1]
      mean(e)
    })
  })
  colnames(t) <- gsub("_", " ", colnames(t))
  rownames(t) <- gsub("_", " ", rownames(t))
  t
})
t <- apply(simplify2array(expr), c(1,2), mean)
heatmap <- Heatmap(t, name = 'Z-Score\nexpression',
          column_title = "Averaged",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          column_names_side = c("top"), 
          row_names_side = c("left"),
          width = unit(ncol(t)*0.8, "lines")
)
heatmaps <- lapply(expr, function(t){
  Heatmap(t, name = 'Z-Score\nexpression',
          column_title = gsub("donor", "Donor ", d),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          column_names_side = c("top"), 
          row_names_side = c("left"),
          width = unit(ncol(t)*0.8, "lines")
  )
})
heatmaps <- Reduce('+',heatmaps)
heatmaps <- heatmap + heatmaps

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
  expr <- brainExprNorm[[d]]
  heatmaps <- lapply(names(networks)[3:4], function(n){
    print(paste(d, "; ", n))
    i <- sample_info[[d]][, n]
    
    # Expression of thalmuscholin DEGs
    t <- expr[tc_degs, i == 1]
    rownames(t) <- entrezId2Name(rownames(t))
    col_order <- order(ontology[match(colnames(t), ontology$id), c("graph_order")])
    t <- t[, col_order]
    colnames(t) <- ontology[match(colnames(t), ontology$id), c("acronym")]
    region_groups <- unique(colnames(t))
    
    # Mean across groups of regions based on acronym
    t <- sapply(region_groups, function(rg){
      cols <- which(colnames(t) == rg)
      if(length(cols) > 1) apply(t[, cols], 1, mean) else t[, cols]
      # 
    })
    
    df <- data.frame(ahba_color = colnames(t))
    ha <- HeatmapAnnotation(df, col = col_color, show_legend = FALSE)
    ht <- if (n == "Network_C") 
      Heatmap(t, name = 'Z-Score\nexpression',
              column_title = gsub("_", " ", n),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_names_side = c("top"), 
              row_names_side = c("left"),
              clustering_method_columns = "single",
              top_annotation = ha,
              width = unit(ncol(t)*0.8, "lines")
      ) 
    else
      Heatmap(t, name = 'Z-Score\nexpression',
              column_title = gsub("_", " ", n),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_names_side = c("top"), 
              row_names_side = c("left"),
              clustering_method_columns = "single",
              top_annotation = ha,
              width = unit(ncol(t), "lines"),
              show_row_names = FALSE
      )
    
  })
  heatmaps <- Reduce('+',heatmaps)
  heatmaps <- draw(heatmaps, column_title = gsub("donor", "Donor ", d), column_title_side = "top")
})

pdf("output/heatmap_celltypes_thalamuscholin_degs.pdf", 60, 3.5)
heatmaps
dev.off()
