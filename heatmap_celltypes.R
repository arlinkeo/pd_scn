# Expression heatmap of cell-type markers

library(ComplexHeatmap)
library(RColorBrewer)

# Normalize expression across samples
brainExprNorm <- lapply(brainExpr, function(x) t(scale(t(x))))

#Heatmaps

# Row (marker) annotations
row_annot <- rep(names(markerlist), sapply(markerlist, length))#[-dup]
# row_annot <- data.frame(celltype = celltype)
# row_color <- rainbow(length(markerlist))
# names(row_color) <- names(markerlist)
# row_color <- row_color[row_annot]

# Column (sample) annotation
ahba_color <- ontology$color_hex_triplet
ahba_color <- paste0("#", ahba_color)
names(ahba_color) <- ontology$acronym
col_color <- list(ahba_color = ahba_color)

heatmaps <- lapply(donorNames, function(d){
  expr <- brainExprNorm[[d]]
  networks <- c("Network_C", "Network_D")
  heatmaps <- lapply(networks, function(n){
    
    i <- sample_info[[d]][, n]
    
    # Mean expression of cell-types
    t <- t(sapply(markerlist, function(markers){
      e <- expr[markers, i==1]
      if(length(markers) > 1) apply(e, 2, mean) else e
    }))
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

pdf("output/heatmap_celltypes.pdf", 18, 7)
heatmaps
dev.off()
