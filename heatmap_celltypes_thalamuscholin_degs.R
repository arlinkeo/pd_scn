# Heatmap DEGs

# ThalamusCholin DEGs
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
  heatmaps <- lapply(names(networks), function(n){
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
    ht <- if (n == "Network_A") 
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