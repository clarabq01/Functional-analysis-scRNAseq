#### FUNCTIONAL ANALYSIS #####
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratData)
library(scDblFinder)
library(presto)
library(SingleCellExperiment)
library(dplyr)
library(scran)
library(patchwork)
library(viridis)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ggforce)
library(gghalves)
library(ggridges)
library(RCurl)
library(glmGamPoi)
# library(BPCells) # BPCells is an R package that allows for computationally efficient single-cell analysis.BPCells allows us to easily analyze these large datasets in memory (ERROR)
library(AnnotationHub)
library(ensembldb)
library(reticulate)
library(sctransform)
library(SingleR)
library(scales)
library(cowplot)

10. FUNCTIONAL ANALYSIS #############NOT FINISHED
{
  10.1 DIFFERENTIAL EXPRESSION IN CLUSTERS THAT VARY ACCORDING TO SAMPLE
  10.1.1 SEURAT WORKFLOW
  {
    # Pseudobulk the counts based on donor-condition-celltype 
    pseudo_H <- AggregateExpression(SCT_objH, assays = "SCT", return.seurat = T, group.by = c("sample","dataset", "combined_label"))
    View(pseudo_H@meta.data)
    # Each 'cell' is a donor-condition-celltype pseudobulk profile
    tail(Cells(pseudo_H))
    pseudo_H$sample.label <- paste(pseudo_H$sample, pseudo_H$combined_label, sep = "_")
    Idents(pseudo_H) <- "sample.label" #we can identify cells according to sample, dataset, region, etc. and each cell is detected individually
    View(pseudo_H@meta.data)
    library(DESeq2)
    bulkDEGsample <- FindMarkers(object = pseudo_H, 
                                 ident.1 = "ctrl_Neuron", 
                                 ident.2 = "inj_Neuron",
                                 test.use = "DESeq2")
    head(bulkDEGsample, n = 5)
    #               p_val  avg_log2FC pct.1 pct.2    p_val_adj
    # Elane  9.213295e-09  -5.977280 0.000 0.333 0.0002311616
    # Lgals3 1.308263e-06  -6.714246 0.667 1.000 0.0328243161
    #  v  5.829756e-06  -6.622052 0.333 1.000 0.1462685814
    # Stab1  7.458110e-06  -6.209453 0.333 1.000 0.1871239742
    # Ifi30  2.030241e-05  -5.262095 0.667 1.000 0.5093874684
    
    top10Neuron <- bulkDEGsample %>% dplyr::slice_min(get(grep("^p_val_adj", colnames(bulkDEGsample), value = TRUE)), n = 10)
    
    top10Neuron
    #                      p_val  avg_log2FC pct.1 pct.2    p_val_adj
    # Elane         9.213295e-09 -5.9772799 0.000 0.333 0.0002311616
    # Lgals3        1.308263e-06 -6.7142455 0.667 1.000 0.0328243161
    # C3ar1         5.829756e-06 -6.6220518 0.333 1.000 0.1462685814
    # Stab1         7.458110e-06 -6.2094534 0.333 1.000 0.1871239742
    # Ifi30         2.030241e-05 -5.2620948 0.667 1.000 0.5093874684
    # Ncf1          2.059499e-05 -5.3923174 0.333 1.000 0.5167283123
    # Apobec1       3.479760e-05 -5.9772799 0.333 1.000 0.8730717520
    # Cd84          3.837705e-05 -5.4195389 0.667 1.000 0.9628802280
    # Cybb          4.850501e-05 -5.2352165 0.333 1.000 1.0000000000
    
    
    
    bulkDEGsample2 <- FindMarkers(object = pseudo_H, 
                                  ident.1 = "inj_Neuron", 
                                  ident.2 = "ctrl_Neuron",
                                  test.use = "DESeq2")
    head(bulkDEGsample2, n = 5)
    #               p_val avg_log2FC pct.1 pct.2    p_val_adj
    # Elane  2.945808e-09   5.977280 0.333 0.000 7.391032e-05
    # Lgals3 1.308315e-06   6.714246 1.000 0.667 3.282562e-02
    # C3ar1  5.829952e-06   6.622052 1.000 0.333 1.462735e-01
    # Stab1  7.458345e-06   6.209453 1.000 0.333 1.871299e-01
    # Ifi30  2.030304e-05   5.262095 1.000 0.667 5.094033e-01
    
    top10Neuroninj <- bulkDEGsample2 %>% dplyr::slice_min(get(grep("^avg_log2FC", colnames(bulkDEGsample2), value = TRUE)), n = 10)
    
    top10Neuroninj
    # p_val avg_log2FC pct.1 pct.2    p_val_adj
    # Elane         2.945808e-09  5.9772799 0.333 0.000 7.391032e-05
    # Lgals3        1.308315e-06  6.7142455 1.000 0.667 3.282562e-02
    # C3ar1         5.829952e-06  6.6220518 1.000 0.333 1.462735e-01
    # Stab1         7.458345e-06  6.2094534 1.000 0.333 1.871299e-01
    # Ifi30         2.030304e-05  5.2620948 1.000 0.667 5.094033e-01
    # Ncf1          2.059559e-05  5.3923174 1.000 0.333 5.167434e-01
    # Apobec1       3.479875e-05  5.9772799 1.000 0.333 8.731006e-01
    # Cd84          3.837812e-05  5.4195389 1.000 0.667 9.629071e-01
    # Cybb          4.850638e-05  5.2352165 1.000 0.333 1.000000e+00
    
    
    VlnPlot(pseudo_H, features <- c('Elane','Lgals3'), idents = c("ctrl_Neuron", "inj_Neuron")) 
    
    # Create a basic volcano plot
    ggplot(data = bulkDEGsample, aes(x = avg_log2FC, y = -log10(p_val_adj)), label= TRUE) + geom_point()
    
    FeaturePlot(SCT_objH, 
                reduction = "umap.harmony", split.by = c("sample"),
                features = 'Ppbp', 
                order = TRUE,
                min.cutoff = 'q10', 
                label = TRUE)
    
  }
  10.1.2 BEST PRACTICES
  # Two tools for DE analysis: edgeR with a quasi-likelihood test and MAST with random effects.
  EDGER ALL MARKERS ADDING DATASET
  # edgeR is based on a negative binomial model of gene expression and uses a generalized linear model (GLM) framework, 
  # the enables us to include other factors such as batch to the model.
  {
  library(edgeR)
  # We will need to work with raw counts 
  View(SCT_objH_labels@meta.data)
  
  SCT_objH_labels
  # An object of class Seurat 
  # 53205 features across 56158 samples within 2 assays 
  # Active assay: SCT (25090 features, 3000 variable features)
  # 3 layers present: counts, data, scale.data
  # 1 other assay present: RNA
  # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
  
  DefaultAssay(SCT_objH_labels) <- "RNA"
  SCT_objH_labels[["SCT"]] <- NULL
  SCT_objH_labels
  # An object of class Seurat 
  # 28115 features across 56158 samples within 1 assay 
  # Active assay: RNA (28115 features, 0 variable features)
  # 1 layer present: counts
  # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
  
  
  ### Create pseudo-bulk samples: Seurat2PB
  # Pseudo-bulk samples are created by aggregating read counts together 
  # for all the cells with the same combination
  y <- Seurat2PB(SCT_objH_labels, sample="sample.dataset", cluster="combined_label") #for creating a pseudo-bulk DGEList object from a Seurat object.
  dim(y)
  # 28115    80
  
  counts <- y$counts
  samples <- y$samples
  
  group <- factor(c('ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',  'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj'))
  y$samples$group <- group
  samples <- y$samples
  
  ### Filtering and normalization
  # filter out library size below the threshold of 50,000
  summary(y$samples$lib.size)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # 158    61765   192546  1872731   719799 49993207 
  keep.samples <- y$samples$lib.size > 5e4 #not needed if libraries are bigger than this but run anyways
  keep.samples
  table(keep.samples)
  y <- y[, keep.samples]
  
  # filter out lowly expressed genes
  keep.genes <- filterByExpr(y, group=y$samples$cluster,
                             min.count=10, min.total.count=20)
  table(keep.genes)
  # FALSE  TRUE 
  # 14817 13298 # eliminamos esos genes FALSE y nos quedamos solo con los TRUE
  y <- y[keep.genes, , keep=FALSE]
  
  ### TMM normalization is performed to estimate effective library sizes.
  y <- normLibSizes(y)
  samples <- y$samples
  
  ### Data exploration: explore the pseudo-bulk profiles
  # using a multi-dimensional scaling (MDS) plot
  cluster <- as.factor(y$samples$cluster)
  plotMDS(y, pch=16, col=c(1:14)[cluster], main="MDS")
  legend("top", legend=paste0("cluster", levels(cluster)),
         pch=16, col=1:14, y.intersp = 0.6)
  
  group2 <- as.factor(y$samples$group)
  plotMDS(y, pch=16, col=c(1:2)[group2], main="MDS")
  legend("top", legend=paste0("group2", levels(group2)),
         pch=16, col=1:14, y.intersp = 0.6)
  
  
  ### Dispersion estimation with a specified design
  # Design a matrix
  type <- factor(y$samples$group)
  design <- model.matrix(~cluster + type)
  
  # estimate dispersion using the estimateDisp function, and the GLM method
  y <- estimateDisp(y, design, robust = TRUE) # common, trend, tagwise
  plotBCV(y) # visualize with plotBCV
  
  # estimate quasi-likelihood (QL) dispersion using the glmQLFit function
  fit <- glmQLFit(y, design, robust = TRUE) # Generalized linear models
  plotQLDisp(fit) # visualize with plotQLDisp
  
  ### Marker genes identification
  ncls <- nlevels(cluster)
  contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                  + matrix(0, ncol(design)-ncls, ncls) )
  diag(contr) <- 1
  contr[1,] <- 0
  rownames(contr) <- colnames(design)
  colnames(contr) <- paste0("cluster", levels(cluster))
  
  qlf <- list()
  for(i in 1:ncls){
    qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
    qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
  }
  
  Astro <- qlf[[1]][["table"]]
  Bcells <- qlf[[2]][["table"]]
  Ependymal <- qlf[[3]][["table"]]
  Erythrocytes <- qlf[[4]][["table"]]
  Granulocytes <- qlf[[5]][["table"]]
  Macrophages <- qlf[[6]][["table"]]
  Microglia <- qlf[[7]][["table"]]
  Monocytes <- qlf[[8]][["table"]]
  Neuron <- qlf[[9]][["table"]]
  NKcells <- qlf[[10]][["table"]]
  Oligodendrocytes <- qlf[[11]][["table"]]
  OPC <- qlf[[12]][["table"]]
  Tcells <- qlf[[13]][["table"]]
  Vacular <- qlf[[14]][["table"]]
  
  # show the numbers of DE genes under each comparison
  dt <- lapply(lapply(qlf, decideTestsDGE), summary)
  dt.all <- do.call("cbind", dt)
  dt.all
  
  #                     clusterAstrocyte_vs_others clusterB cells_vs_others clusterEpendymal_vs_others clusterErythrocytes_vs_others
  # Down                         1060                      413                       1205                             9
  # NotSig                       8321                    11738                       9005                         12328
  # Up                           3917                     1147                       3088                           961
  #                       clusterGranulocytes_vs_others clusterMacrophages_vs_others clusterMicroglia_vs_others clusterMonocytes_vs_others
  # Down                            1510                          193                        417                        303
  # NotSig                          9574                        11565                      10397                      11347
  # Up                              2214                         1540                       2484                       1648
  #                     clusterNeuron_vs_others clusterNK cells_vs_others clusterOligodendrocytes_vs_others clusterOPC_vs_others
  # Down                       346                       301                               768                 1040
  # NotSig                    8273                     12290                              8605                 9298
  # Up                        4679                       707                              3925                 2960
  #                      clusterT cells_vs_others clusterVascular_vs_others
  # Down                        308                       356
  # NotSig                    12170                     10799
  # Up                          820                      2143
  
  # select the top 20 marker (up-regulated) genes for each cluster.
  top <- 20
  topMarkers <- list()
  for(i in 1:ncls) {
    ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
    up <- qlf[[i]]$table$logFC > 0
    topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
  }
  topMarkers <- unique(unlist(topMarkers))
  topMarkers
  lcpm <- cpm(y, log=TRUE)
  topMarkers20 <- lcpm[topMarkers, ]
  
  # Heatmap
  # Create a color palette with enough colors for all clusters
  cluster_colors <- turbo(length(levels(cluster)))
  
  # Create the annotation dataframe with cluster labels
  annot <- data.frame(cluster = factor(cluster))
  rownames(annot) <- colnames(y)
  
  # Assign unique colors to each cluster label
  ann_colors <- list(cluster = cluster_colors)
  names(ann_colors$cluster) <- levels(cluster)
  
  # Plot with unique colors for each cluster
  pheatmap::pheatmap(topMarkers20, 
                     breaks = seq(-2, 2, length.out = 101),
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     scale = "row", 
                     cluster_cols = TRUE, 
                     border_color = "NA", 
                     fontsize_row = 5,
                     treeheight_row = 70, 
                     treeheight_col = 70, 
                     cutree_cols = 7,
                     clustering_method = "ward.D2", 
                     show_colnames = FALSE,
                     annotation_col = annot, 
                     annotation_colors = ann_colors)
  
  }
  EDGER ALL MARKERS NOT ADDING DATASET
  {
    library(edgeR)
    # We will need to work with raw counts 
    View(SCT_objH_labels@meta.data)
    
    SCT_objH_labels
    # An object of class Seurat 
    # 53205 features across 56158 samples within 2 assays 
    # Active assay: SCT (25090 features, 3000 variable features)
    # 3 layers present: counts, data, scale.data
    # 1 other assay present: RNA
    # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
    
    DefaultAssay(SCT_objH_labels) <- "RNA"
    SCT_objH_labels[["SCT"]] <- NULL
    SCT_objH_labels
    # An object of class Seurat 
    # 28115 features across 56158 samples within 1 assay 
    # Active assay: RNA (28115 features, 0 variable features)
    # 1 layer present: counts
    # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
    
    
    ### Create pseudo-bulk samples: Seurat2PB
    # Pseudo-bulk samples are created by aggregating read counts together 
    # for all the cells with the same combination
    y <- Seurat2PB(SCT_objH_labels, sample="sample", cluster="combined_label") #for creating a pseudo-bulk DGEList object from a Seurat object.
    dim(y)
    # 28115    28
    
    counts <- y$counts
    samples <- y$samples
    
    group <- factor(c('ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',  'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj'))
    y$samples$group <- group
    samples <- y$samples
    
    ### Filtering and normalization
    # filter out library size below the threshold of 50,000
    summary(y$samples$lib.size)
    # Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
    # 98424   320117   786996  5350659  3724328 63581689  
    keep.samples <- y$samples$lib.size > 5e4 #not needed if libraries are bigger than this but run anyways
    keep.samples
    table(keep.samples)
    y <- y[, keep.samples]
    
    # filter out lowly expressed genes
    keep.genes <- filterByExpr(y, group=y$samples$cluster,
                               min.count=10, min.total.count=20)
    table(keep.genes)
    # FALSE  TRUE 
    # 12816 15299  # eliminamos esos genes FALSE y nos quedamos solo con los TRUE
    y <- y[keep.genes, , keep=FALSE]
    
    ### TMM normalization is performed to estimate effective library sizes.
    y <- normLibSizes(y)
    samples <- y$samples
    
    ### Data exploration: explore the pseudo-bulk profiles
    # using a multi-dimensional scaling (MDS) plot
    # Generate colors using turbo palette
    cluster_colors <- viridis(length(unique(cluster)))
    # Plot MDS with viridis colors
    cluster <- as.factor(y$samples$cluster)
    plotMDS(y, pch=16, col=cluster_colors[cluster], main="MDS")
    legend("top", legend=paste0("cluster", levels(cluster)),
           pch=16, col=cluster_colors, y.intersp = 0.6)
    
    group2 <- as.factor(y$samples$group)
    plotMDS(y, pch=16, col=c(1:2)[group2], main="MDS")
    legend("top", legend=paste0("group2", levels(group2)),
           pch=16, col=1:14, y.intersp = 0.6)
    
    
    ### Dispersion estimation with a specified design
    # Design a matrix
    # type <- factor(y$samples$group)
    # design <- model.matrix(~cluster + type)
    design <- model.matrix(~cluster)
    
    # estimate dispersion using the estimateDisp function, and the GLM method
    y <- estimateDisp(y, design, robust = TRUE) # common, trend, tagwise
    plotBCV(y) # visualize with plotBCV
    
    # estimate quasi-likelihood (QL) dispersion using the glmQLFit function
    fit <- glmQLFit(y, design, robust = TRUE) # Generalized linear models
    plotQLDisp(fit) # visualize with plotQLDisp
    
    ### Marker genes identification
    ncls <- nlevels(cluster)
    contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                    + matrix(0, ncol(design)-ncls, ncls) )
    diag(contr) <- 1
    contr[1,] <- 0
    rownames(contr) <- colnames(design)
    colnames(contr) <- paste0("cluster", levels(cluster))
    
    qlf <- list()
    for(i in 1:ncls){
      qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
      qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
    }
    
    Astro <- qlf[[1]][["table"]]
    Bcells <- qlf[[2]][["table"]]
    Ependymal <- qlf[[3]][["table"]]
    Erythrocytes <- qlf[[4]][["table"]]
    Granulocytes <- qlf[[5]][["table"]]
    Macrophages <- qlf[[6]][["table"]]
    Microglia <- qlf[[7]][["table"]]
    Monocytes <- qlf[[8]][["table"]]
    Neuron <- qlf[[9]][["table"]]
    NKcells <- qlf[[10]][["table"]]
    Oligodendrocytes <- qlf[[11]][["table"]]
    OPC <- qlf[[12]][["table"]]
    Tcells <- qlf[[13]][["table"]]
    Vacular <- qlf[[14]][["table"]]
    
    # show the numbers of DE genes under each comparison
    dt <- lapply(lapply(qlf, decideTestsDGE), summary)
    dt.all <- do.call("cbind", dt)
    dt.all
    
    #                   clusterAstrocyte_vs_others clusterB cells_vs_others clusterEpendymal_vs_others clusterErythrocytes_vs_others
    # Down                         1034                      327                       1405                            12
    # NotSig                      10688                    13649                      10624                         13838
    # Up                           3577                     1323                       3270                          1449
    #                   clusterGranulocytes_vs_others clusterMacrophages_vs_others clusterMicroglia_vs_others clusterMonocytes_vs_others
    # Down                            1391                          183                        329                        260
    # NotSig                         11307                        12381                      12011                      12827
    # Up                              2601                         2735                       2959                       2212
    #                   clusterNeuron_vs_others clusterNK cells_vs_others clusterOligodendrocytes_vs_others clusterOPC_vs_others
    # Down                      1914                       299                              1119                 1267
    # NotSig                    8680                     13964                             10691                11304
    # Up                        4705                      1036                              3489                 2728
    #                    clusterT cells_vs_others clusterVascular_vs_others
    # Down                        301                       343
    # NotSig                    13770                     12655
    # Up                         1228                      2301
    # 
    # select the top 20 marker (up-regulated) genes for each cluster.
    top <- 20
    topMarkers <- list()
    for(i in 1:ncls) {
      ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
      up <- qlf[[i]]$table$logFC > 0
      topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
    }
    topMarkers <- unique(unlist(topMarkers))
    topMarkers
    lcpm <- cpm(y, log=TRUE)
    topMarkers20 <- lcpm[topMarkers, ]
    
    # Heatmap
    # Create a color palette with enough colors for all clusters
    cluster_colors <- turbo(length(levels(cluster)))
    
    # Create the annotation dataframe with cluster labels
    annot <- data.frame(cluster = factor(cluster))
    rownames(annot) <- colnames(y)
    
    # Assign unique colors to each cluster label
    ann_colors <- list(cluster = cluster_colors)
    names(ann_colors$cluster) <- levels(cluster)
    
    # Plot with unique colors for each cluster
    pheatmap::pheatmap(topMarkers20, 
                       breaks = seq(-2, 2, length.out = 101),
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       scale = "row", 
                       cluster_cols = TRUE, 
                       border_color = "NA", 
                       fontsize_row = 5,
                       treeheight_row = 70, 
                       treeheight_col = 70, 
                       cutree_cols = 7,
                       clustering_method = "ward.D2", 
                       show_colnames = FALSE,
                       annotation_col = annot, 
                       annotation_colors = ann_colors)
    
  }
  EDGER IN SUBSETS #ERROR 
  {
  neuron_subset <- subset(SCT_objH_labels, subset = combined_label == "Neuron")
  View(neuron_subset@meta.data)
  neuron_subset
  
  ### Create pseudo-bulk samples: Seurat2PB
  # Pseudo-bulk samples are created by aggregating read counts together 
  # for all the cells with the same combination
  y <- Seurat2PB(neuron_subset, sample="dataset", cluster="sample") #for creating a pseudo-bulk DGEList object from a Seurat object.
  dim(y)
  # 28115    6
  
  counts <- y$counts
  samples <- y$samples
  
  group <- factor(c('ctrl', 'inj', 'ctrl', 'inj', 'ctrl', 'inj'))
  y$samples$group <- group
  samples <- y$samples
  
  ### Filtering and normalization
  # filter out library size below the threshold of 50,000
  summary(y$samples$lib.size)
  # Min.  1st Qu.   Median     Mean   3rd Qu.     Max. 
  #  98725  293336  820510 2126050 3247330 6842676  
  keep.samples <- y$samples$lib.size > 5e4 #not needed if libraries are bigger than this but run anyways
  keep.samples
  table(keep.samples)
  y <- y[, keep.samples]
  
  # filter out lowly expressed genes
  keep.genes <- filterByExpr(y, group=y$samples$cluster,
                             min.count=10, min.total.count=20)
  table(keep.genes)
  # FALSE  TRUE 
  # 18418  9697  # eliminamos esos genes FALSE y nos quedamos solo con los TRUE
  y <- y[keep.genes, , keep=FALSE]
  
  ### TMM normalization is performed to estimate effective library sizes.
  y <- normLibSizes(y)
  samples <- y$samples
  
  ### Data exploration: explore the pseudo-bulk profiles
  # using a multi-dimensional scaling (MDS) plot
  # Plot MDS with viridis colors
  cluster <- as.factor(y$samples$cluster)
  # Generate colors using turbo palette
  cluster_colors <- viridis(length(unique(cluster)))
  plotMDS(y, pch=16, col=cluster_colors[cluster], main="MDS")
  legend("top", legend=paste0("cluster", levels(cluster)),
         pch=16, col=cluster_colors, y.intersp = 0.6)
  
  group2 <- as.factor(y$samples$group)
  plotMDS(y, pch=16, col=c(1:2)[group2], main="MDS")
  legend("top", legend=paste0("group2", levels(group2)),
         pch=16, col=1:14, y.intersp = 0.6)
  
  
  ### Dispersion estimation with a specified design
  # Design a matrix
  # type <- factor(y$samples$group)
  # design <- model.matrix(~cluster + type)
  design <- model.matrix(~cluster)
  
  # estimate dispersion using the estimateDisp function, and the GLM method
  y <- estimateDisp(y, design, robust = TRUE) # common, trend, tagwise
  plotBCV(y) # visualize with plotBCV
  
  # estimate quasi-likelihood (QL) dispersion using the glmQLFit function
  fit <- glmQLFit(y, design, robust = TRUE) # Generalized linear models
  plotQLDisp(fit) # visualize with plotQLDisp
  
  ### Marker genes identification
  ncls <- nlevels(cluster)
  contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                  + matrix(0, ncol(design)-ncls, ncls) )
  diag(contr) <- 1
  contr[1,] <- 0
  rownames(contr) <- colnames(design)
  colnames(contr) <- paste0("cluster", levels(cluster))
  
  qlf <- list()
  for(i in 1:ncls){
    qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
    qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
  }
  
  CtrlNeuron <- qlf[[1]][["table"]]
  InjNeuron <- qlf[[2]][["table"]]
  
  # show the numbers of DE genes under each comparison
  dt <- lapply(lapply(qlf, decideTestsDGE), summary)
  dt.all <- do.call("cbind", dt)
  dt.all
  
  # clusterctrl_vs_others clusterinj_vs_others
  # Down                       0                    0
  # NotSig                  9697                 9697
  # Up                         0                    0
  
  # select the top 20 marker (up-regulated) genes for each cluster.
  top <- 20
  topMarkers <- list()
  for(i in 1:ncls) {
    ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
    up <- qlf[[i]]$table$logFC > 0
    topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
  }
  topMarkers <- unique(unlist(topMarkers))
  topMarkers
  lcpm <- cpm(y, log=TRUE)
  topMarkers20 <- lcpm[topMarkers, ]
  
  # Heatmap
  # Create a color palette with enough colors for all clusters
  # Create the annotation dataframe with cluster labels
  annot <- data.frame(cluster = factor(cluster))
  rownames(annot) <- colnames(y)
  
  # Assign unique colors to each cluster label
  ann_colors <- list(cluster = cluster_colors)
  names(ann_colors$cluster) <- levels(cluster)
  
  # Plot with unique colors for each cluster
  pheatmap::pheatmap(topMarkers20, 
                     breaks = seq(-2, 2, length.out = 101),
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     scale = "row", 
                     cluster_cols = TRUE, 
                     border_color = "NA", 
                     fontsize_row = 5,
                     treeheight_row = 70, 
                     treeheight_col = 70, 
                     cutree_cols = 7,
                     clustering_method = "ward.D2", 
                     show_colnames = FALSE,
                     annotation_col = annot, 
                     annotation_colors = ann_colors)
  
  }
  EDGER GENERATING A COLUMN TO CONCATENATES LABEL AND CONDITION 
  {
    library(edgeR)
    # We will need to work with raw counts 
    View(SCT_objH_labels@meta.data)
    
    SCT_objH_labels
    # An object of class Seurat 
    # 53205 features across 56158 samples within 2 assays 
    # Active assay: SCT (25090 features, 3000 variable features)
    # 3 layers present: counts, data, scale.data
    # 1 other assay present: RNA
    # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
    
    DefaultAssay(SCT_objH_labels) <- "RNA"
    SCT_objH_labels[["SCT"]] <- NULL
    SCT_objH_labels
    # An object of class Seurat 
    # 28115 features across 56158 samples within 1 assay 
    # Active assay: RNA (28115 features, 0 variable features)
    # 1 layer present: counts
    # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
    
    SCT_objH_labels$label.condition <- paste(SCT_objH_labels$sample, SCT_objH_labels$combined_label, sep = "_")
    ### Create pseudo-bulk samples: Seurat2PB
    # Pseudo-bulk samples are created by aggregating read counts together 
    # for all the cells with the same combination
    y <- Seurat2PB(SCT_objH_labels, sample="dataset", cluster="label.condition") #for creating a pseudo-bulk DGEList object from a Seurat object.
    dim(y)
    # 28115    80
    
    counts <- y$counts
    samples <- y$samples
    
    group <- factor(c('ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',
                      'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj',
                      'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',
                      'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj',
                      'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',
                      'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj'))
    y$samples$group <- group
    samples <- y$samples
    
    ### Filtering and normalization
    # filter out library size below the threshold of 50,000
    summary(y$samples$lib.size)
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # 158    61765   192546  1872731   719799 49993207 
    keep.samples <- y$samples$lib.size > 5e4 #not needed if libraries are bigger than this but run anyways
    keep.samples
    table(keep.samples)
    y <- y[, keep.samples]
    
    # filter out lowly expressed genes
    keep.genes <- filterByExpr(y, group=y$samples$cluster,
                               min.count=10, min.total.count=20)
    table(keep.genes)
    # FALSE  TRUE 
    # 12592 15523  # eliminamos esos genes FALSE y nos quedamos solo con los TRUE
    y <- y[keep.genes, , keep=FALSE]
    
    ### TMM normalization is performed to estimate effective library sizes.
    y <- normLibSizes(y)
    samples <- y$samples
    
    ### Data exploration: explore the pseudo-bulk profiles
    # using a multi-dimensional scaling (MDS) plot
    cluster <- as.factor(y$samples$cluster)
    # Generate colors using turbo palette
    cluster_colors <- turbo(length(unique(cluster)))
    plotMDS(y, pch=16, col=cluster_colors[cluster], main="MDS")
    legend("top", legend=paste0("cluster", levels(cluster)),
           pch=16, col=cluster_colors, y.intersp = 0.6)
    
    group2 <- as.factor(y$samples$group)
    plotMDS(y, pch=16, col=c(1:2)[group2], main="MDS")
    legend("top", legend=paste0("group2", levels(group2)),
           pch=16, col=1:14, y.intersp = 0.6)
    
    
    ### Dispersion estimation with a specified design
    # Design a matrix
    type <- factor(y$samples$group)
    design <- model.matrix(~cluster)
    
    # estimate dispersion using the estimateDisp function, and the GLM method
    y <- estimateDisp(y, design, robust = TRUE) # common, trend, tagwise
    plotBCV(y) # visualize with plotBCV
    
    # estimate quasi-likelihood (QL) dispersion using the glmQLFit function
    fit <- glmQLFit(y, design, robust = TRUE) # Generalized linear models
    plotQLDisp(fit) # visualize with plotQLDisp
    
    ### Marker genes identification
    ncls <- nlevels(cluster)
    contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                    + matrix(0, ncol(design)-ncls, ncls) )
    diag(contr) <- 1
    contr[1,] <- 0
    rownames(contr) <- colnames(design)
    colnames(contr) <- paste0("cluster", levels(cluster))
    
    qlf <- list()
    for(i in 1:ncls){
      qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
      qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
    }
    
    Astro_ctrl <- qlf[[1]][["table"]]
    Bcells_ctrl <- qlf[[2]][["table"]]
    Ependymal_ctrl <- qlf[[3]][["table"]]
    Erythrocytes_ctrl <- qlf[[4]][["table"]]
    Granulocytes_ctrl <- qlf[[5]][["table"]]
    Macrophages_ctrl <- qlf[[6]][["table"]]
    Microglia_ctrl <- qlf[[7]][["table"]]
    Monocytes_ctrl <- qlf[[8]][["table"]]
    Neuron_ctrl <- qlf[[9]][["table"]]
    NKcells_ctrl <- qlf[[10]][["table"]]
    Oligodendrocytes_ctrl <- qlf[[11]][["table"]]
    OPC_ctrl <- qlf[[12]][["table"]]
    Tcells_ctrl <- qlf[[13]][["table"]]
    Vacular_ctrl <- qlf[[14]][["table"]]
    Astro_inj <- qlf[[15]][["table"]]
    Bcells_inj <- qlf[[16]][["table"]]
    Ependymal_inj <- qlf[[17]][["table"]]
    Erythrocytes_inj <- qlf[[18]][["table"]]
    Granulocytes_inj <- qlf[[19]][["table"]]
    Macrophages_inj <- qlf[[20]][["table"]]
    Microglia_inj <- qlf[[21]][["table"]]
    Monocytes_inj <- qlf[[22]][["table"]]
    Neuron_inj <- qlf[[23]][["table"]]
    NKcells_inj <- qlf[[24]][["table"]]
    Oligodendrocytes_inj <- qlf[[25]][["table"]]
    OPC_inj <- qlf[[26]][["table"]]
    Tcells_inj <- qlf[[27]][["table"]]
    Vacular_inj <- qlf[[28]][["table"]]
    
    
    # show the numbers of DE genes under each comparison
    dt <- lapply(lapply(qlf, decideTestsDGE), summary)
    dt.all <- do.call("cbind", dt)
    dt.all
    
    # clusterctrl_Astrocyte_vs_others clusterctrl_B cells_vs_others clusterctrl_Ependymal_vs_others
    # Down                               325                            57                             268
    # NotSig                           11629                         14646                           12772
    # Up                                3569                           820                            2483
    # clusterctrl_Erythrocytes_vs_others clusterctrl_Granulocytes_vs_others clusterctrl_Macrophages_vs_others
    # Down                                    9                                314                                 1
    # NotSig                              15068                              13796                             14772
    # Up                                    446                               1413                               750
    # clusterctrl_Microglia_vs_others clusterctrl_Monocytes_vs_others clusterctrl_Neuron_vs_others
    # Down                                37                              28                          111
    # NotSig                           12783                           14307                        11187
    # Up                                2703                            1188                         4225
    # clusterctrl_NK cells_vs_others clusterctrl_Oligodendrocytes_vs_others clusterctrl_OPC_vs_others
    # Down                               11                                     74                       137
    # NotSig                          15073                                  12528                     13916
    # Up                                439                                   2921                      1470
    # clusterctrl_T cells_vs_others clusterctrl_Vascular_vs_others clusterinj_Astrocyte_vs_others
    # Down                               7                             66                            134
    # NotSig                         15006                          13583                          12714
    # Up                               510                           1874                           2675
    # clusterinj_B cells_vs_others clusterinj_Ependymal_vs_others clusterinj_Erythrocytes_vs_others
    # Down                             15                            180                                 1
    # NotSig                        14574                          12756                             14273
    # Up                              934                           2587                              1249
    # clusterinj_Granulocytes_vs_others clusterinj_Macrophages_vs_others clusterinj_Microglia_vs_others
    # Down                                 177                               25                             53
    # NotSig                             13399                            13159                          12675
    # Up                                  1947                             2339                           2795
    # clusterinj_Monocytes_vs_others clusterinj_Neuron_vs_others clusterinj_NK cells_vs_others
    # Down                               15                          41                            46
    # NotSig                          13790                       11681                         14915
    # Up                               1718                        3801                           562
    # clusterinj_Oligodendrocytes_vs_others clusterinj_OPC_vs_others clusterinj_T cells_vs_others
    # Down                                     127                      231                           34
    # NotSig                                 12339                    12676                        14814
    # Up                                      3057                     2616                          675
    # clusterinj_Vascular_vs_others
    # Down                              19
    # NotSig                         13647
    # Up                              1857
    # select the top 20 marker (up-regulated) genes for each cluster.
    top <- 20
    topMarkers <- list()
    for(i in 1:ncls) {
      ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
      up <- qlf[[i]]$table$logFC > 0
      topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
    }
    topMarkers <- unique(unlist(topMarkers))
    topMarkers
    lcpm <- cpm(y, log=TRUE)
    topMarkers20 <- lcpm[topMarkers, ]
    
    # Heatmap
    # Create the annotation dataframe with cluster labels
    annot <- data.frame(cluster = factor(cluster))
    rownames(annot) <- colnames(y)
    
    # Assign unique colors to each cluster label
    ann_colors <- list(cluster = cluster_colors)
    names(ann_colors$cluster) <- levels(cluster)
    
    # Plot with unique colors for each cluster
    pheatmap::pheatmap(topMarkers20, 
                       breaks = seq(-2, 2, length.out = 101),
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       scale = "row", 
                       cluster_cols = TRUE, 
                       border_color = "NA", 
                       fontsize_row = 5,
                       treeheight_row = 70, 
                       treeheight_col = 70, 
                       cutree_cols = 7,
                       clustering_method = "ward.D2", 
                       show_colnames = FALSE,
                       annotation_col = annot, 
                       annotation_colors = ann_colors)
    
  }
  EDGER exactTest
  {
    ### 1. Create pseudo-bulk samples: Seurat2PB
    # Pseudo-bulk samples are created by aggregating read counts together 
    # for all the cells with the same combination
    y <- Seurat2PB(SCT_objH_labels, sample="sample", cluster="combined_label") #for creating a pseudo-bulk DGEList object from a Seurat object.
    dim(y)
    # 28115    28
    counts <- y$counts
    samples <- y$samples
    
    group <- factor(c('ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl', 'ctrl',  'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj', 'inj'))
    y$samples$group <- group
    samples <- y$samples
    
    ### 2. Filtering and normalization
    # filter out library size below the threshold of 50,000
    summary(y$samples$lib.size)
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    # 98424   320117   786996  5350659  3724328 63581689  
    keep.samples <- y$samples$lib.size > 5e4 #not needed if libraries are bigger than this but run anyways
    keep.samples
    table(keep.samples)
    y <- y[, keep.samples]
    # Filter out lowly expressed genes
    keep.genes <- filterByExpr(y, group=y$samples$cluster,min.count=10, min.total.count=20)
    table(keep.genes)
    # FALSE  TRUE 
    # 12816 15299  # eliminamos esos genes FALSE y nos quedamos solo con los TRUE
    y <- y[keep.genes, , keep=FALSE]
    ### TMM normalization is performed to estimate effective library sizes.
    y <- normLibSizes(y)
    samples <- y$samples
    
    ### 3. Data exploration: explore the pseudo-bulk profiles
    # using a multi-dimensional scaling (MDS) plot
    # Plot MDS with viridis colors
    # Generate colors using turbo palette
    cluster <- as.factor(y$samples$cluster)
    cluster_colors <- turbo(length(unique(cluster)))
    plotMDS(y, method="logFC", col=cluster_colors[cluster])
    legend("bottomleft", as.character(unique(y$samples$group)), col=1:2, pch=20)
    
    plotMDS(y, pch=16, method= "logFC", col=cluster_colors[cluster], main="MDS")
    legend("bottomleft", as.character(unique(y$samples$cluster)), col=cluster_colors,
           pch=16, y.intersp = 0.6)
    
    group2 <- as.factor(y$samples$group)
    plotMDS(y, method="logFC", pch=16, col=c(1:2)[group2], main="MDS")
    legend("bottomleft", as.character(unique(y$samples$group)), col=1:2, pch=20)
    
  
    
    ### 4. Estimating the dispersion of the data
    # The first major step in the analysis of DGE data using the NB model is to estimate 
    # the dispersion parameter for each tag, a measure of the degree of inter-library variation for that tag. 
    # Estimating the common dispersion gives an idea of overall variability across the genome for this dataset.
    
    
    # 4.1 Naive dispersion
    y1 <- estimateCommonDisp(y, verbose=T)
    # Disp = 1.52906 , BCV = 1.2366
    names(y1)
    # [1] "counts"            "samples"           "genes"             "common.dispersion" "pseudo.counts"    
    # [6] "pseudo.lib.size"   "AveLogCPM"  
    y1 <- estimateTagwiseDisp(y1)
    names(y1)
    # [1] "counts"             "samples"            "genes"              "common.dispersion"  "pseudo.counts"     
    # [6] "pseudo.lib.size"    "AveLogCPM"          "prior.df"           "prior.n"            "tagwise.dispersion"
    # [11] "span" 
    plotBCV(y1)
    # Here we see that a single estimate for the coefficient of variation is a bad model since tagwise dispersion 
    # does not follow the model.
    
    # 4.2 GLM dispersion
    # Design a matrix
    design.mat <- model.matrix(~ 0 + y$samples$group)
    colnames(design.mat) <- levels(y$samples$group)
  
    # Estimate dispersion using the estimateDisp function, and the GLM method
    y2 <- estimateGLMCommonDisp(y,design.mat)
    y2 <- estimateGLMTrendedDisp(y2,design.mat, method="power")
    # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
    # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
    y2 <- estimateGLMTagwiseDisp(y2,design.mat)
    plotBCV(y2) # visualize with plotBCV
    
    # we model the tagwise dispersions based on the model derived from the glm model that we choose.
    # If we change the method power above to something else, 
    # the tagwise errors change to reflect that the method is different
    
    
    ### 5. Differential expression
    # 5.1 Exact Tests for Differences between Two Groups of Negative-Binomial Counts
    
    et21 <- exactTest(y1, pair=c(2,1))
    topTags(et21, n=10) #TOMANDO TODO SIN MIRAR POR POBLACION
    # Comparison of groups:  ctrl-inj 
    #                 gene     logFC   logCPM     PValue          FDR
    # Nudc-ps1   Nudc-ps1 -7.849146 2.066713 1.070110e-07 0.0009341291
    # Gm10736     Gm10736 -3.541407 5.166750 1.221164e-07 0.0009341291
    # Uts2           Uts2  7.314761 5.609979 5.557791e-07 0.0028342882
    # Cxcl13       Cxcl13 -5.864729 3.889212 8.194757e-07 0.0031342898
    # Serpina3n Serpina3n -4.970033 6.224948 6.884017e-06 0.0210637149
    # Gm50022     Gm50022 -9.237060 3.298198 1.452226e-05 0.0370293378
    # Gm33576     Gm33576  3.867439 1.370251 2.719743e-05 0.0594419237
    # Tk1             Tk1 -3.189917 3.339740 3.243943e-05 0.0620363632
    # Uts2b         Uts2b  5.283009 3.195219 3.759609e-05 0.0639091790
    # Mir5136     Mir5136  2.482683 3.379835 6.541880e-05 0.1000842242
    
    de1 <- decideTestsDGE(et21, adjust.method="BH", p.value = 0.05)
    summary(de1)
    # ctrl-inj
    # Down          5
    # NotSig    15293
    # Up            1
    
    et12 <- exactTest(y1, pair=c(1,2))
    topTags(et12, n=10) #TOMANDO TODO SIN MIRAR POR POBLACION
    # Comparison of groups:  inj-ctrl 
    #               gene     logFC   logCPM       PValue          FDR
    # Nudc-ps1   Nudc-ps1  7.849146 2.066713 1.070110e-07 0.0009341291
    # Gm10736     Gm10736  3.541407 5.166750 1.221164e-07 0.0009341291
    # Uts2           Uts2 -7.314761 5.609979 5.557791e-07 0.0028342882
    # Cxcl13       Cxcl13  5.864729 3.889212 8.194757e-07 0.0031342898
    # Serpina3n Serpina3n  4.970033 6.224948 6.884017e-06 0.0210637149
    # Gm50022     Gm50022  9.237060 3.298198 1.452226e-05 0.0370293378
    # Gm33576     Gm33576 -3.867439 1.370251 2.719743e-05 0.0594419237
    # Tk1             Tk1  3.189917 3.339740 3.243943e-05 0.0620363632
    # Uts2b         Uts2b -5.283009 3.195219 3.759609e-05 0.0639091790
    # Mir5136     Mir5136 -2.482683 3.379835 6.541880e-05 0.1000842242
    
    de2 <- decideTestsDGE(et12, adjust.method="BH", p.value = 0.05)
    summary(de2)
    # inj-ctrl
    # Down          1
    # NotSig    15293
    # Up            5
    
    de2tags12 <- rownames(y1)[as.logical(de2)] 
    plotSmear(et12, de.tags=de2tags12)
    abline(h = c(-2, 2), col = "blue")
    
    
    
    # 5.2 Genewise Negative Binomial Generalized Linear Models
    design.mat
    fit <- glmFit(y2, design.mat)
    lrt12 <- glmLRT(fit, contrast=c(1,-1))
    topTags(lrt12, n=10)
    # Coefficient:  1*ctrl -1*inj 
    #               gene     logFC   logCPM       LR       PValue          FDR
    # Nudc-ps1   Nudc-ps1 -7.849177 2.065741 37.46705 9.297022e-10 1.422351e-05
    # Gm10736     Gm10736 -3.540537 5.166610 29.69436 5.058182e-08 3.869256e-04
    # Uts2           Uts2  7.359668 5.610100 20.16890 7.089639e-06 2.969132e-02
    # Gm50022     Gm50022 -9.237742 3.297947 19.80831 8.560970e-06 2.969132e-02
    # Cxcl13       Cxcl13 -5.933886 3.887744 19.56888 9.703681e-06 2.969132e-02
    # Serpina3n Serpina3n -4.977420 6.224621 18.58127 1.628125e-05 4.151448e-02
    # Mir5136     Mir5136  2.479540 3.380339 18.02825 2.176511e-05 4.756921e-02
    # Tk1             Tk1 -3.185350 3.339874 17.25133 3.274678e-05 6.262412e-02
    # Wdr87-ps   Wdr87-ps  2.587380 2.381162 16.56905 4.691025e-05 7.974222e-02
    # Gm33576     Gm33576  3.904907 1.372197 15.69044 7.460009e-05 1.141307e-01
    
    deGLM <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
    summary(deGLM)
    # 1*ctrl -1*inj
    # Down               5
    # NotSig         15292
    # Up                 2
    
    deGLMtags12 <- rownames(y2)[as.logical(deGLM)]
    plotSmear(lrt12, de.tags=deGLMtags12)
    abline(h = c(-2, 2), col = "blue")
    
    
    #### CHINO
     ### Marker genes identification
    ncls <- nlevels(cluster)
    contr <- rbind( matrix(1/(1-ncls), ncls, ncls),
                    + matrix(0, ncol(design)-ncls, ncls) )
    diag(contr) <- 1
    contr[1,] <- 0
    rownames(contr) <- colnames(design)
    colnames(contr) <- paste0("cluster", levels(cluster))
    
    qlf <- list()
    for(i in 1:ncls){
      qlf[[i]] <- glmQLFTest(fit, contrast=contr[,i])
      qlf[[i]]$comparison <- paste0("cluster", levels(cluster)[i], "_vs_others")
    }
    
    CtrlNeuron <- qlf[[1]][["table"]]
    InjNeuron <- qlf[[2]][["table"]]
    
    # show the numbers of DE genes under each comparison
    dt <- lapply(lapply(qlf, decideTestsDGE), summary)
    dt.all <- do.call("cbind", dt)
    dt.all
    
    # clusterctrl_vs_others clusterinj_vs_others
    # Down                       0                    0
    # NotSig                  9697                 9697
    # Up                         0                    0
    
    # select the top 20 marker (up-regulated) genes for each cluster.
    top <- 20
    topMarkers <- list()
    for(i in 1:ncls) {
      ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
      up <- qlf[[i]]$table$logFC > 0
      topMarkers[[i]] <- rownames(y)[ord[up][1:top]]
    }
    topMarkers <- unique(unlist(topMarkers))
    topMarkers
    lcpm <- cpm(y, log=TRUE)
    topMarkers20 <- lcpm[topMarkers, ]
    
    # Heatmap
    # Create a color palette with enough colors for all clusters
    # Create the annotation dataframe with cluster labels
    annot <- data.frame(cluster = factor(cluster))
    rownames(annot) <- colnames(y)
    
    # Assign unique colors to each cluster label
    ann_colors <- list(cluster = cluster_colors)
    names(ann_colors$cluster) <- levels(cluster)
    
    # Plot with unique colors for each cluster
    pheatmap::pheatmap(topMarkers20, 
                       breaks = seq(-2, 2, length.out = 101),
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       scale = "row", 
                       cluster_cols = TRUE, 
                       border_color = "NA", 
                       fontsize_row = 5,
                       treeheight_row = 70, 
                       treeheight_col = 70, 
                       cutree_cols = 7,
                       clustering_method = "ward.D2", 
                       show_colnames = FALSE,
                       annotation_col = annot, 
                       annotation_colors = ann_colors)
    
  }
}
  
  
  MAST
  # MAST is based on a zero-inflated negative binomial model. 
  # It tests for differential expression using a hurdle model to combine tests of discrete (0 vs not zero) 
  # and continuous (non-zero values) aspects of gene expression
  BiocManager::install('MAST')
  library(MAST)
  SCT_objH_labels
  # An object of class Seurat 
  # 53205 features across 56158 samples within 2 assays 
  # Active assay: SCT (25090 features, 3000 variable features)
  # 3 layers present: counts, data, scale.data
  # 1 other assay present: RNA
  # 5 dimensional reductions calculated: pca, umap.unintegrated, harmony, umap.harmony, umap.harmony3
  DefaultAssay(SCT_objH_labels) <- "RNA"
  SCT_objH_labels  
  countsH <- GetAssayData(object = SCT_objH_labels, layer = "counts", assay="RNA")
  
  sce <- SingleCellExperiment(list(countsH), colData = SCT_objH_labels@meta.data)
  
  sce <- as.SingleCellExperiment(SCT_objH_labels, layer = "counts", assay="RNA")
  logcounts(sce) <- log2(counts(sce) + 1)
  sca <- SceToSingleCellAssay(sce)
  freq_expressed = 0.2
  expressed_genes <- freq(sca) > freq_expressed
  sca <- sca[expressed_genes,]
  cdr <-colSums(assay(sca)>0)
  sca@colData[["cngeneson"]] <- scale(cdr)
  cond <-factor(colData(sca)[[group]])
  cond <-relevel(cond,level_group)
  sca@colData[["cond"]] <- cond
  sca
  
  
  10.2 GENE SET ENRICHMENT AND PATHWAY ANALYSIS
  10.3 CELL TO CELL COMMUNICATION
  