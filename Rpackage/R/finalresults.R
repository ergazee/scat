## This function use the optimal value of two parameters in Seurat (nFeatures, resolution) to output results
## when using a reference scRNA-seq data to identify the cell type of unmark data.

#' @param datalist a list of reference and unmark data in Seurat object format
#' @param Normalize_Method choose different normalize method by set "Standard" or "SCT"
#' @param nFeatures the best value of Seurat parameter 'nfeatures'
#' @param Resolution the best value of Seurat parameter 'resolution'
#' @param Dims number of selected PCA components
#'
#'
finalresults <- function(data.list,
                   Normalize_Method = "Standard",
                   nFeatures,
                   Resolution,
                   Dims = 30) {

    if(Normalize_Method == "Standard"){
      for (m in 1:length(data.list)){
        data.list[[m]] <- FindVariableFeatures(data.list[[m]], selection.method = "vst", nfeatures = nFeatures)
      }
      data.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:Dims)
      data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:Dims)
      DefaultAssay(data.integrated) <- "integrated"
      data.integrated <- ScaleData(data.integrated)
    }
    else if(Normalize_Method == "SCT"){
      data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = nFeatures)
      data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, verbose = FALSE)
      data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
      data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)
    }
    else{
      stop("please set a normalize method: Standard or SCT")
    }

    data.integrated <- RunPCA(data.integrated, npcs = Dims)
    data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:Dims)
    data.integrated <- RunTSNE(data.integrated, reduction = "pca", dims = 1:Dims)
    data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:Dims)

    data.integrated <- FindClusters(data.integrated, resolution = Resolution)


    for(i in c("umap", "tsne")){
      pdf(as.character(paste("1.1.Cluster",Normalize_Method, nFeatures, Resolution, i,".pdf",sep = "_")),width = 10, height = 6)
      # Samples
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "type",
                    pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "type", split.by = "type",
                    ncol = 2, pt.size =0.5))

      # refLabel
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "refLabel",
                    label = TRUE, repel = TRUE, pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "refLabel", split.by ="type",
                    label = FALSE, repel = TRUE, ncol = 2, pt.size =0.05))
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "refLabel", split.by ="refLabel",
                    label = FALSE, repel = TRUE, ncol = 4, pt.size =0.05))

      # Seurat
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "seurat_clusters",
                    label = TRUE, repel = TRUE, pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "seurat_clusters", split.by = "type",
                    label = TRUE, repel = TRUE, ncol = 2, pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i,
                    group.by = "seurat_clusters", split.by = "seurat_clusters",
                    label = FALSE, repel = TRUE, ncol = 4, pt.size =0.5))

      dev.off()
    }

    dt <- data.frame("ID" = names(data.integrated@active.ident),
                     "original" = unlist(lapply(data.list, function(x) as.character(x@meta.data$refLabel))),
                     "seurat" = data.integrated@meta.data$seurat_clusters)

    datanames <- unlist(lapply(data.list, function(x) if(sum(!duplicated(x@meta.data$refLabel))==1) x@project.name))
    dt1 <- as.data.frame(unclass(t(table(dt[,2:3]))))
    dt1 <- dt1 %>% select(-datanames,datanames)

    dt2 <- dt1
    for(m in 1:dim(dt2)[1]){
      for(n in 1:dim(dt2)[2]){
        if (dt2[m,n] < max(dt1[m,]))
          dt2[m,n] = 0
      }
    }
    ratio <- round(sum(dt2[,!(names(dt2) %in% datanames)])/sum(dt1[,!(names(dt1) %in% datanames)]),4)
    dt3 <- dt1
    dt3["MaxSum",] <- as.numeric(colSums(dt2))
    dt3["AllSum",] <- as.numeric(colSums(dt1))
    dt3["cRatio",] <- round(dt3["MaxSum",]/dt3["AllSum",],2)

    rowsums1 <- as.numeric(rowSums(dt2))
    rowsums2 <- as.numeric(rowSums(dt1))

    length(rowsums1) <- dim(dt3)[1]
    length(rowsums2) <- dim(dt3)[1]

    dt3[,"MaxSum"] <- rowsums1
    dt3[,"AllSum"] <- rowsums2

    rRatio <- round(dt3[,"MaxSum"]/dt3[,"AllSum"],2)
    dt3[,"rRatio"] <- rRatio
    rm(rowsums1,rowsums2,rRatio)

    dt3[dim(dt3)[1],dim(dt3)[2]] <- ratio
    write.csv(dt3, as.character(paste("Final_table",nFeatures,Resolution,".csv",sep = "_")))


    NewLabel <- apply(dt2, 1, function(t) colnames(dt2)[which.max(t)])
    NewLabel <- replace(NewLabel, names(NewLabel), make.names(NewLabel,unique=T))
    data.integrated <- RenameIdents(data.integrated, NewLabel)
    data.integrated@meta.data$celltype <- as.character(data.integrated@active.ident)


    for(i in c("umap", "tsne")){
      pdf(as.character(paste("2.Rename",Normalize_Method, nFeatures, Resolution, i,".pdf",sep = "_")),width = 10, height = 6)
      print(DimPlot(data.integrated, reduction = i, label = TRUE, pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i, label = TRUE, split.by = "type", ncol = 2, pt.size =0.5) + NoLegend())
      print(DimPlot(data.integrated, reduction = i, label = FALSE, group.by = c("type", "refLabel"), ncol = 2 ,pt.size =0.5))
      print(DimPlot(data.integrated, reduction = i, label = FALSE, group.by = c("seurat_clusters", "celltype"), ncol = 2 ,pt.size =0.5))
      dev.off()
    }


    CellType <- data.frame("barcodes" = names(data.integrated@active.ident),
                             "type" = data.integrated@meta.data$type,
                             "seurat_cluster" = data.integrated@meta.data$seurat_clusters,
                             "original_labels" = data.integrated@meta.data$refLabel,
                             "celltype" = as.character(data.integrated@active.ident),
                             stringsAsFactors = FALSE)
    rownames(CellType) <- CellType$barcodes
    write.csv(CellType, as.character(paste("CellType",nFeatures,Resolution,".csv",sep = "_")))


    return(data.integrated)
  }

