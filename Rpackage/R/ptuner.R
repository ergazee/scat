## This function trys to find the optimal value of two parameters in Seurat Standard Workflow (nFeatures, resolution)
## when using a reference scRNA-seq data to identify the cell type of unmark data.

#' @param datalist a list of reference and unmark data in Seurat object format
#' @param Normalize_Method choose different normalize method by set "Standard" or "SCT"
#' @param Nmin the minimum of Seurat parameter 'nFeatures'
#' @param Nmax the maximum of Seurat parameter 'nFeatures'
#' @param Nstep step size of the Seurat parameter 'nFeatures' in the for loop
#' @param Rmin the minimum of Seurat parameter 'resolution'
#' @param Rmax the maximum of Seurat parameter 'resolution'
#' @param Rstep step size of the Seurat parameter 'resolution' in the for loop
#' @param Dims number of selected PCA components
#'
#'
ptuner <- function(data.list,
                   Normalize_Method = "Standard",
                   Nmin          = 2000,
                   Nmax          = 10000,
                   Nstep         = 1000,
                   Rmin          = 0.2,
                   Rmax          = 1.0,
                   Rstep         = 0.1,
                   Dims          = 30) {

  ratios <- data.frame()

  for (i in seq(Nmin, Nmax, Nstep)){


    if(Normalize_Method == "Standard"){
      for (m in 1:length(data.list)){
        data.list[[m]] <- FindVariableFeatures(data.list[[m]], selection.method = "vst", nfeatures = i)
      }
      data.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:Dims)
      data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:Dims)
      DefaultAssay(data.integrated) <- "integrated"
      data.integrated <- ScaleData(data.integrated)
    }
    else if(Normalize_Method == "SCT"){
      data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = i)
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

    for (j in seq(Rmin, Rmax, Rstep)){
      data.integrated <- FindClusters(data.integrated, resolution = j)

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
      ratios[as.character(i),as.character(j)] <- ratio
    }
  }
  write.csv(ratios, as.character(paste("Coincidence",Nmin,Nmax,Rmin,Rmax,".csv",sep = "_")))
  bestp <- as.list(as.numeric(which(ratios == max(ratios), arr.ind = TRUE)))
  bestp[[1]] <- as.numeric(rownames(ratios)[bestp[[1]]])
  bestp[[2]] <- as.numeric(colnames(ratios)[bestp[[2]]])

  data.integrated <- finalresults(data.list = data.list,
                                  Normalize_Method = Normalize_Method,
                                  nFeatures = bestp[[1]],
                                  Resolution = bestp[[2]],
                                  Dims = Dims)
  return(data.integrated)
  message(paste("The best paramater for your datasets is: nfeatures=", bestp[[1]], ", resolution=", bestp[[2]], ". the max coincidence is: ", max(ratios), sep = ''))
}

