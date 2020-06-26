## This function will process the expression matrix data into Seurat objects
## and add a column of 'refLabel' labels to it.

#' @param mat a data.frame object of the cell expression matrix
#' @param name sample name
#' @param refLabel a data.frame object with two columns, one for cell index and one for CellType
#'
#'
datainput <- function(mat,
                      name,
                      refLabel = name){
  data <- CreateSeuratObject(counts = mat, project = name, min.cells = 5 )
  data$type <- name
  if(class(refLabel) != "character"){
    data@meta.data$refLabel <- refLabel[rownames(data@meta.data),"refLabel"]
    if(sum(is.na(data@meta.data$refLabel)) > 0){
      stop("Your cell matrix does not match the cell index of refLabel.")
    }
  } else{
    data@meta.data$refLabel <- name
    }
  return(data)
}
