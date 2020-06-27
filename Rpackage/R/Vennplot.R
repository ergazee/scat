#' @param data a list of reference and unmark data in Seurat object format
#' @param name the best value of Seurat parameter 'nfeatures'
#'

VN2 <- function(data, name){
  if (is.list(data) == TRUE)
    venn.plot <- venn.diagram(
      x = data,
      filename = name, imagetype = "png",
      lwd = 3,
      fill = c("cornflowerblue", "darkorchid1"),
      alpha = 0.6,
      label.col = "white",
      cex = 1.5,
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("cornflowerblue", "darkorchid1"),
      cat.cex = 2,
      cat.fontfamily = "serif",
      cat.fontface = "bold",
      margin = 0.05,
      cat.dist = c(0.03, 0.03),
      cat.pos = c(-20, 20)
    )
  else warning("data input must be a list like this: list(A = M$A, B = M$B)")
}


VN3 <- function(data, name){
  if (is.list(data) == TRUE)
    venn.plot <- venn.diagram(
      x = data,
      filename = name, imagetype = "png",
      col = "transparent",
      fill = c("red", "blue", "green"),
      alpha = 0.5,
      label.col = c("darkred", "white", "darkblue", "white",
                    "white", "white", "darkgreen"),
      cex = 2.5,
      fontfamily = "serif",
      fontface = "bold",
      cat.default.pos = "text",
      cat.col = c("darkred", "darkblue", "darkgreen"),
      cat.cex = 2.5,
      cat.fontfamily = "serif",
      cat.dist = c(0.06, 0.06, 0.03),
      cat.pos = 0
    )
  else warning("data input must be a list like this: list(A = M$A, B = M$B)")
}


VN4 <- function(data, name){
  if (is.list(data) == TRUE)
    venn.plot <- venn.diagram(
      x = data,
      filename = name, imagetype = "png",
      col = "black",
      lty = "dotted",
      lwd = 3,
      fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
      alpha = 0.50,
      label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                    "white", "white", "darkblue", "white",
                    "white", "white", "white", "darkgreen", "white"),
      cex = 2.0,
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
      cat.cex = 1.8,
      cat.fontface = "bold",
      cat.fontfamily = "serif"
    )
  else warning("data input must be a list like this: list(A = M$A, B = M$B)")
}


VN5 <- function(data, name){
  if (is.list(data) == TRUE)
    venn.plot <- venn.diagram(
      x = data,
      filename = name, imagetype = "png",
      lty = "dotted",
      lwd = 2,
      col = "black",
      fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
      alpha = 0.60,
      cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
      cat.cex = 0.8,
      cat.fontface = "bold",
      margin = 0.07,
      cex = 0.8
    )
  else warning("data input must be a list like this: list(A = M$A, B = M$B)")
}
