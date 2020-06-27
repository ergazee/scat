## This function use Sankey plot show the label change of cells.

#' @param data
#' @param label1
#' @param label2
#' @param label3
#'
Sankeyplot <- function(data,
                   label1,
                   label2,
                   label3) {
  Sankey_1 <- data.frame("cells" = rep(1,ncol(data)),
                         "freq" = seq(1,ncol(data),1),
                         "method" = rep(label1,ncol(data)),
                         "CellType" = data@meta.data[,label1])
  Sankey_2 <- data.frame("cells" = rep(1,ncol(data)),
                         "freq" = seq(1,ncol(data),1),
                         "method" = rep(label2,ncol(data)),
                         "CellType" = data@meta.data[,label2])
  Sankey_3 <- data.frame("cells" = rep(1,ncol(data)),
                         "freq" = seq(1,ncol(data),1),
                         "method" = rep(label3,ncol(data)),
                         "CellType" = data@meta.data[,label3])


  Sankey <- rbind(Sankey_1, Sankey_2, Sankey_3)

  ggplot(Sankey,
               aes(x = method, stratum = CellType, alluvium = freq,
                   y = cells,
                   fill = CellType, label = CellType)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12))

  }

