#' @param data
#' @param label1
#' @param label2
#'
#'
Barplot <- function(data,
                   label1,
                   label2){

  pdata <- melt(data@meta.data[,c(label1,label2)])
  colnames(pdata) <- c(label1,"variable",label2)
  ggplot(pdata, aes_string(x=label1, y=label2, fill=label1)) +
  geom_boxplot(alpha=0.8) +
  geom_point(shape=16, size=0.6, position = position_jitter(0.4), color="black", alpha=1) +
  theme(plot.title=element_text(hjust=0.5, size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


  }

