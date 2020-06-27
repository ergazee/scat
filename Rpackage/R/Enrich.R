#' @param data.markers
#' @param orgdb
#' @param database
#' @param organism
#'
Enrich <- function(data.markers,
                   database = "KEGG",
                   orgdb = "org.Hs.eg.db",
                   organism = 'hsa'
                  ){

  if(database == "GO"){
    for (i in names(table(data.markers$cluster))) {
      x <- (data.markers[data.markers$cluster==as.character(i),]$gene)
      x <- bitr(x, fromType="SYMBOL", toType = c("ENTREZID","ENSEMBL"), OrgDb = orgdb)
      go <- enrichGO(x$ENTREZID, OrgDb = orgdb, ont='ALL', pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,keyType = 'ENTREZID')
      write.csv(summary(go),paste(i,"GO_Pathway_enrichment.csv", sep = "_"),row.names =FALSE)
      pdf(paste(i,"GO_Pathway_enrichment.pdf", sep = "_"),width = 8,height = 6)
      print(dotplot(go,showCategory=10))
      dev.off()
    }
  }
  else if(database == "KEGG"){
    for (i in names(table(data.markers$cluster))) {
      x <- (data.markers[data.markers$cluster==as.character(i),]$gene)
      x <- bitr(x, fromType="SYMBOL", toType = c("ENTREZID","ENSEMBL"), OrgDb = orgdb)
      kegg <- enrichKEGG(x$ENTREZID, organism = organism, keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                         minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
      write.csv(as.data.frame(kegg),paste(i,"KEGG_Pathway_enrichment.csv", sep = "_"),row.names =FALSE)
      pdf(paste(i,"KEGG_Pathway_enrichment.pdf", sep = "_"),width = 8,height = 6)
      print(dotplot(kegg,showCategory=10))
      dev.off()
    }
  }
  else{
    stop("please choose a database: GO or KEGG")
  }
}


