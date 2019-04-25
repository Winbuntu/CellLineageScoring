
CellLineageScoring <- function (object, TE.genes, PE.genes, EPI.genes, set.ident = FALSE) 
{
  enrich.name <- "Cell lineage"
  genes.list <- list(TE.Score = TE.genes, PE.Score = PE.genes, EPI.Score = EPI.genes)
  object.cc <- AddModuleScore(object = object, genes.list = genes.list, 
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list, 
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "TE", second = "PE", third = "EPI",null = "Untypical") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      return(c(first, second, third)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "TE.Score", "PE.Score", "EPI.Score",
                               "Lineage")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("TE.Score", "PE.Score", "EPI.Score","Lineage")]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = "old.ident")
    object <- SetAllIdent(object = object, id = "Lineage")
  }
  return(object)
}
