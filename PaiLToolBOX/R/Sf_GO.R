#' @title GO enrichment of Spodoptera frugiperda
#
#' @param genelist should be factors
#' @param method will be set to "holm" as default if not chosen
#' @param withanc means that genes will be enriched into ancestor GO terms
#' @description To enrich genes into Sf genome (from BIPAA).
#' @author Pai Li
#' @return GO enrichment Table
#' @export
#
#
#


Sf_GO<- function(genelist, method=p.adjust.methods, withanc= TRUE)
{
  if (missing(genelist) ||length(genelist)<= 1 ||class(genelist)!= "factor")
  {
    stop("Genelist needs to be factor with no less than 2 genes")
  }
  method= match.arg(method)

  if (withanc) {
    g<- Sf_with_anc[Sf_with_anc$Sf %in% genelist, ]

    if (nrow(g)<=0) {
      stop("Genelist is invalid, check whether there is sapce in gene names")
    }

    hyperdf<- data.frame()

    for (i in 1:length(unique(g[,3]))) {
      hyperdf[i,1]<- length(unique(g[g[,3]== unique(g[,3])[i],1]))
      hyperdf[i,2]<- length(unique(Sf_with_anc[Sf_with_anc[,3]==unique(g[,3])[i],1]))
      hyperdf[i,3]<- length(unique(g[,1]))
      hyperdf[i,4]<- length(unique(Sf_with_anc[,1]))
    }

    p_value<- data.frame()
    for (i in 1:nrow(hyperdf)) {
      p_value[i,1]<- phyper((hyperdf[i,1]-1), hyperdf[i,2], (hyperdf[i,4]-hyperdf[i,2]), hyperdf[i,3], lower.tail = FALSE)
    }

    q_value<- p.adjust(p_value[,1], method = method)

  } else {

    g<- Sf_GOdb[Sf_GOdb$Sf %in% genelist, ]

    if (nrow(g)<=0) {
      stop("Genelist is invalid, check whether there is sapce in gene names")
    }
    hyperdf<- data.frame()

    for (i in 1:length(unique(g[,3]))) {
      hyperdf[i,1]<- length(unique(g[g[,3]== unique(g[,3])[i],1]))
      hyperdf[i,2]<- length(unique(Sf_GOdb[Sf_GOdb[,3]==unique(g[,3])[i],1]))
      hyperdf[i,3]<- length(unique(g[,1]))
      hyperdf[i,4]<- length(unique(Sf_GOdb[,1]))
    }

    p_value<- data.frame()
    for (i in 1:nrow(hyperdf)) {
      p_value[i,1]<- phyper((hyperdf[i,1]-1), hyperdf[i,2], (hyperdf[i,4]-hyperdf[i,2]), hyperdf[i,3], lower.tail = FALSE)
    }

    q_value<- p.adjust(p_value[,1], method = method)

  }

  GOs<- unique(g[,3])
  Rich_score<- hyperdf[,1]/hyperdf[,2]
  Fold_enrichment<- hyperdf[,1]*hyperdf[,4]/(hyperdf[,2]*hyperdf[,3])
  genes_in_GO<- data.frame()
  for (i in 1:length(GOs)) {
    for (j in 1:length(g[g[,3]==GOs[i], 1])) {
      genes_in_GO[i,j]<- g[g[,3]==GOs[i], 1][j]

    }


  }


  result<- data.frame(GOs, p_value, q_value, hyperdf, Rich_score, Fold_enrichment, genes_in_GO)
  result<- merge(GO_term_ID, result, by.y= "GOs", by.x= "GO_ID", all.y= TRUE)
  names(result)[1:9]<- c("GO_ID", "GO_term", "Ontology", "p_value", "q_value", "List Hits", "Pop Hits", "List Total", "Pop Total")
  result<- result[order(result[,5]),]
}
