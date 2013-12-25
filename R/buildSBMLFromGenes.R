##construct an SBML model including reactions involving (a) given gene(s) 
buildSBMLFromGenes <- function(query, database, logical.fun="any"){
  fun <- match.fun(logical.fun)
  all.gas <- extractGeneAssociations(database)
  relevant.reaction.ids <- unique(names(all.gas[which(apply(sapply(query, grepl, all.gas), 1, fun))]))
  if (length(relevant.reaction.ids)==0){
    cat("No reactions found for query ", query, "\n")
    return(NULL)
  }
  .buildSubModel(database@model, relevant.reaction.ids)                                                
}
