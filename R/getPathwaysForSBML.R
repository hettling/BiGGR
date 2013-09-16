##A function to extract all relevant pathways from a model given a specific database
getPathwaysForSBML <- function(model, database){
  all.pathways <- extractPathways(database)
  unique(sapply(lapply(model@reactions, function(x)x@id), function(y)all.pathways[[y]]))
}
