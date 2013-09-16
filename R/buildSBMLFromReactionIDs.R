##construct an SBML model from given reaction IDs
buildSBMLFromReactionIDs <- function(reaction.ids, database){
  model <- database@model
  notfound.ids <- reaction.ids[which(! reaction.ids %in% lapply(model@reactions, id))]
  if (length(notfound.ids>0))
    warning("Warning! Omitting the following reactions not found in database : ", notfound.ids, "\n")
  .buildSubModel(model, reaction.ids)  
}
