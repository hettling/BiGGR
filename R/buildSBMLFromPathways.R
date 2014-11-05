##construct an SBML model for (a) given pathway(s) 
buildSBMLFromPathways <- function(query, database, match.exact=TRUE){
model <- database@model
all.pathways <- extractPathways(database)   
##Get all reactions in specific pathway(s)
relevant.reaction.ids <- unique(unlist(lapply(query, function(x){
if (match.exact)
  names(model@reactions[names(which(all.pathways==x))])
else
  names(model@reactions[grep(x, all.pathways, ignore.case=TRUE)])
})))         
if (length(relevant.reaction.ids)==0){
message("No Recon Pathways found for query ", query, "\n")
return(NULL)
}
.buildSubModel(model, relevant.reaction.ids)
}
