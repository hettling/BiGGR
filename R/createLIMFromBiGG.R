createLIMFromBiGG <- function(reactions.filename, ...){
  model <- buildSBMLFromBiGG(reactions.filename)
  createLIMFromSBML(model, ...)
}
