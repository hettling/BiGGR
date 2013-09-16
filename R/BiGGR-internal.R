##Get a subset of a model in which all reactions, species, compartments are associated
##  with the reactions in argument reaction.ids (a list or vector of character objects).
##  Returns a rsbml Model object. This is an internal function 
.buildSubModel <- function(model, reaction.ids){
  ##Get all reactions
  relevant.reactions <- model@reactions[reaction.ids]
  ##Get all species that are associated with given reactions
  relevant.species.names <- unique(unlist(sapply(relevant.reactions, function(r){
    lapply(c(r@reactants, r@products), function(sr){
      model@species[[sr@species]]@id
    })
  })))
  ##remove duplicates
  relevant.species <- model@species[relevant.species.names]
    ##Get the compartments from the species
  relevant.compartment.names <- unique(sapply(relevant.species, function(s)s@compartment))
  relevant.compartments <- model@compartments[relevant.compartment.names]

  ##Build SBML model from reactions, species and compartments
  new("Model", species=relevant.species, reactions=relevant.reactions, compartments=relevant.compartments)
}    
