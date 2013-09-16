
buildSBMLFromBiGG <- function(reactions.filename, model.id=character(0), model.name=character(0)){
  tab <- read.delim(reactions.filename, header=TRUE, sep="\t")
  
  ##parse equations and create reaction objects
  all.reactions <- list()
  all.compartment.names <- vector()
  all.species.names <- vector()
  for (i in seq_len(nrow(tab))){
    row <- tab[i,]  
    ##get reaction equation
    eq <- as.character(unlist(row$equation))
    
    ## There are two different formats for reaction equations:
    ## -- "[c] : ru5p-D <==> xu5p-D"   ##one compartment for all species
    ## -- "(2) ficytC[m] + (2) h[m] + q10h2[m] --> (2) focytC[m] + (4) h[c] + q10[m]" ##Compartment for each species
    
    if (substr(eq, 1, 1)=="["){
      ##if equation starts with "[compartment] : "  
      reactants <- gsub(".*?:", "", gsub("[-|<].+", "", eq))
      products <- gsub(".+>\\s?(.+)", "\\1", eq) 
      compartments <- gsub("\\[(.*)\\].+", "\\1", eq)
    } else {  
      ##if equation does not start with "[compartment] : "
      s <- gsub("[\\[].+?[\\]]", "", eq, perl=TRUE) ##equation without brackets at the beginning
      reactants <- gsub("\\s+[-|<|>].*"  , "", s)  
      products <- gsub(".*[-|<|>]"  , "", s)   
      compartments <- unlist(regmatches(eq, gregexpr("\\[.+?\\]", eq)))
      compartments <- gsub("\\]|\\[", "", compartments)
    }
    
    ##reactant and product vectors still containing possible stoichiometries
    reactants <- unlist(strsplit(reactants, "\\+"))
    products <- unlist(strsplit(products, "\\+"))
    st.list <- regmatches(c(reactants, products), gregexpr("\\(.+?\\)", c(reactants, products)))
    stoichiometries <- sapply(st.list, function(x){
      if(length(x)==0)
        1 ##if there is no stichiometry given, set it to 1
      else
        as.numeric(gsub("\\(|\\)", "", x))
    })
    
    ##Clear reactant and product vectors from stoichiometries
    reactants <- as.vector(sapply(reactants, function(x)gsub("\\s", "", gsub("\\(.+?\\)", "", x))))
    products <- as.vector(sapply(products, function(x)gsub("\\s", "", gsub("\\(.+?\\)", "", x))))
    
    
    ##Assign compartment to each reactant and product if it is not yet assigned to all reactants and products
    if (length(compartments)==1 & length(reactants) + length(products)>0) 
      compartments <- rep(compartments, length(reactants) + length(products))
    
    all.compartment.names <- append(all.compartment.names, compartments)
    all.species.names <- append(all.species.names, c(reactants, products))
    ##get reversibility
    reversible <- length(grep("<==>|<-->", eq))>0
     
    
    ##Build SBML objects
    ##species references for reactants and products
    reactant.refs <- mapply(function(r, c, s)new("SpeciesReference", id=paste(r, c, sep="_"), species=paste(r, c, sep="_"), stoichiometry=s)
                            , reactants, compartments[seq_along(reactants)], stoichiometries[seq_along(reactants)])
    
    product.refs <- mapply(function(r, c, s)new("SpeciesReference", id=paste(r, c, sep="_"), species=paste(r, c, sep="_"), stoichiometry=s)
                           , products, compartments[-seq_along(reactants)], stoichiometries[-seq_along(reactants)])
    
    ##reaction object
  all.reactions[[i]] <- new("Reaction", id=as.character(row$abbreviation), name=as.character(row$name), reactants=reactant.refs, products=product.refs, reversible=reversible)
  }
  names(all.reactions) <- sapply(all.reactions, id)
  
  ##Build SBML objects for species and compartments
  all.compartments <- lapply(unique(all.compartment.names), function(c)new("Compartment", id=c, name=c))
  all.species <- apply(unique(cbind(all.species.names, all.compartment.names)), 1,
                       function(x)new("Species", id=paste(x[1], x[2], sep="_"), name=paste(x[1], x[2], sep="_"), compartment=x[2]))
  names(all.species) <- sapply(all.species, name)
  
  return(new("Model", species=all.species, reactions=all.reactions, compartments=all.compartments,
             id=ifelse(length(model.id)==0, reactions.filename, model.id),
             name=ifelse(length(model.name)==0, reactions.filename, model.name)))
}


