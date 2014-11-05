#Generate a LIM model file for FBA(flux balance analysis) where model file is an SBML file consisting of reaction and metabolite lists
#maximize is a character string giving name of reaction tag to be maximized, 
#equation_var is a character string indicating the start of the reaction
# equation_value is a numeric value the initial value of the equation_var
# constraints is a character string in the format for example "[0,1000]" where 0 is the minimal and 1000 is maximum
# externals is a character vector provided by user

createLIMFromSBML <- function(model, maximize, equations=NULL, inequalities=NULL, constraints=NULL, externals=NULL, file.name="model.lim"){  
##REACTIONS
extract.eq <- function(r){
reacts <- sapply(r@reactants, function(p)p@species)
reac.stoich <- sapply(r@reactants, function(p)p@stoichiometry)
lhs.term <- paste(mapply(function(p, s){if(s!=1)paste(s, p, sep=" * ") else p}, reacts, reac.stoich ), collapse=" + ")
prods <- sapply(r@products, function(p)p@species)
prod.stoich <- sapply(r@products, function(p)p@stoichiometry)
rhs.term <- paste(mapply(function(p, s){if(s!=1)paste(s, p, sep=" * ") else p}, prods, prod.stoich ), collapse=" + ")
paste(r@id, ":", paste(lhs.term, rhs.term, sep=ifelse(r@reversible, " <-> ", " -> ")))
}

reactions <- as.vector(sapply(model@reactions, extract.eq))      
write("###   REACTION", file=file.name)
  
##check reactions containing  hyphens
if(isTRUE(grep("[a-z]-", reactions))==TRUE){
rr <- grep("[a-z]-", reactions)
react <- reactions[-rr]
react1 <- reactions[rr]
rr1 <- sub("-", "", react1)
reactions <- c(rr1, react)
}

write(reactions, file=file.name, append=TRUE)
write("### END REACTION", file=file.name, append=TRUE)
write("\n\n", file=file.name, append=TRUE)


##OPTIMIZE MINIMUM OR MAXIMUM REACTION
write("### MAXIMISE", file=file.name, append=TRUE)
write(maximize, file=file.name, append=TRUE) 
write("### END MAXIMISE", file=file.name, append=TRUE)
write("\n\n", file=file.name, append=TRUE)

##CONSTRAINTS
if (hasArg(constraints)){
if (length(constraints)!=3){
  warning("Warning: 'constraints' argument must be a list of length three. See ?createLIMFromSBML for details!\n")
} else {      
  write("###   CONSTRAINTS", file=file.name, append=TRUE)
  write(paste(constraints[[1]], " = [", constraints[[2]], ",", constraints[[3]], "]", sep=""), file=file.name, append=TRUE)
  write("###   END CONSTRAINTS",file=file.name,append=TRUE) 
  write("\n\n", file=file.name, append=TRUE)
}
}

##EQUATIONS
if (hasArg(equations)){
if (length(equations) != 2){
  warning("Warning: 'equations' argument must be a list of length two. See ?createLIMFromSBML for details!\n")
} else {
  write("###   EQUATIONS", file=file.name,append=TRUE)
  write(paste(equations[[1]], equations[[2]], sep=" = "), file=file.name, append=TRUE) 
  write("###   END EQUATIONS", file=file.name, append=TRUE)     
  write("\n\n", file=file.name, append=TRUE)              
}
}

##INEQUALITIES
if (hasArg(inequalities)){
if (length(inequalities) != 3){
  warning("Warning: 'inequalities' argument must be a list of length three. See ?createLIMFromSBML for details!\n")
} else {
  write("###   INEQUALITIES", file=file.name,append=TRUE)
  write(paste(inequalities[[1]], inequalities[[2]], sep=inequalities[[3]]), file=file.name, append=TRUE) 
  write("###   END INEQUALITIES", file=file.name, append=TRUE)     
  write("\n\n", file=file.name, append=TRUE)              
}
}

##COMPONENTS
all.components <- as.vector(sapply(model@species, function(s)s@id))
write("###   COMPONENTS", file=file.name, append=TRUE)    
components <- setdiff(all.components, externals)
write(components, file=file.name, append=TRUE)
write("###   END COMPONENTS", file=file.name, append=TRUE)    
write("\n\n", file=file.name, append=TRUE)


##Give warning if user-supplied externals are not in model
if (!all(externals %in% all.components)){
warning("Warning! Following metabolites provided in argument 'externals' are not in the metabolites of given model: ",
	paste(externals[which(!externals %in% all.components)], collapse=","), "\n")
}

##EXTERNALS
if (length(externals) > 0){
write("###   EXTERNALS", file=file.name, append=TRUE)
write(externals, file=file.name, append=TRUE)   #User specified externals for a system
write("###   END EXTERNALS", file=file.name, append=TRUE) 
}
message("Lim file written to ", paste(file.name, sep="/"), "\n")
}
