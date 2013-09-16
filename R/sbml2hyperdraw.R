sbml2hyperdraw <- function(sbml.model,
                           rates=NULL,
                           relevant.species=unname(sapply(sbml.model@species, id)),
                           relevant.reactions=unname(sapply(sbml.model@reactions, id)),
                           layoutType="dot", lwd.max=3, lwd.min=0.5, plt.margins=c(150, 150, 150, 150)){
  

  
  ##Build hyperedges from  SBML reactions
  hyperedges <- unlist(sapply(sbml.model@reactions, function(r){
    if(r@id %in% relevant.reactions){
      reactants <- intersect(sapply(r@reactants, species), relevant.species)
      products <- intersect(sapply(r@products, species), relevant.species)
      if (length(reactants)>0 & length(products)>0){
        ##Change direction into direction of flux if rate is  < 0
        my.label <- ifelse(r@id %in% names(rates), paste(r@id, round(rates[r@id], 2), sep=":"), r@id)
        if (! is.null(rates) && rates[r@id]<0)
          DirectedHyperedge(products, reactants, label=my.label)
        else      
          DirectedHyperedge(reactants, products, label=my.label)
      }
    }
  }, simplify=T))


  ##If rates are given, append them  to edge names
  if(hasArg(rates)){
    names(hyperedges) <- paste(names(hyperedges), ":", round(rates, 2)[names(hyperedges)], sep="")
    names(rates) <- paste(names(rates), ":", round(rates, 2), sep="")
  } else {
    rates <- sapply(sbml.model@reactions, function(r)1) ##if no rates given, set them to one
  }
  
  ##Build graph object
  node.names <- unique(unlist(c(lapply(hyperedges, function(x)x@head), lapply(hyperedges, function(x)x@tail))))
  hg <- Hypergraph(node.names, hyperedges)
  testbph <- graphBPH(hg)
  my.graph <- graphLayout(testbph, layoutType=layoutType)
  
  nodeDataDefaults(my.graph, "shape") <- "box"
  nodeDataDefaults(my.graph, "margin") <- 'unit(3, "mm")'  
  edgeDataDefaults(my.graph, "lwd") <- 1
  graphDataDefaults(my.graph, "arrowLoc") <- "end"
    
  ##calculate line widths from rates
  lwds <- abs(rates[my.graph@edgeNodes]) / max(abs(rates[my.graph@edgeNodes])) * lwd.max
  lwds[which(lwds<lwd.min)] <- lwd.min
  for (rxn.id in names(rates)){
    ##set grey color for reaction edges with rate 0     
    if (unname(rates[rxn.id])==0){
      lapply(my.graph@edgeNodeIO$outgoing[[rxn.id]], function(x)edgeData(my.graph, rxn.id, x, "color") <- "grey")
      lapply(my.graph@edgeNodeIO$incoming[[rxn.id]], function(x)edgeData(my.graph, x, rxn.id, "color") <- "grey")    
    }
    ##set red color for reaction edges with rate < 0 
    ##also reverse direction of arrow for reactions edges with rate < 0
    if (unname(rates[rxn.id])<0){
      lapply(my.graph@edgeNodeIO$outgoing[[rxn.id]], function(x)edgeData(my.graph, rxn.id, x, "color") <- "red")
      lapply(my.graph@edgeNodeIO$incoming[[rxn.id]], function(x)edgeData(my.graph, x, rxn.id, "color") <- "red")    
    }
    
    ##Set line widths
    lwd <- unname(ifelse(lwds[rxn.id]==2.8,  as.character(min(lwds)), as.character(lwds[rxn.id]))) ##set lwd to min lwd if rate is zero
    lapply(my.graph@edgeNodeIO$outgoing[[rxn.id]], function(x)edgeData(my.graph, rxn.id, x, "lwd") <- as.character(lwds[rxn.id]))
    lapply(my.graph@edgeNodeIO$incoming[[rxn.id]], function(x)edgeData(my.graph, x, rxn.id, "lwd") <- as.character(lwds[rxn.id]))
  }
  
  ##Set plot margins
  my.graph@graph@boundBox@botLeft@y <- my.graph@graph@boundBox@botLeft@y-plt.margins[1] ##bottom
  my.graph@graph@boundBox@botLeft@x <- my.graph@graph@boundBox@botLeft@x-plt.margins[2] ##left 
  my.graph@graph@boundBox@upRight@y <- my.graph@graph@boundBox@upRight@y+plt.margins[3] ##top
  my.graph@graph@boundBox@upRight@x <- my.graph@graph@boundBox@upRight@x+plt.margins[4] ##right  
  ##  pushViewport(viewport(gp=gpar(cex=cex.lab)))

  return(my.graph)
}
