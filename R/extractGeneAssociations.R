##extract gene information from database
extractGeneAssociations <- function(database){
  lapply(database@model@reactions, function(r){
    gene.ass <- NA
    if (length(r@notes)>0){
      gene.ass <- gsub(".*GENE.ASSOCIATION:(.*?)</.*>.*$", "\\1", r@notes)
      gene.ass <- ifelse(gene.ass==r@notes, NA, gene.ass)
    }
    gene.ass
  })
}
