##extract pathay information from database
extractPathways <- function(database){
lapply(database@model@reactions, function(r){
pathway <- NA
if (length(r@notes)>0){
  pathway <- gsub(".*SUBSYSTEM: (.*?)</.*>.*$", "\\1", r@notes)
  pathway <- ifelse(pathway==r@notes, NA, pathway)
}
pathway
})
}
