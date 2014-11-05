rmvSpliceVariant<-function(gene.info){
gene.info$GPR<-gsub(" ","",gene.info$GPR,fixed=TRUE)
gene.info$GPR<-gsub("(","",gene.info$GPR,fixed=TRUE)
gene.info$GPR<-gsub(")","",gene.info$GPR,fixed=TRUE)
gene.info$GPR<-gsub("\\.[0-9]+","",gene.info$GPR)
return(gene.info)
}