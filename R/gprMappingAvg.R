
gprMappingAvg<-function(gene_express,react_gene_map){
file1<-NULL

for (i in 1:dim(react_gene_map)[1]){
t<-as.integer(unlist(strsplit(react_gene_map[i,2],"[andor]")))
t<-as.vector(na.omit(t))
dt<-NULL
dmm<-NULL
for (k in 1:length(t)){
dt1<-gene_express[which(t[k]==gene_express[,2]),3]
dt<-cbind(dt,dt1)
}
dt<-cbind(react_gene_map[i,1:2],mean(dt))
colnames(dt)<-c("react_id","GPR","average") 
file1<-rbind(file1,dt)
ff<-file1

}
return(ff)
#write.csv(file1,file="Lewis_Reaction_GPR_complete.csv")
}














