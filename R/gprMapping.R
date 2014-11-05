
gprMapping<-function(gene_express,react_gene_map,OR=c("mean","median","min","max"),AND=c("min","max","mean","median")){
match.arg(OR)
match.arg(AND)
OR<-match.fun(OR)
AND<-match.fun(AND)

gpr.map<-NULL
for (j in 1:dim(react_gene_map)[1]){
al<-str_detect(react_gene_map[j,2],"and")
ol<-str_detect(react_gene_map[j,2],"or")
# This is a check for OR !!!
if(al==FALSE){
  t<-unlist(strsplit(react_gene_map[j,2],split="or",fixed=TRUE))
  dt<-NULL
  dmm<-NULL
 for (k in 1:length(t)){
   dt1<-gene_express[which(t[k]==gene_express[,2]),3]
   if (length(dt1) == 0) dt1<-0  # to remove the 0 values
	dt<-cbind(dt,dt1)
 }
 dmm<-paste(dt,collapse="or")
 #dt<-cbind(react_gene_map[j,1:2],dmm,mean(dt))
  dt<-cbind(react_gene_map[j,1:2],mean(dt))
  colnames(dt)<-c("react_id","formula","GPR")
 gpr.map<-rbind(gpr.map,dt)
 # Remember this is a check for AND !!!
}else if(ol==FALSE){
 t<-unlist(strsplit(react_gene_map[j,2],split="and",fixed=TRUE))
  dt<-NULL
 dmm<-NULL
  for (k in 1:length(t)){
   dt1<-gene_express[which(t[k]==gene_express[,2]),3]
	dt<-cbind(dt,dt1)
  }
 dmm<-paste(dt,collapse="and")
  #dt<-cbind(react_gene_map[j,1:2],dmm,min(dt))
  dt<-cbind(react_gene_map[j,1:2],min(dt))
 colnames(dt)<-c("react_id","formula","GPR")
  gpr.map<-rbind(gpr.map,dt)
 # Remember this is a check for both AND and OR !!!
}else if (ol==TRUE & al==TRUE) {
  t<-unlist(strsplit(react_gene_map[j,2],split="or",fixed=TRUE))
  dt<-NULL
  dmm<-NULL
 dt1<-NULL
  for (k in 1:length(t)){
	u<-unlist(strsplit(t[k],split="and",fixed=TRUE))
	dt1<-cbind(dt1,min(na.omit(gene_express[match(u,gene_express[,2]),3])))
 }
 #dt<-cbind(react_gene_map[j,1:2],"This is not important",mean(dt1))
 dt<-cbind(react_gene_map[j,1:2],mean(dt1))
  colnames(dt)<-c("react_id","formula","GPR")
  gpr.map<-rbind(gpr.map,dt) 
}
}
return(gpr.map)
#write.csv(file1,file="Lewis_Reaction_GPR_complete.csv")
}
















