getRates<-function(modelFile){
LP<-Linp(modelFile)
unlist(as.data.frame(LP$X))   #optimized rates from the model file
}
