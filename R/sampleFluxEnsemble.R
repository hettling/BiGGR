sampleFluxEnsemble <- function(model, uncertain.vars=NULL, iter=3000, ...){
lim <- NULL
if (class(model) == "character"){
lim <- Setup(model)
} else if (class(model) == "lim"){
lim <- model
} else if (class(model) == "Model"){
limfile.path <- tempfile()
createLIMFromSBML(model, file.name=limfile.path)
lim <- Setup(limfile.path)
}
if (is.null(uncertain.vars)){
Xsample(lim, ...)
} else {
A <- t(apply(uncertain.vars, 1, function(x){
  str <- x[1]
  splitted <- strsplit(str, " ")[[1]]
  if (length(splitted) %% 2) splitted <- c('+', splitted)
  signs <- splitted[c(TRUE, FALSE)]
  vars  <- splitted[c(FALSE, TRUE)]
  a.vec <- rep(0, length(lim$Unknowns))
  names(a.vec) <- lim$Unknowns
  a.vec[vars] <- ifelse(signs == "+", -1, 1)
  a.vec
}))
B <- -uncertain.vars[,2]
sdB <- uncertain.vars[,3]    
##remove equality constraints that are approximate equalities in A and B from  E and F
`%inm%` <- function(x, matrix){
  test <- apply(matrix, 1, `==`, x)
  any(apply(test, 2, all))
}
remove.idx <- which(apply(lim$A, 1, '%inm%', A))
E <- lim$A
F <- lim$B
if (length(remove.idx)>0){
  E <- E[-remove.idx,]
  F <- F[-remove.idx]
}  
xsample(A=A, B=B, E=E, F=F, G=lim$G, H=lim$H, sdB=sdB, iter=iter, ...)$X    
}
}
