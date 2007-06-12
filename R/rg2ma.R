"rg2ma" <-
function(r,g){
  x <- log2(cbind(r,g))
  x <- x %*% cbind(M=c(1,-1), A=c(1,1))
  x[,2] <- x[,2]/2 
  return(x)
}

