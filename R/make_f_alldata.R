
##makes a kinship matrix when there is no missing data. Takes frequency data, not # of copies. Later would be good to add a check for this. This function also drops an individual to deal with losing a degree of freedom when mean centering. The input table has rows as individuals, loci as columns.

make_f_alldata <- function(myG){
  myEs = 1/(colMeans(myG)*(1-colMeans(myG)))
  myEsr = replace(myEs, myEs==Inf, 0)
  myS = matrix(0,nrow=dim(myG)[2], ncol=dim(myG)[2])
  diag(myS) = myEsr
  myM = dim(myG)[1]
  myT = matrix(data = -1/myM, nrow=myM-1, ncol=myM)
  diag(myT) = (myM-1)/myM
  myK = dim(myG)[2]
  myF = (1/(myK-1)) * myT %*% myG %*% myS %*% t(myG) %*% t(myT)
  return(myF)
}


