
##makes a kinship matrix using the cov function. input is a matrix where the rows are individuals/populations and the columns are loci and the values are the allele frequency.

make_k_complete <- function(myG){
scaleFactor = 1/sqrt(colMeans(myG) * (1 - colMeans(myG)))
myS = matrix(0, nrow = dim(myG)[2], ncol = dim(myG)[2])
diag(myS) = scaleFactor
myM = dim(myG)[1]
myT = matrix(data = -1/myM, nrow = myM - 1, ncol = myM)
diag(myT) = (myM - 1)/myM
myGstand = myT %*% myG %*% myS
myK = cov(t(myGstand))
return(myK)
}

