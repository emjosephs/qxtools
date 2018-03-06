#runs Qx on a set of individuals or popualtions. Inouts are the genotype table where rows are individuals/populations and columns are loci and the values are the allele frequency within individuals/populations. The second input in the kinship matrix for the individuals inputed (one individual has been dropped due to centering). The third input is the effect sizes of the loci described in the genotyping table.

run_qx_basic <- function(myGenos, myK, myBetas){  #myK is kinship just for genotyping set, should be already centered.
  myM = dim(myK)[1]
  
  #get breeding values
  myZ = myGenos %*% myBetas
  
  #center breeding values
  myT = matrix(data = rep(-1/(myM+1), (myM+1)*(myM)), nrow=myM, ncol=myM+1)
  diag(myT) = (myM)/(myM+1)
  myZprime = myT %*% myZ
  
  #calculate Va
  myVa = calcVa(colMeans(myGenos)/2 ,myBetas)
  
  #actually calculate things
  qxOld = t(myZprime) %*% solve(2*myVa*myK) %*% (myZprime)
  
  return(list(qxOld, myZprime, myVa))
}

