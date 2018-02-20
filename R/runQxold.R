### runs a simple Qx analysis (no conditional analysis) using genotype information (# old allele copies not frequencies), a kinship matrix, and a list of effect sizes.


function(myGenos, myK, myBetas){  #myK is kinship just for genotyping set, should be already centered.
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

