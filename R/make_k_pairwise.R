
##makes a kinship matrix using pairwise complete data. input is a matrix where the rows are loci and the columns are individuals or populations and the values are the number of allele copies in each individual (ie 2x allele frequency).

make_k_pairwise <- function(myMatrix){
mygenos = as.matrix(myMatrix)/2 #divide to get frequencies
scaleFactor = 1/sqrt(rowMeans(mygenos, na.rm=T)*(1-rowMeans(mygenos, na.rm=T)))  #scale (rows are locus so divde each by sqrt or mean(1-mean))
myGcent = scale(mygenos, scale=F)
myGstand =  t(sapply(1:nrow(myGcent), function(x) {myGcent[x,]/scaleFactor[x]}))
myK = cov(myGstand, use="pairwise.complete.obs")
return(myK)
}

