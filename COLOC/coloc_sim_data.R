library(coloc)
setClass("simdata", representation(df1="data.frame",df2="data.frame"))
setValidity("simdata", function(object) {
  n <- nrow(object@df1)
  if(nrow(object@df2)!=n)
    return("nrow of '@df1' should equal nrow of '@df2'")
})

setMethod("show", signature="simdata", function(object) {
  cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n")
})

sim.data <- function(nsnps=50,nsamples=200,causals=1:2,nsim=1) {
  cat("Generate",nsim,"small sets of data\n")
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y2 <- rnorm(nsamples,rowSums(X2[,causals]),2)
  colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
  df1 <- cbind(Y=Y1,X1)
  df2 <- cbind(Y=Y2,X2)
  if(nsim==1) {
    return(new("simdata",
               df1=as.data.frame(df1),
               df2=as.data.frame(df2)))
  } else {
    index <- split(1:(nsamples * nsim), rep(1:nsim, nsamples))
    objects <- lapply(index, function(i) new("simdata", df1=as.data.frame(df1[i,]),
                                             df2=as.data.frame(df2[i,])))
    return(objects)
  }
}

## simulate some data and load the coloc library
set.seed(46411)
data <- sim.data(nsamples=1000,nsim=1)

#df <- as.data.frame(cbind(Y=Y, as(X,"numeric")))

if(requireNamespace("snpStats")) {
  library(snpStats)
  
  Y1 <- data@df1$Y
  Y2 <- data@df2$Y
  
  X1 <- new("SnpMatrix",as.matrix(data@df1[,-1]))
  X2 <- new("SnpMatrix",as.matrix(data@df2[,-1]))
  
  p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
  p2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
  
  maf <- col.summary(X2)[,"MAF"]
  
  get.beta <- function(...) {
    tmp <- snpStats::snp.rhs.estimates(..., family="gaussian")
    beta <- sapply(tmp,"[[","beta")
    varbeta <- sapply(tmp, "[[", "Var.beta")
    return(list(beta=beta,varbeta=varbeta))
  }
  b1 <- get.beta(Y1 ~ 1, snp.data=X1)
  b2 <- get.beta(Y2 ~ 1, snp.data=X2)
} else {
  Y1 <- data@df1$Y
  Y2 <- data@df2$Y
  
  X1 <- as.matrix(data@df1[,-1])
  X2 <- as.matrix(data@df2[,-1])
  
  tests1 <- lapply(1:ncol(X1), function(i) summary(lm(Y1 ~ X1[,i]))$coefficients[2,])
  tests2 <- lapply(1:ncol(X2), function(i) summary(lm(Y2 ~ X2[,i]))$coefficients[2,])
  
  p1 <- sapply(tests1,"[",4)
  p2 <- sapply(tests2,"[",4)
  
  maf <- colMeans(X2)/2
  
  get.beta <- function(x) {
    beta <- sapply(x,"[",1)
    varbeta <- sapply(x, "[", 2)^2
    return(list(beta=beta,varbeta=varbeta))
  }
  b1 <- get.beta(tests1)
  b2 <- get.beta(tests2)
} 

my.res <- finemap.abf(dataset=list(beta=b1$beta, varbeta=b1$varbeta, N=nrow(X1),sdY=sd(Y1),type="quant"))
head(my.res)
tail(my.res) 
if(requireNamespace("snpStats")) {
  library(snpStats)
  Y1 <- data@df1$Y
  Y2 <- data@df2$Y
  
  X1 <- new("SnpMatrix",as.matrix(data@df1[,-1]))
  X2 <- new("SnpMatrix",as.matrix(data@df2[,-1]))
  
  p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
  p2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
  
  maf <- col.summary(X2)[,"MAF"]
  
  get.beta <- function(...) {
    tmp <- snpStats::snp.rhs.estimates(..., family="gaussian")
    beta <- sapply(tmp,"[[","beta")
    varbeta <- sapply(tmp, "[[", "Var.beta")
    return(list(beta=beta,varbeta=varbeta))
  }
  b1 <- get.beta(Y1 ~ 1, snp.data=X1)
  b2 <- get.beta(Y2 ~ 1, snp.data=X2)
} else {
  Y1 <- data@df1$Y
  Y2 <- data@df2$Y
  
  X1 <- as.matrix(data@df1[,-1])
  X2 <- as.matrix(data@df2[,-1])
  
  tests1 <- lapply(1:ncol(X1), function(i) summary(lm(Y1 ~ X1[,i]))$coefficients[2,])
  tests2 <- lapply(1:ncol(X2), function(i) summary(lm(Y2 ~ X2[,i]))$coefficients[2,])
  
  p1 <- sapply(tests1,"[",4)
  p2 <- sapply(tests2,"[",4)
  
  maf <- colMeans(X2)/2
  
  get.beta <- function(x) {
    beta <- sapply(x,"[",1)
    varbeta <- sapply(x, "[", 2)^2
    return(list(beta=beta,varbeta=varbeta))
  }
  b1 <- get.beta(tests1)
  b2 <- get.beta(tests2)
}
my.res <- coloc.abf(dataset1=list(beta=b1$beta, varbeta=b1$varbeta, N=nrow(X1),sdY=sd(Y1),type="quant"),
                    dataset2=list(beta=b2$beta, varbeta=b2$varbeta, N=nrow(X2),sdY=sd(Y2),type="quant"),
                    MAF=maf)
print(my.res[[1]]) 
