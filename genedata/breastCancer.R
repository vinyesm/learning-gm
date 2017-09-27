## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("breastCancerNKI")
biocLite("Biobase")
biocLite("impute")
biocLite("preprocessCore")
# 
## load Biobase package
library(breastCancerNKI)
library(Biobase)
library(impute)
library(preprocessCore)
## load the dataset
data(nki)
## show the first 5 rows and columns of the expression data
exprs(nki)[1:5,1:5]
## show the first 6 rows of the phenotype data
head(pData(nki))
## show first 20 feature names
featureNames(nki)[1:20]
## show the experiment data summary
experimentData(nki)
## show the used platform
annotation(nki)
## show the abstract for this dataset
abstract(nki)

## keep genes with less than 10% missing values
N = dim(exprs(nki))[1]
P = dim(exprs(nki))[2]
ex = exprs(nki)
filter = 1.*is.na(ex)
filter = apply(filter,1,sum)
filter = filter<.01*P
exprs2 = ex[filter,]
dim(exprs2)
sum(1.*is.na(exprs2))

## impute missing values
nki.imputed<-impute.knn(exprs2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
sum(1.*is.na(nki.imputed))

## quantile normalization
#nki.qn <- normalize.quantiles(nki.imputed,copy=TRUE)

## Save in Matlab v6 format with 'writeMat'
library(R.matlab)
writeMat("Rosetta.mat", exprs = nki.imputed, names = featureNames(nki.imputed))


# to intall rJava Package
# in console:
#/usr/libexec/java_home -v 1.8
#/Library/Java/JavaVirtualMachines/jdk1.8.0_121.jdk/Contents/Home
# in R:
#options("java.home"="/Library/Java/JavaVirtualMachines/jdk1.8.0_121.jdk/Contents/Home/jre")
#options('java.home')
# Sys.setenv(JAVA_HOME='/Library/Java/JavaVirtualMachines/jdk1.8.0_121.jdk/Contents/Home')
# Sys.setenv(LD_LIBRARY_PATH='$JAVA_HOME/jre/lib/server')
# Sys.getenv()
# library(xlsx)
# write.xlsx(nki, "Rosetta.xlsx")
