#
# Example 1
#

set.seed(100)

x<-matrix(rnorm(50*20),ncol=20)
s<-var(x)
a<-glasso(s, rho=.01)
aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)

# example with structural zeros and no regularization,
# from Whittaker's Graphical models book  page xxx.
s=c(10,1,5,4,10,2,6,10,3,10)
S=matrix(0,nrow=4,ncol=4)
S[row(S)>=col(S)]=s
S=(S+t(S))
diag(S)<-10
zero<-matrix(c(1,3,2,4),ncol=2,byrow=TRUE)
a<-glasso(S,0,zero=zero)

plot.matrix.fn <- function(m,maxval) {
  #enter answer to this question here
  image(abs(t(m[nrow(m):1,] )), axes=FALSE, zlim=c(0,maxval))
}

par(mfrow = c(3,2))
plot.matrix.fn(S,10)
title("true covariance")
plot.matrix.fn(a$w,10)
title("estimated covariance")
plot.matrix.fn(solve(S),max(solve(S)))
title("inverse covariance")
plot.matrix.fn(a$wi,max(solve(S)))
title("estimated inverse")
plot.matrix.fn(abs(solve(S))<1e-10,1)
title("inverse support")
plot.matrix.fn(abs(a$wi)<1e-10,1)
title("estimated support")

#----------------------------------------------------------------------------
#
# Example 2
#


data(HumanPw)

#g<-graph_from_adjacency_matrix(dataHuman$TDag)
g1<-graph.adjacency(dataHuman$TDag)
g1$layout <- layout_in_circle

par(mfrow = c(1,2))
plot.matrix.fn(dataHuman$TDag,1)
title("adjacency matrix")
plot(g1)
title("dag")

