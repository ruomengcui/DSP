library('numDeriv')
library('BB')

solveFixedPoint<-function(alpha.beta, S, Lambda, L.Ls=c(0,1)){
#Find fixed-point-solution impulse response matrix, given exogenous parameters, alpha, beta, Sigma, and Lambda. Routine searches across matrix A.
	H=sqrt(length(S))
	fixedPointError<-function(A){
		A=matrix(A, sqrt(length(A)))
		c(A-iterate.A(A, S, Lambda, alpha.beta, L.Ls))
	}
	o=BBsolve(c(diag(H))/2, fixedPointError, quiet=TRUE)
	o$A=matrix(o$par,sqrt(length(o$par)))
	o
}

iterate.A<-function(A, S, Lambda, alpha.beta, L.Ls){
	get.A.from.marginals(sqrt(length(S)), get.marginals.from.A(A, S, Lambda, alpha.beta, L.Ls), L.Ls)
}

get.A.from.marginals<-function(H, a.b, L.Ls=c(0,1)){
#a.b are the relative marginal costs of increasing order stdiv and supplier inv std div.
	C=makeC(H)
	D=makeD(H, L.Ls[1])
	G=makeG(H, L.Ls[1])
	B=makeB(H, L.Ls[1])
	g=pmin(t(D)%*%t(C)%*%C%*%D+a.b[1]*diag(H)+a.b[2]*t(C)%*%makeIj(L.Ls[2],H)%*%C, 10^12)
	f=B+G%*%(solve(t(G)%*%g%*%G, t(G)%*%t(D)%*%t(C)%*%C-t(G)%*%g%*%B))
	f
}

get.marginals.from.A<-function(A, S, Lambda, alpha.beta, L.Ls){
#calc marginal costs of increasing order stdiv and supplier inv std div.
	iota=c(1,alpha.beta)/getStdDiv(A, S, Lambda, L.Ls)
	iota[2:3]/iota[1]
}

getStdDiv<-function(A, S, Lambda, L.Ls){
#the first element is the standard deviation of inventory, the second is the standard deviation of production, and the third the standard deviation of supplier inventories.
	n=dim(A)[1]
	C=makeC(n)
	D=makeD(n, L.Ls[1])
	std=array(dim=3)
	std[1]=sqrt(sum(diag(C%*%(D%*%A-diag(n))%*%S%*%t(D%*%A-diag(n))%*%t(C)+C%*%D%*%Lambda%*%t(D)%*%t(C))))
	std[2]=sqrt(sum(diag(A%*%S%*%t(A)+Lambda)))
	Ij=makeIj(L.Ls[2],n)
	std[3]=sqrt(sum(diag(Ij%*%C%*%(A%*%S%*%t(A)+Lambda)%*%t(C)%*%Ij)))
	std
}

makeC <- function(n){
	C=matrix(0, n, n)
	C[lower.tri(C, diag = TRUE)] <-1
	C
}

makeD <- function(n, L){
	D=matrix(0, n, n)
	D[seq(L+1,n^2-n*L,n+1)]<-1
	D
}

makeB<-function(H, L){
	o=matrix(0, H, H)
	o[H-L,]<-1
	o
}

makeG<-function(H, L){
	rbind(diag(H-L-1), rep(-1, H-L-1), matrix(0, L, H-L-1))
}

makeIj <- function(j, n){diag(c(array(1,j),array(0,n-j)))}

makeR <- function(n, j){R=t(makeD(n, j)); R[1,seq(j)]=1; R}

make.el<-function(n, l){ x=rep(0, n); x[l]=1; x}