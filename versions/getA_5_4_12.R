library('numDeriv')
library('BB')

solveFixedPoint<-function(alpha.beta, S, Lambda, L.Ls=c(0,1)){
#Find fixed-point-solution impulse response matrix, given exogenous parameters, alpha, beta, Sigma, and Lambda. Routine searches across matrix A.tilde, which has dimentions (H-L)*(H+1).
	H=dim(S)[1]-1
	G=makeG(H, L.Ls[1])
	K=makeK(H, L.Ls[1])
	fixedPointError<-function(A.tilde){
		A=K+G%*%matrix(A.tilde, nrow=H-L.Ls[1])
		c((A-iterate.A(A, S, Lambda, alpha.beta, L.Ls))[seq(H-L.Ls[1]),])
	}
	A.tilde.0=makeR(H+1, L.Ls[1])[seq(H-L.Ls[1]),]
	o=BBsolve(c(A.tilde.0), fixedPointError, quiet=TRUE)
	o$A=K+G%*%matrix(o$par, nrow=H-L.Ls[1])
	o
}

iterate.A<-function(A, S, Lambda, alpha.beta, L.Ls=c(0,1)){
	H=dim(S)[1]-1
	get.A.from.marginals(H, get.marginals.from.A(A, S, Lambda, alpha.beta, L.Ls), L.Ls)
}

get.A.from.marginals<-function(H, a.b, L.Ls=c(0,1)){
#a.b are the relative marginal costs of increasing order stdiv and supplier inv std div.
	C=makeC(H+1)
	D=makeD(H+1, L.Ls[1])
	G=makeG(H, L.Ls[1])
	K=makeK(H, L.Ls[1])
	g=t(D)%*%t(C)%*%C%*%D+a.b[1]*diag(H+1)+a.b[2]*t(C)%*%makeIj(H+1,L.Ls[2])%*%C
	g=pmin(g, 10^12)
	f=K+G%*%(solve(t(G)%*%g%*%G, t(G)%*%t(D)%*%t(C)%*%C-t(G)%*%g%*%K))
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
	Ij=makeIj(n,L.Ls[2])
	std[3]=sqrt(sum(diag(Ij%*%C%*%(A%*%S%*%t(A)+Lambda)%*%t(C)%*%Ij)))
	pmax(std, 10^-20)
}

makeC <- function(H){
	C=matrix(0, H, H)
	C[lower.tri(C, diag = TRUE)] <-1
	C
}

makeD <- function(H, L){
	x=rbind(matrix(0, abs(L), H), diag(1, H-abs(L), H))
	if(L>=0) x
	else t(x)
}

makeK<-function(H, L){
	o=matrix(0, H+1, H+1)
	o[H+1-L,]<-1
	o
}

makeG<-function(H, L){
	rbind(diag(H-L), rep(-1, H-L), matrix(0, L, H-L))
}

makeIj <- function(H, L){diag(c(array(1,L),array(0,H-L)))}

makeR <- function(H, L){R=t(makeD(H, L)); R[1,seq(L)]=1; R}

make.el<-function(H, l){ x=rep(0, H); x[l]=1; x}


#######Sound#########
get.A.from.marginals.sound<-function(H, theta, L=0){
	C=makeC(H+1)
	D=makeD(H+1, L)
	G=makeG(H, L)
	K=makeK(H, L)

	sum.theta=aaply(laply(seq(theta), function(j){theta[j]*t(make.Delta.j(H, j-1))%*%(make.Delta.j(H, j-1))}), 2:3, sum)
	K+G%*%solve(t(G)%*%t(D)%*%t(C)%*%C%*%D%*%G+t(G)%*%sum.theta%*%G, t(G)%*%t(D)%*%t(C)%*%C%*%(diag(H+1)-t(D)%*%K)-t(G)%*%sum.theta%*%K)
}

make.Delta.j<-function(H, j){
	((diag(H+1+j)-makeD(H+1+j,1))%^%j)[,seq(H+1)]
}


########Convolution#########
conv<-function(x1, x2){ 
	vapply(seq(x1)-ceiling(length(x1)/2), function(m){ x1%*%makeD(length(x2), m)%*%x2[length(x2):1]}, 0)
}