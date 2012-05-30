library('numDeriv')
library('BB')

getA<-function(T, alpha, beta, L=0, Ls=1){
	#the first element in gamma is gamma_p, and the remaining are gamma_j
	C=makeC(T)
	D=makeD(T, L)
	G=makeG(T, L)
	B=makeB(T, L)
	g=pmin(t(G)%*%(t(D)%*%t(C)%*%C%*%D+alpha*diag(T)+beta*t(C)%*%makeIj(Ls,T)%*%C), 10^12)
	f=B+G%*%(solve(g%*%G, t(G)%*%t(D)%*%t(C)%*%C-g%*%B))
	f
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

makeB<-function(T, L){
	o=matrix(0, T, T)
	o[T-L,]<-1
	o
}

makeG<-function(T, L){
	rbind(diag(T-L-1), rep(-1, T-L-1), matrix(0, L, T-L-1))
}

makeIj <- function(j, n){diag(c(array(1,j),array(0,n-j)))}

makeR <- function(n, j){R=t(makeD(n, j)); R[1,seq(j)]=1; R}

make.el<-function(n, l){ x=rep(0, n); x[l]=1; x}

getGamma<-function(A, S, Lambda, k, L){
#the first element is gamma_i, the second is gamma_p, and the remaining are gamma_c^j.
	iota = c(1,k)/getStdDiv(A, S, Lambda, L)
	iota[2:length(iota)]/iota[1]
}

getStdDiv<-function(A, S, Lambda, L){
#the first element is the standard deviation of inventory, the second is the standard deviation of production, and the remaining are the standard deviation of the supplier's inventories.
	n=dim(A)[1]
	C=makeC(n)
	D=makeD(n, L)
	std = array(0,n+1)
	std[1] = sqrt(sum(diag(C%*%(D%*%A-diag(n))%*%S%*%t(D%*%A-diag(n))%*%t(C)+C%*%D%*%Lambda%*%t(D)%*%t(C))))
	std[2] =sqrt(sum(diag(A%*%S%*%t(A)+Lambda)))
	for(j in 1:n){
		Ij=makeIj(j,n)
		std[j+2] =sqrt(sum(diag(Ij%*%C%*%(A%*%S%*%t(A)+Lambda)%*%t(C)%*%Ij)))
	} 
	std
}

calcObjective<-function(A, S, Lambda, k, L){	
#k_i is normalized to 1; the first element of k is k_p, and the remaining are theta^j*k_c^{j-1}.
	A=matrix(A, sqrt(length(A)), sqrt(length(A)))
	sum(c(1,k)*getStdDiv(A, S, Lambda, L))	
}

rawOptimization<-function(k, S=diag(length(k)-1), Lambda=diag(length(k)-1), L, startA=c(diag(length(k)-1))/2){
	o=optim(startA, calcObjective, gr = NULL, S, k, L, method = "BFGS", control = list(maxit=5000))
	o$A=matrix(o$par,sqrt(length(o$par)),sqrt(length(o$par)))
	o
}

solveFixedPoint<-function(k, S=diag(length(k)-1), Lambda=diag(length(k)-1), L=0){
	fixedPointError<-function(A, S, k, L){
		A=matrix(A, sqrt(length(A)), sqrt(length(A)))
		c(A-getA(getGamma(A, S, Lambda, k, L), L))
	}
	o=BBsolve(c(diag(length(k)-1))/2, fixedPointError, method=c(2,3,1), control=list(), quiet=FALSE, S, k, L)
	o$A=matrix(o$par,sqrt(length(o$par)),sqrt(length(o$par)))
	o
}

checkFixedPoint<-function(k, S=diag(length(k)-1), Lambda=diag(length(k)-1), L=0){
#first solve fixed point, then optimize from that fixed point to see if I can do any better
	o1=solveFixedPoint(k, S, Lambda, L)
	o2=rawOptimization(k, S, Lambda, L, o1$par)
	list(o1$A, sum(o1$A-o2$A))
}
