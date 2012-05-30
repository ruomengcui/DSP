solve.for.A<-function(p){
#Find fixed-point-solution impulse response matrix, given exogenous parameters p={alpha, which.alpha, beta, which.beta, Sigma, Lambda, Gamma, L}. Items which.alpha and which.beta tell which parameters are non-zero; if blank we assume it's the first ones that are nonzero. The routine searches across matrix A.tilde, which has dimentions (H-L)*(H+1).
  if(is.null(p$which.alpha)) p$which.alpha=sq(p$alpha)-1  
  if(is.null(p$which.beta)) p$which.beta=sq(p$beta)-1
	attach(p)
	H=dim(Sigma)[1]-1
	G=makeG(H, L)
	K=makeK(H, L)	
	fixedPointError<-function(A.tilde){
		A=K+G%*%matrix(A.tilde, nrow=H-L)
		c((A-iterate.A(A))[seq(H-L),])
	}
	A.tilde.0=.9*makeR(H+1, L)[seq(H-L),]
	o=BBsolve(c(A.tilde.0), fixedPointError, quiet=TRUE)
	o$A=K+G%*%matrix(o$par, nrow=H-L)
	detach(p)
	o
}

iterate.A<-function(A){
	get.A.from.marginals(get.marginals.from.A(A))
}

get.A.from.marginals<-function(a.b){
#a.b is a list of two marginal cost vectors. The elements corrospond to which.alpha and which.beta entries.
	H=dim(Sigma)[1]-1
	C=makeC(H+1)
	D=makeD(H+1, L)
	G=makeG(H, L)
	K=makeK(H, L)
	a.terms=laply(sq(a.b$a), function(n){a.b$a[n]*t(make.S(H, which.alpha[n]))%*%make.S(H, which.alpha[n])})
	b.terms=laply(sq(a.b$b), function(n){a.b$b[n]*t(C)%*%makeIj(H+1,which.beta[n]+1)%*%C})
	if(length(dim(a.terms))==3) a.terms=apply(a.terms, c(2,3), sum)
	if(length(dim(b.terms))==3) b.terms=apply(b.terms, c(2,3), sum)
	
	g=t(D)%*%t(C)%*%C%*%D+a.terms+b.terms
	f=K+G%*%solve(t(G)%*%g%*%G, t(G)%*%t(D)%*%t(C)%*%C-t(G)%*%g%*%(K+Gamma%*%solve(Sigma)))
	f
}

get.marginals.from.A<-function(A){
#calc marginal costs of increasing various standard deviations.
	s=getStdDiv(A)	
	a=laply(sq(alpha), function(x){alpha[x]*s$std.i/s$std.o[x]})
	b=laply(sq(beta), function(x){beta[x]*s$std.i/s$std.u[x]})
	llply(list(a=a, b=b), pmin, 10^10)
}

getStdDiv<-function(A){
	H=dim(Sigma)[1]-1
	C=makeC(H+1)
	D=makeD(H+1, L)	
	std.i=sqrt(sum(diag(C%*%(D%*%A-diag(H+1))%*%Sigma%*%t(D%*%A-diag(H+1))%*%t(C)+2*C%*%(D%*%A-diag(H+1))%*%Gamma%*%t(D)%*%t(C)+C%*%D%*%Lambda%*%t(D)%*%t(C))))
	std.o=laply(sq(alpha), function(j){
		Sn=make.S(H, which.alpha[1])
		sqrt(sum(diag(Sn%*%A%*%Sigma%*%t(A)%*%t(Sn)+2*Sn%*%A%*%Gamma%*%t(Sn)+Sn%*%Lambda%*%t(Sn))))
	})
	std.u=laply(sq(beta), function(j){
		D.l=makeD(H+1, which.beta[j]+1)
    R.l=makeR(H+1, which.beta[j]+1)
    sqrt(sum(diag(C%*%(D.l%*%R.l-diag(H+1))%*%(A%*%Sigma%*%t(A)+A%*%Gamma+t(Gamma)%*%t(A)+Lambda)%*%t(D.l%*%R.l-diag(H+1))%*%t(C))))
	})
	std=list(std.i=std.i, std.o=std.o, std.u=std.u)
	llply(std, pmax, 10^-16)
}