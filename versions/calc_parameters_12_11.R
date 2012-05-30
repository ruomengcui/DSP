calc.firm<-function(id){
	H=(dim(D1)[2]-2)/2
	x=2*(1:H)+1; y=x+1; 
	A=array(dim=c(H, H, 2))	
	Lambda=A
		
	r1=(D1$samp==id[1]) & (D1$key==id[2])
	r2=(D2$samp==id[1]) & (D2$key==id[2])
	S.cols=substring(names(D2), 2, 2)%in%as.character(0:9)
	E=data.matrix(D1[r1, x])
	Eo=data.matrix(D1[r1, y]) 
	
	A[,,1]=est.AF(Eo, E)
	W=get.weights(Eo, E)
	gmm<-function(Ahat){
		N=t(t(Eo)-Ahat%*%t(E))
		M=c(t(E)%*%N)
		M%*%W%*%M/dim(E)[1]
	}

	o.TN=est.TN(gmm, H)		
	if(o.TN$convergence==0) A[,,2]=o.TN$A	
	
	Lambda[,,1]=est.Lambda(Eo, E, A[,,1])
	Lambda[,,2]=est.Lambda(Eo, E, A[,,2])
	
	sic=D2$sic[r2]
	biz=1+(sic<5200)*(sic>=5000)+2*(sic<4000)
	
	theta=array(dim=2)
	S=matrix(as.real(D2[r2,S.cols]), sqrt(sum(S.cols)));
	a.b=o.TN$a.b
	Omega=array(c(D2$vari[r2], D2$covii[r2], D2$covii[r2], D2$vari[r2]), dim=c(2,2))
	Xi=array(c(D2$varp[r2], D2$covoo[r2], D2$covoo[r2], D2$varp[r2]), dim=c(2,2))
	Psi=array(c(D2$vard[r2], D2$covds[r2], D2$covds[r2], D2$vars[r2]), dim=c(2,2))
	C=makeC(H)
	theta[1]=sqrt(Xi[1,1])/sqrt(Omega[1,1])*a.b[1]
	theta[2]=sqrt(sum(diag(C%*%(A[,,2]%*%S%*%t(A[,,2])+Lambda[,,2])%*%t(C))[1]))/sqrt(Omega[1,1])*a.b[2]
	
	print(id)
	list(id=id, A=A, S=S, Lambda=Lambda, W=W, a.b=a.b, theta=theta, Omega=Omega, Xi=Xi, Psi=Psi, m=D2$m[r2], sic=sic, biz=biz, r=substr(id[2], nchar(id[2]), nchar(id[2])))
}

#All-feasible estimator
est.AF<-function(Eo, E) t(Eo)%*%E%*%solve(t(E)%*%E)

#tripple-newsvendor estimator
est.TN<-function(gmm, H){
	objective.TN<-function(log.a.b){
		norm.penalty=(sum(exp(log.a.b))/10^4)^2		
		gmm(getA(H, exp(log.a.b), c(0,1)))+norm.penalty
	}	
	count=0
	while(count<5){
		o=optim(log(.8^(1:2))+rnorm(2)/3, objective.TN, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$a.b=exp(o$par)
	o$A=getA(H, o$a.b, c(0,1))
	o
}

est.Lambda<-function(Eo, E, A){
	N=t(t(Eo)-A%*%t(E))
	t(N)%*%N/dim(E)[1]
}

get.weights<-function(Eo, E){
#solve with all-feasible estimator to get initial A, then find moment conditions under that A
	N=t(t(Eo)-est.AF(Eo, E)%*%t(E))
	S=array(0, dim=c(dim(E)[2]^2, dim(E)[2]^2))
	for(i in 1:dim(E)[1]){
		m=c(cbind(E[i,])%*%rbind(N[i,]))
		S=S+cbind(m)%*%rbind(m)
	}
	#Add 1% of the mean diagonal element to all the elements, to make S invertable
	S=S+.01*diag(dim(S)[1])*mean(diag(S))
	solve(S/dim(E)[1])
}

logit.A<-function(log.A){	
#creates a matrix in which all elements are positive and sum to less than one
	k=sqrt(length(log.A))
	m=array(exp(log.A), dim=c(k,k))
	t(t(m)/(apply(m, 2, sum)+1))
}