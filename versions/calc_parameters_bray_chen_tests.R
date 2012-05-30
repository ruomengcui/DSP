calc.firm<-function(id){
	k=6
	x=2*(1:k)+1; y=x+1; 
	A=array(dim=c(k, k, 5))	
	Lambda=A
		
	r1=(D1$samp==id[1]) & (D1$gvkey==id[2])
	r2=(D2$samp==id[1]) & (D2$gvkey==id[2])
	E=data.matrix(D1[r1, x])
	Eo=data.matrix(D1[r1, y]) 
	
	A[,,1]=est.AF(Eo, E)
	W=get.weights(Eo, E)
	gmm<-function(Ahat){
		N=t(t(Eo)-Ahat%*%t(E))
		M=c(t(E)%*%N)
		M%*%W%*%M/dim(E)[1]
	}

	o.ES=est.ES(gmm, k)
	if(o.ES$convergence==0) A[,,2]=o.ES$A
	o.TN=est.TN(gmm, k)		
	if(o.TN$convergence==0) A[,,3]=o.TN$A	
	o.bray=est.bray(gmm, k)		
	if(o.bray$convergence==0) A[,,4]=o.bray$A	
	o.chen=est.chen(gmm, k)		
	if(o.chen$convergence==0) A[,,5]=o.chen$A		
	
	objectives=array()
	for(i in 1:dim(A)[3]){
		Lambda[,,i]=est.Lambda(Eo, E, A[,,i])
		objectives[i]=gmm(A[,,i])
	}	
	
	sic=D2$sic[r2]
	biz=1*(sic>=5200)+2*(sic<5200)*(sic>=5000)+3*(sic<4000)
	
	theta=array(dim=k+1)
	S=array(data.matrix(D2[r2,5:40]), dim=c(6,6))[1:k, 1:k]
	gamma=o.TN$gamma
	Omega=array(c(D2$vari[r2], D2$covii[r2], D2$covii[r2], D2$vari[r2]), dim=c(2,2))
	Xi=array(c(D2$varp[r2], D2$covoo[r2], D2$covoo[r2], D2$varp[r2]), dim=c(2,2))
	Psi=array(c(D2$vard[r2], D2$covds[r2], D2$covds[r2], D2$vars[r2]), dim=c(2,2))
	C=makeC(k)
	theta[1]=sqrt(Xi[1,1])/sqrt(Omega[1,1])*gamma[1]
	for(j in 1:length(gamma)-1){			
		theta[j+1]=sqrt(sum(diag(C%*%(A[,,3]%*%S%*%t(A[,,3])+Lambda[,,3])%*%t(C))[1:j]))/sqrt(Omega[1,1])*gamma[j+1]
	}
	print(id)
	list(id=id, A=A, S=S, Lambda=Lambda, W=W, gamma=gamma, theta=theta, Omega=Omega, Xi=Xi, Psi=Psi, m=D2$m[r2], sic=sic, biz=biz, objectives=objectives)
}

#All-feasible estimator
est.AF<-function(Eo, E) t(Eo)%*%E%*%solve(t(E)%*%E)

#economically-sound estimator
est.ES<-function(gmm, k){
	objective.ES<-function(log.A) gmm(logit.A(log.A))
	count=0
	while(count<5){
		o=optim(log(c(.8*diag(k))+.05)+rnorm(k^2)/2, objective.ES, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$A=logit.A(o$par)
	o
}

#tripple-newsvendor estimator
est.TN<-function(gmm, k){
	objective.TN<-function(log.gamma){
		norm.penalty=(sum(exp(log.gamma)^2)/10^4)^2
		gmm(getA(exp(log.gamma), 0))+norm.penalty
	}	
	count=0
	while(count<5){
		o=optim(log(.8^(1:(k+1)))+rnorm(k+1)/3, objective.TN, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$gamma=exp(o$par)
	o$A=getA(o$gamma, 0)
	o
}

#bray estimator
est.bray<-function(gmm, k){
	objective.bray<-function(log.gamma){
		norm.penalty=(sum(exp(log.gamma)^2)/10^4)^2
		gmm(getA(c(exp(log.gamma), rep(0, k)), 0))+norm.penalty
	}	
	count=0
	while(count<5){
		o=optim(log(.8+rnorm(1)/3), objective.bray, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$gamma=exp(o$par)
	o$A=getA(c(o$gamma, rep(0, k)), 0)
	o
}

#chen estimator
est.chen<-function(gmm, k){
	objective.chen<-function(log.gamma){
		norm.penalty=(sum(exp(log.gamma)^2)/10^4)^2
		gmm(getA(c(0, exp(log.gamma), rep(0, k-1)), 0))+norm.penalty
	}	
	count=0
	while(count<5){
		o=optim(log(.8+rnorm(1)/3), objective.chen, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$gamma=exp(o$par)
	o$A=getA(c(0, o$gamma, rep(0, k-1)), 0)
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