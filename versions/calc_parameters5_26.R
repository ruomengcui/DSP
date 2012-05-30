calc.policies<-function(){
	attach(D.eps)
	gv=unique(gvkey)
	x=2*(1:5); y=x+1;
	A=array(dim=c(length(x), length(x), length(gv), 3))
	W=array(dim=c(length(x)^2, length(x)^2, length(gv)))
	gamma=array(dim=c(length(gv), length(x)+1))
	for (i in 1:length(gv)){
		Y=data.matrix(D.eps[D.eps[,1]==gv[i], y]) 
		X=data.matrix(D.eps[D.eps[,1]==gv[i], x])		
		W[,,i]=get.weights(Y, X)
		A[,,i,1]=est.AF(Y, X)
		o.ES=est.ES(Y, X, W[,,i])
		if(o.ES$convergence==0) A[,,i,2]=o.ES$A
		o.TN=est.TN(Y, X, W[,,i])		
		if(o.TN$convergence==0){
			A[,,i,3]=o.TN$A
			gamma[i,]=o.TN$gamma
		}	
		print(c(i, o.ES$convergence, o.TN$convergence))
	}
	detach(D.eps)
	list(A=A, gamma=gamma, W=W)
}

calc.Lambda<-function(A){
	attach(D.eps)
	gv=unique(gvkey)
	x=2*(1:dim(A)[1]); y=x+1;
	Lambda=array(dim=c(length(y), length(x), length(gv), 3))
	for (i in 1:length(gv)){
		Y=data.matrix(D.eps[D.eps[,1]==gv[i], y]) 
		X=data.matrix(D.eps[D.eps[,1]==gv[i], x])		
		for (j in 1:3) Lambda[,,i,j]=est.Lambda(Y, X, A[,,i,j])
	}
	detach(D.eps)
	Lambda
}

calc.SigmaCov<-function(){
	attach(D.stats)
	S=array(dim=c(5, 5, length(gvkey)))
	Omega=array(dim=c(2, 2, length(gvkey)))
	Xi=array(dim=c(2, 2, length(gvkey)))
	Psi=array(dim=c(2, 2, length(gvkey)))
	for (i in 1:length(gvkey)){
		S[,,i]=array(data.matrix(D.stats[i,4:39]), dim=c(6,6))[1:5, 1:5]
		Omega[,,i]=c(vari[i], covii[i], covii[i], vari[i])
		Xi[,,i]=c(varp[i], covoo[i], covoo[i], varp[i])
		Psi[,,i]=c(vard[i], covds[i], covds[i], vars[i])
	}
	o=list(gvkey=gvkey, sic=sic, S=S, Omega=Omega, Xi=Xi, Psi=Psi, m=m)	
	detach(D.stats)	
	o
}

#All-feasible estimator
est.AF<-function(Eo, E) t(Eo)%*%E%*%solve(t(E)%*%E)

#economically-sound estimator
est.ES<-function(Eo, E, W){
	k=dim(E)[2]
	gmm=def.gmm(Eo, E, W)
	objective.ES<-function(log.A) gmm(logit.A(log.A))
	count=0
	while(count<5){
		o=optim(log(c(.8*diag(k))+.05)+rnorm(k^2)/2, objective.ES, NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$A=logit.A(o$par)
	o
}

#tripple-newsvendor estimator
est.TN<-function(Eo, E, W){
	k=dim(E)[2]
	gmm=def.gmm(Eo, E, W)
	objective.TN<-function(log.gamma){
		norm.penalty=(sum(exp(log.gamma)^2)/10^4)^2
		gmm(getA(exp(log.gamma), 0))+norm.penalty
	}	
	count=0
	while(count<5){
		o=optim(log(.8^(1:(k+1)))+rnorm(k+1)/3, objective.TN, NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+100*(o$convergence==0)
	}
	o$gamma=exp(o$par)
	o$A=getA(o$gamma, 0)
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
		S=S+cbind(m)%*%rbind(m)/dim(E)[1]
	}
	#Add 1% of the mean diagonal element to all the elements, to make S invertable
	S=S+.01*diag(dim(S)[1])*mean(diag(S))
	solve(S/dim(E)[1])
}

def.gmm<-function(Eo, E, W){
	gmm<-function(Ahat){
		N=t(t(Eo)-Ahat%*%t(E))
		M=c(t(E)%*%N)
		M%*%W%*%M
	}
	gmm
}

logit.A<-function(log.A){	
#creates a matrix in which all elements are positive and sum to less than one
	k=sqrt(length(log.A))
	m=array(exp(log.A), dim=c(k,k))
	t(t(m)/(apply(m, 2, sum)+1))
}

calc.theta<-function(o){
	attach(o)
	theta=array(dim=dim(gamma))
	C=makeC(dim(gamma)[2]-1)
	for(i in 1:dim(gamma)[1]){
		theta[i,1]=sqrt(Xi[1,1,i])/sqrt(Omega[1,1,i])*gamma[i,1]
		for(j in 1:dim(gamma)[2]-1){
			theta[i,j+1]=sqrt(sum(diag(C%*%(A[,,i,3]%*%S[,,i]%*%t(A[,,i,3])+Lambda[,,i,3])%*%t(C))[1:j]))/sqrt(Omega[1,1,i])*gamma[i,j+1]
		}
	}
	detach(o)
	theta
}

get.biz<-function(sic){
	b=array(dim=length(dat$sic))
	b[dat$sic<4000]<-'Manufacturing'
	b[(dat$sic<5200)&(dat$sic>=5000)]<-'Wholesale'
	b[dat$sic>=5200]<-'Retail'
	as.factor(b)
}