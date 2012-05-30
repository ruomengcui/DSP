calc.policies<-function(){
	attach(D.eps)
	gv=unique(gvkey)
	x=2*(1:4); y=x+1;
	A=array(dim=c(length(y), length(x), length(gv), 3))
	gamma=array(dim=c(length(gv), length(x)+1))
	for (i in 1:length(gv)){
		Y=data.matrix(D.eps[D.eps[,1]==gv[i], y]) 
		X=data.matrix(D.eps[D.eps[,1]==gv[i], x])		
		W=get.weights(Y, X)
		A[,,i,1]=est.A1(Y, X)
		#o=est.A2(Y, X, W)
		#if(o$convergence==0) A[,,i,2]=o$par
		o2=est.TN(Y, X, W)		
		if(o2$convergence==0){
			A[,,i,3]=o2$A
			gamma[i,]=o2$gamma
		}	
		print(c(i, o$convergence, o2$convergence))
	}
	detach(D.eps)
	list(A=A, gamma=gamma)
}

calc.Lambda<-function(A){
	attach(D.eps)
	gv=unique(gvkey)
	x=2*(1:4); y=x+1;
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
	S=array(dim=c(4, 4, length(gvkey)))
	Omega=array(dim=c(2, 2, length(gvkey)))
	Xi=array(dim=c(2, 2, length(gvkey)))
	Psi=array(dim=c(2, 2, length(gvkey)))
	for (i in 1:length(gvkey)){
		S[,,i]=array(data.matrix(D.stats[i,4:39]), dim=c(6,6))[1:4, 1:4]
		Omega[,,i]=c(vari[i], covii[i], covii[i], vari[i])
		Xi[,,i]=c(varp[i], covoo[i], covoo[i], varp[i])
		Psi[,,i]=c(vard[i], covds[i], covds[i], vars[i])
	}
	o=list(gvkey=gvkey, sic=sic, S=S, Omega=Omega, Xi=Xi, Psi=Psi, m=m)	
	detach(D.stats)	
	o
}

#All-feasible estimator
est.A1<-function(Eo, E) t(Eo)%*%E%*%solve(t(E)%*%E)

#economically-sound estimator
est.A2<-function(Eo, E, W){
	k=dim(E)[2]
	c1=diag(k^2)
	for(i in 1:k){
		t=array(0, dim=c(k,k))
		t[,i]=-1
		c1=rbind(c1, c(t))
	}
	c2=c(rep(0,k^2), rep(-1,k))
	gmm=def.gmm(Eo, E, W)
	count=0
	while(count<4){
		starting=(.5 + .5*runif(1))*(c(.95*diag(k))+rep(exp(-5), k^2))
		o=constrOptim(starting, gmm, NULL, c1, c2, outer.iterations=3000, outer.eps=1e-05, control=list(maxit=3000))
		count=count+1+100*(o$convergence==0)
	}
	o
}

#tripple-newsvendor estimator
est.TN<-function(Eo, E, W){
	k=dim(E)[2]
	gmm=def.gmm(Eo, E, W)
	objective.TN<-function(log.gamma){
		norm.penalty=(sum(exp(log.gamma)^2)/10^4)^2
		gmm(c(getA(exp(log.gamma), 0)))+norm.penalty
	}	
	count=0
	while(count<5){
		starting=log(c(.3, .2, .1, .05, .05))+rnorm(5)/3
		o=optim(starting, objective.TN, NULL, method="Nelder-Mead", control=list(maxit=3000))
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
	N=t(t(Eo)-est.A1(Eo, E)%*%t(E))
	S=array(0, dim=c(dim(E)[2]^2, dim(E)[2]^2))
	for(i in 1:dim(E)[1]){
		m=c(cbind(E[i,])%*%rbind(N[i,]))
		S=S+cbind(m)%*%rbind(m)
	}
	solve(S/dim(E)[1])
}

def.gmm<-function(Eo, E, W){
	gmm<-function(Ahat){
		k=sqrt(length(Ahat))
		N=t(t(Eo)-array(c(Ahat), dim=c(k,k))%*%t(E))
		M=c(t(E)%*%N)
		M%*%W%*%M
	}
	gmm
}
	
#Dubuggers
debug.code<-function(){
	source('getA'); 
	A=getA(c(.5,.5,.25,.1,.1,0,0,0), 0);
	E=array(rnorm(700), dim=c(1000, 7))
	N=array(rnorm(700), dim=c(1000, 7))
	Eo=t(A%*%t(E)+t(N))
	get.weights(Eo, E)
}