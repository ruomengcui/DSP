gmm.chi.squared<-function(my.pols){
	attach(my.pols)
	attach(D.eps)
	gv=unique(gvkey)
	x=2*(1:dim(A)[2]); y=x+1;
	diff.stat=array(dim=dim(A)[3])
	for (i in 1:dim(A)[3]){
		Y=data.matrix(D.eps[D.eps[,1]==gv[i], y]) 
		X=data.matrix(D.eps[D.eps[,1]==gv[i], x])
		diff.stat[i]=calc.diff.stat(Y, X, A[,,i,3], W[,,i])
	}
	detach(my.pols)
	detach(D.eps)
	diff.stat
}

calc.diff.stat<-function(Eo, E, A, W){
	N=t(t(Eo)-A%*%t(E))
	m=array(0, dim=c(dim(E)[2]^2))
	for(i in 1:dim(E)[1]){
		m=m+c(cbind(E[i,])%*%rbind(N[i,]))
	}
	m=m/dim(E)[1]
	m%*%W%*%m
}