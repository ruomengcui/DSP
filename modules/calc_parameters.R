calc.parameters<-function(data.in, data.out, which.alpha, which.beta, L=0, mc=T){
	load(paste(varSave, data.in, sep=""))
	num.firms=dim(unique(d[,c("samp", "cid")]))[1]
	if(mc) o=mclapply(seq(num.firms), calc.firm, mc.preschedule=FALSE, d, L, which.alpha, which.beta)
	else o=llply(seq(num.firms), calc.firm, d, L, which.alpha, which.beta)
	save(o, file=paste(varSave, data.out, sep=''))
}

calc.firm<-function(id, d, L, which.alpha, which.beta){
	my.firm=unique(d[,c("samp", "cid")])[id,]
	print(my.firm)
	f=d[d$samp==my.firm$samp & d$cid==my.firm$cid,]
	E=as.matrix(f[,which(substr(names(f), 1, 2)=="e_")])
	Eo=as.matrix(f[,which(substr(names(f), 1, 2)=="o_")])
	Z=as.matrix(f[,which(substr(names(f), 1, 2)=="z_")])
	r.ests=c(list(firm=my.firm), reduced.form.ests(E, Eo, Z))
	gmm.objective<-function(log.alpha.beta){
		param=unlog.param(log.alpha.beta, length(which.alpha))
		A=solve.for.A(list(alpha=param$alpha, which.alpha=which.alpha, beta=param$beta, which.beta=which.beta, Sigma=r.ests$Sigma, Lambda=r.ests$Lambda, Gamma=r.ests$Gamma, L=L))$A
		norm.penalty=sum((exp(log.alpha.beta)/8)^10)
		sum(cov(get.eta(A, E, Eo), Z)^2)+norm.penalty
	}
  s.ests=unlog.param(structural.ests(gmm.objective, length(which.alpha)+length(which.beta)),length(which.alpha))
	c(r.ests, s.ests)
}

reduced.form.ests<-function(E, Eo, Z){	
	A=t(Eo)%*%Z%*%t(Z)%*%E%*%solve(t(E)%*%Z%*%t(Z)%*%E)
	N=get.eta(A, E, Eo)
	list(A=A, Sigma=cov(E), Lambda=cov(N), Gamma=cov(E, N))
}

structural.ests<-function(gmm.objective, num.params){
	num.iterations=5; count=0
	while(count<num.iterations){
		o=optim(log(runif(num.params)), gmm.objective, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+num.iterations*(o$convergence==0)
	}
	o$par
}

get.eta<-function(A, E, Eo) t(t(Eo)-A%*%t(E))

unlog.param<-function(log.alpha.beta, num.alpha){
	alpha=exp(log.alpha.beta[seq(num.alpha)])
	beta=exp(log.alpha.beta[(num.alpha+1):length(log.alpha.beta)])
	list(alpha=alpha, beta=beta)
}