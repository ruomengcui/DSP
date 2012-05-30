calc.parameters<-function(data.name, L=0, mc=T, param.count=c(1,1)){
	load(paste(varSave, data.name, sep=""))
	num.firms=dim(unique(d[,c("samp", "cid")]))[1]
	if(mc) o=mclapply(seq(num.firms), calc.firm, mc.preschedule=FALSE, d, L, param.count)
	else o=llply(seq(num.firms), calc.firm, d, L, param.count)
	o
}

calc.firm<-function(id, d, L, param.count){
	my.firm=unique(d[,c("samp", "cid")])[id,]
	print(my.firm)
	f=d[d$samp==my.firm$samp & d$cid==my.firm$cid,]
	E=as.matrix(f[,which(substr(names(f), 1, 2)=="e_")])
	Eo=as.matrix(f[,which(substr(names(f), 1, 2)=="o_")])
	Z=as.matrix(f[,which(substr(names(f), 1, 2)=="z_")])
	ests=c(list(firm=my.firm), reduced.form.ests(E, Eo, Z))
	W=get.weights(ests$N, Z)
	gmm.objective<-function(log.alpha.beta){
		norm.penalty=(sum(exp(log.alpha.beta))/10^3)^3
		param=unlog.alpha.beta(log.alpha.beta, param.count)
		A=solve.for.A(list(alpha=param$alpha, beta=param$beta, Sigma=ests$Sigma, Lambda=ests$Lambda, Gamma=ests$Gamma, L=L))$A
		moments=c(t(get.eta(A, E, Eo))%*%Z)
		moments%*%W%*%moments/dim(Z)[1]+norm.penalty
	}
	c(ests, structural.ests(gmm.objective, param.count))
}

reduced.form.ests<-function(E, Eo, Z){	
	A=t(Eo)%*%Z%*%t(Z)%*%E%*%solve(t(E)%*%Z%*%t(Z)%*%E)
	N=get.eta(A, E, Eo)
	list(A=A, Sigma=cov(E), Lambda=cov(N), Gamma=cov(E, N), N=N)
}

structural.ests<-function(gmm.objective, param.count){
	num.iterations=5; count=0
	while(count<num.iterations){
		o=optim(log(runif(sum(param.count))), gmm.objective, gr=NULL, method="Nelder-Mead", control=list(maxit=50000))
		count=count+1+num.iterations*(o$convergence==0)
	}
	c(unlog.alpha.beta(o$par, param.count), list(convergence=o$convergence))
}

get.weights<-function(N, Z){
#Get gmm weighting matrix
	moment.draws=laply(seq(dim(N)[1]), function(m){
		c(cbind(N[m,])%*%rbind(Z[m,]))
	})
	W=solve(cov(moment.draws))
	#Add 1% of the mean diagonal element to all the elements, to make W invertable
	W+.01*diag(dim(W)[1])*mean(diag(W))
}

get.eta<-function(A, E, Eo) t(t(Eo)-A%*%t(E))

unlog.alpha.beta<-function(log.alpha.beta, param.count){
	alpha=exp(log.alpha.beta[seq(param.count[1])])
	beta=exp(log.alpha.beta[param.count[1]+seq(param.count[2])])
	list(alpha=alpha, beta=beta)
}