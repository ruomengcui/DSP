hausman.test<-function(f){
	f=as.data.frame(do.call(rbind, f))
	ids=as.data.frame(do.call(rbind, f$id))
	gvs=unique(ids$gvkey)
	haus<-function(my.firm){
		draws=f[f$id[2]==my.firm,]
		get.delta.theta<-function(my.draw){
			a=c(draws$A[my.draw][[1]][1]-draws$A[my.draw][[1]][2])
			b=c(draws$A[my.draw][[1]][1]-draws$A[my.draw][[1]][3])
			cbind(a,b)
		}
		m=c(cbind(E[i,])%*%rbind(N[i,]))
		S=S+cbind(m)%*%rbind(m)
	}
	mclapply(gvs, haus)	
}