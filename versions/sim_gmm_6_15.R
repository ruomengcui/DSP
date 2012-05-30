sim.gmm<-function(k, sim.count, gamma.scaler){	
#k[1]=length of signals; k[2]=number of parameters to estimate
	samps=unique(D1$samp)
	sic2=sics; sic2$sic=floor(sics$sic/100); s=unique(sic2$sic)
	sim.samp=load.sim.samp(s, sic2, gamma.scaler, k, sim.count)
	lapply(samps, sim.samp)
}

load.sim.samp<-function(s, sic2, gamma.scaler, k, sim.count){
	sim.samp<-function(my.samp){
		group=list()
		for(i in seq(s)){
			g=sic2[sic2$sic==s[i],1]
			group[[i]]=g[g%in%unique(D1[D1$samp==my.samp, "gvkey"])]
		}
		sim.group=load.sim.group(my.samp, gamma.scaler, group, k, sim.count)
		o=mclapply(seq(group), sim.group, mc.preschedule=FALSE)
		save(o, file=paste(varSave, 'sim_samp_', my.samp, '.txt', sep=''))
		o
	}
}

load.sim.group<-function(my.samp, gamma.scaler, group, k, sim.count){
	sim.group<-function(group.id){
		x=2*(1:k[1])+1; y=x+1; 
		my.group=group[[group.id]];	firm.count=length(my.group)
		demand=list(); order=list(); gs=array(dim=c(firm.count, k[2]))
		print(c(my.samp, my.group))		
		for(i in seq(my.group)){
			r1=(D1$samp==my.samp) & (D1$gvkey==my.group[i])
			r2=(gamma.scaler[,1]==my.samp) & (gamma.scaler[,2]==my.group[i])
			demand[[i]]=D1[r1, x]
			order[[i]]=D1[r1, y]
			gs[i,]=gamma.scaler[r2, seq(k[2])+2]
		}
		draws=array(runif(sim.count*(k[2])*length(my.group)), dim=c(k[2], sim.count, firm.count))
		gmm=load.gmm(draws, demand, order, gs, k, firm.count, sim.count)
		o=optim(log(c(rep(1,k[2]), rep(.2, k[2]))), gmm, gr=NULL, W=diag(k[1]^2), out.W=FALSE, control=list(maxit=50000))
		W=gmm(o$par, diag(k[1]^2), out.W=TRUE)
		optim(o$par, gmm, gr=NULL, W=W, out.W=FALSE, control=list(maxit=50000))
	} 
}

load.gmm<-function(draws, demand, order, gs, k, firm.count, sim.count){
	gmm<-function(param, W, out.W){
		firm.gamma<-function(this.firm){
			g=array(dim=c(k[2], sim.count))
			for(i in seq(k[2])) g[i,]=qgamma(draws[i,,this.firm], exp(param[i]), exp(param[i+k[2]]))
			g[is.na(g)]<-10^-5
			gs[this.firm,]*g
		}
		gamma=vapply(seq(firm.count), firm.gamma, array(0, dim=c(k[2], sim.count)))
		get.firm.moments<-function(firm.num){
			E=t(demand[[firm.num]]); Eo=t(order[[firm.num]])
			mean.A=apply(vapply(seq(sim.count), function(x) getA(c(gamma[,x,firm.num], rep(0, k[1]-k[2]+1))), diag(k[1])), c(1,2), mean)
			N=Eo-mean.A%*%E		
			t(vapply(seq(dim(N)[2]), function(x) c(E[,x]%*%t(N[,x])), array(0, dim=k[1]^2)))
		}
		moments=do.call(rbind, lapply(seq(firm.count), get.firm.moments))
		norm.penalty=(mean(gamma)/1000)^3
		val=colMeans(moments)%*%W%*%colMeans(moments)+norm.penalty
		if(!out.W) return(val)
		solve(t(moments)%*%moments/dim(moments)[1])
	}	
}	

get.gamma.scale<-function(i, k){
	attach(firms[[i]])
	C=makeC(k[1])
	gs=array(dim=c(1, k[2]+2))
	gs[1:2]=id
	gs[3]=sqrt(Omega[1,1])/sqrt(Xi[1,1])
	for(j in seq(k[2]-1)){
		gs[3+j]=sqrt(Omega[1,1])/sqrt(sum(diag(C%*%(A[,,1]%*%S%*%t(A[,,1])+Lambda[,,1])%*%t(C))[1:j]))
	}
	detach(firms[[i]])
	gs
}