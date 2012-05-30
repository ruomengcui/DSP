sim.mle<-function(k, sim.count, gamma.scaler){	
#k[1]=length of signals; k[2]=number of theta parameters
	samps=unique(D1$samp)
	samps=samps[samps]
	biz=sics; biz$biz=1*(sics$sic>=5200)+2*(sics$sic<5200)*(sics$sic>=5000)+3*(sics$sic<4000)
	est.samp=load.est.samp(biz, gamma.scaler, k, sim.count)
	mclapply(samps, est.samp, mc.preschedule=FALSE)
}

load.est.samp<-function(biz, gamma.scaler, k, sim.count){
	est.samp<-function(my.samp){
		group=list()
		for(i in 1:3){
			g=biz[biz$biz==i,1]
			group[[i]]=g[g%in%unique(D1[D1$samp==my.samp, "key"])]
		}
		est.group=load.est.group(my.samp, gamma.scaler, group, k, sim.count)
		o=lapply(seq(group), est.group)
		save(o, file=paste(varSave, 'simMLE_', my.samp, d.tag, '.txt', sep=''))
	}
}

load.est.group<-function(my.samp, gamma.scaler, group, k, sim.count){
	est.group<-function(group.id){
		x=2*(1:k[1])+1; y=x+1; 
		my.group=group[[group.id]];	firm.count=length(my.group)
		demand=list(); order=list(); gs=array(dim=c(firm.count, k[2]))
		print(c(my.samp, my.group))		
		for(i in seq(my.group)){
			r1=(D1$samp==my.samp) & (D1$key==my.group[i])
			r2=(gamma.scaler[,1]==my.samp) & (gamma.scaler[,2]==my.group[i])
			demand[[i]]=D1[r1, x]
			order[[i]]=D1[r1, y]
			gs[i,]=gamma.scaler[r2, seq(k[2])+2]
		}
		draws=array(rnorm(sim.count*k[2]*firm.count), dim=c(k[2], sim.count, firm.count))
		mle.objective=load.mle.objective(draws, demand, order, gs, k, firm.count, sim.count)
		optim(c(log(rep(.3, k[2])), c(diag(k[2])[lower.tri(diag(k[2]), diag=TRUE)]), c(diag(k[1])[lower.tri(diag(k[1]), diag=TRUE)])), mle.objective, gr=NULL, method="L-BFGS-B", lower=c(rep(-11, k[2]+k[2]*(k[2]+1)/2), rep(-Inf, k[1]*(k[1]+1)/2)), upper=c(rep(11, k[2]+k[2]*(k[2]+1)/2), rep(Inf, k[1]*(k[1]+1)/2)), control=list(maxit=50000))
	} 
}

load.mle.objective<-function(draws, demand, order, gs, k, firm.count, sim.count){
	mle.objective<-function(param){
		print(param[1:6])
		p1=k[2]; p2=k[2]*(k[2]+1)/2; p3=k[1]*(k[1]+1)/2
		theta.LD=diag(k[2]); Lambda.LD=diag(k[1])
		theta.LD[lower.tri(diag(k[2]), diag=TRUE)]<-param[seq(p2)+p1]
		Lambda.LD[lower.tri(diag(k[1]), diag=TRUE)]<-param[seq(p3)+p1+p2]
		gamma=vapply(seq(firm.count), function(x) gs[x,]*exp(param[seq(k[2])]+theta.LD%*%draws[,,x]), draws[,,1])
		#limit size of gamma to 10^10+1, so getA doesn't get screwed up
		gamma[gamma>10^10]<-10^10+(gamma-10^10)/(gamma-10^10+1)
		As=apply(gamma, c(2,3), function(x) getA(c(x, rep(0, k[1]-k[2]+1))))
		get.firm.loglikelihood<-function(firm.num){
			E=t(demand[[firm.num]]); Eo=t(order[[firm.num]])
			sim.loglikelihood<-function(my.sim){
				A=matrix(As[, my.sim, firm.num], k[1])
				sum(dmnorm(t(Eo-A%*%E), varcov=Lambda.LD%*%t(Lambda.LD), log=T))
			}
			logs=vapply(seq(sim.count), sim.loglikelihood, 0)
			log(mean(exp(logs-max(logs))))+max(logs)
		}
		-sum(vapply(seq(firm.count), get.firm.loglikelihood, 0))
	}	
}	

load.mle.estimates<-function(k){
	my.files=dir(varSave)
	my.files=my.files[substr(my.files, 1, 6)=="simMLE"]
	for(i in seq(my.files)){
		load(paste(varSave, my.files[i], sep=""))
		if(!exists("param")) param=array(dim=c(3*length(my.files), 2+length(o[[1]]$par)))
		trial=as.integer(substr(my.files[i],8,9))
		param[3*i-2,]=c(trial, 1, o[[1]]$par)
		param[3*i-1,]=c(trial, 2, o[[2]]$par)
		param[3*i,]=c(trial, 3, o[[3]]$par)
	}
	get.dist<-function(stats){
		p1=k[2]; p2=k[2]*(k[2]+1)/2
		theta.LD=diag(k[2]);
		theta.LD[lower.tri(diag(k[2]), diag=TRUE)]<-stats[2+seq(p2)+p1]
		c(stats[seq(2+k[2])], sqrt(diag(theta.LD%*%t(theta.LD))))
	}
	o=as.data.frame(t(apply(param, 1, get.dist)))
	names(o)<-c("samp", "biz", "m0", "m1", "m2", "s0", "s1", "s2")
	o=melt(o, c("samp", "biz"))
 	o$theta=substr(o$variable,2,2)
 	o$variable=substr(o$variable,1,1)
 	o=cast(o, ...~variable)
	list(full=param, short=o)
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