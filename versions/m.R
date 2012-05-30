#d.tag=""#This one has a 5-period horizon.
#d.tag="_r" #This one is for before/after 1995
#d.tag="_short" #This one is for debugging.
d.tag="_h8" #This one is for spectral bw.

load.structures=0
solve.A=0
counterfactual=1
sim=0
spectral=0
test.h=0
make.plots=0
make.tables=0
debugger=0
load.data=0
source('/afs/ir.stanford.edu/users/r/l/rlbray/drives/d1/code/getA/R/header.R')

###Calc Parameters###
if(solve.A){
	firms=mclapply(key, calc.firm, mc.preschedule=FALSE)
	firms=firms[sapply(firms, function(x) length(x)>1)]
	save(firms, file=paste(varSave, 'firms', d.tag, '.txt', sep=''))
	mtx=mclapply(1:length(firms), calc.metrics, firms, mc.preschedule=FALSE)
	save(mtx, file=paste(varSave, 'mtx', d.tag, '.txt', sep=''))
}

###load structures###
if(load.structures){
	load(paste(varSave, 'firms', d.tag, '.txt', sep=''))
	load(paste(varSave, 'mtx', d.tag, '.txt', sep=''))
	load(paste(varSave, 'cf.txt', sep=''))
	load(paste(varSave, 'cf2.txt', sep=''))
	get.sics<-function(s) c(s$id[2], s$sic)
	sics=unique(as.data.frame(t(sapply(firms, get.sics))))
	names(sics)<-c("id.key", "sic")
	f.0=firms[sapply(firms, function(s){ s$id[1]==0})]
	m.0=merge(as.data.frame(t(sapply(mtx[sapply(mtx, function(s){ s$id[1]==0})], unlist))), sics)
	m.all=merge(as.data.frame(t(sapply(mtx, unlist))), sics)
	A.0=vapply(f.0, function(s){ s$A}, f.0[[1]]$A)
	theta.0=as.data.frame(t(vapply(f.0, function(s){ unlist(c(s$id, s$biz, s$sic, s$theta))}, c(0,0,0,0,f.0[[1]]$theta))))
	theta.all=as.data.frame(t(vapply(firms, function(s){ unlist(c(s$id, s$biz, s$sic, s$theta))}, c(0,0,0,0,f.0[[1]]$theta))))
	names(theta.all)=c("samp", "key", "biz", "sic", expression(alpha), expression(beta))
	F=get.fourier.matrix(dim(firms[[1]]$S)[1])
	#gaps=inv.gap()
	#o=load.mle.estimates(c(4,2))
	#mle=o$short; mle.full=o$full
}

###Counterfactual###
if(counterfactual){
	#cf=as.data.frame(do.call(rbind, lapply(lapply(f.0, counter.factual, metric.fun=m.bw, pref.fun=function(x){c(x[1],0)}), unlist)))
	#cf2=lapply(f.0, counter.factual, metric.fun=m.ASA, S.fun=improve.forecast)
	cf2.ASA=melt(vapply(seq(cf2), function(s){ 
		k=dim(cf2[[1]]$o1)[1]
		D=makeD(k, 1)
		S=f.0[[s]]$S
		x=array(0, dim=c(k,2))
		x[,1]=diag(D%*%S%*%t(D)-S)
		x[,2]=diag(cf2[[s]]$o2-cf2[[s]]$o1)
		x
	}, array(0, dim=c(dim(cf2[[1]]$o1)[1],2))))
}

###Simulated MLE###
if(sim){
	g.scaler=do.call(rbind, mclapply(seq(firms), get.gamma.scale, k=c(5,2), mc.preschedule=FALSE))
	sim.mle(k=c(4,2), sim.count=50, gamma.scaler=g.scaler)
}

###Spectral BW###
if(spectral){
	Gamma=lapply(seq(f.0), function(x) get.Gamma(f.0[[x]]))
	G=do.call(rbind, lapply(seq(Gamma), function(x) Gamma[[x]]$Gamma))
}

###Test Hypotheses###
if(test.h){
	test.outs=list()
	h=hausman.test()
	h=h[h[,2]>0,]
	test.outs$hausman=c(length(h[h<=qchisq(.99, 22)])/length(h), length(h[h<=qchisq(.95, 22)])/length(h))
	
	pull.R2<-function(s) var(c(s[,,3]))/(var(c(s[,,3]))+var(c(s[,,2])-c(s[,,3])))
	test.outs$R2=mean(apply(A.0, 4, pull.R2))
}

###Make plots###
if(make.plots){
	#plot.Impulse.Responses(A.0, c(expression(widehat(A)), expression(widehat(A)[alpha][beta])), "impulse_responses")
	plot.Impulse.Responses(cf2.ASA[,1:3,,], c(expression(paste(D, widehat(Sigma),D^"'","-", widehat(Sigma))),expression(paste(A[1],D, widehat(Sigma),D^"'", A[1]^"'","-",A[0], widehat(Sigma), A[0]^"'"))), "cf2", tick=10)
	#plot.Cov.Matricies(f.0)
	#plot.Cov.Matricies(f.0, fourier=1)
	#plot.Inventory.IRFs(A.0)
	#plot.freq.vars()
	#plot.freq.IRF()
	#plot.Inventory.deviations(A.0)
	#plot.BwBreakdown()
	#plot.metrics(m.0)
	#plot.gap.closures(gaps)	
	#plot.frontier(m.0)
	#plot.hausman(h)
	#plot.theta(mle, f.0)
	#plot.newspapers()	
	#plot.F()
	#plot.Preference.Parameters()
	#plot.cf.bw(cf)
	
	#median plots
	#to.plot=melt(m.all[, c("id.samp", "sic", "s1")], c("id.samp", "sic"))
	#names(to.plot)=c("samp", "sic", "stat", "value")
	#plot.metric.medians(to.plot, bound=c(-1, 3), file="u_and_s/s_med.png", tick.length=.5)
	#to.plot=melt(m.all[, c("id.samp", "sic", "u1")], c("id.samp", "sic"))
	#names(to.plot)=c("samp", "sic", "stat", "value")
	#plot.metric.medians(to.plot, bound=c(-.5, 2), file="u_and_s/u_med.png", tick.length=.5)
}

###Make tables###
if(make.tables){
	#theta.mle.table(mle)
	#table.pref.param()
	#table.year95()
	#my.metric.table("u1", "R_tables/u_metric.tex")
	#my.metric.table("s1", "R_tables/s_metric.tex")
	#u.metric.table()
	#table.BW.Break()
	table.cf()
}
