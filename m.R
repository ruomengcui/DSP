source('/afs/ir.stanford.edu/users/r/l/rlbray/drives/d1/code/getA/R/header.R')
mc=1

#Monte Carlo
if(mc) run.monte.carlo()


################################
make.panels=0
load.structures=0
estimate.parameters=0
counterfactual=0
sim=0
spectral=0
test.h=0
make.plots=0
make.tables=0
debugger=0
load.data=0
###Prep Data###
if(make.panels){
	create.wards.panel()
}

###Calc Parameters###
if(estimate.parameters){
	calc.parameters('sim_signal_panel.txt')
	
	firms=mclapply(key, calc.firm, mc.preschedule=FALSE)
	firms=firms[sapply(firms, function(x) length(x)>1)]
	save(firms, file=paste(varSave, 'my_firms', d.tag, '.txt', sep=''))
	mtx=mclapply(1:length(firms), calc.metrics, firms, mc.preschedule=FALSE)
	save(mtx, file=paste(varSave, 'mtx', d.tag, '.txt', sep=''))
}

###load structures###
if(load.structures){
	load(paste(varSave, 'my_firms', d.tag, '.txt', sep=''))
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
	alpha.beta.0=as.data.frame(t(vapply(f.0, function(s){ unlist(c(s$id, s$biz, s$sic, s$alpha.beta))}, c(0,0,0,0,f.0[[1]]$alpha.beta))))
	alpha.beta.all=as.data.frame(t(vapply(firms, function(s){ unlist(c(s$id, s$biz, s$sic, s$alpha.beta))}, c(0,0,0,0,f.0[[1]]$alpha.beta))))
	names(alpha.beta.all)=c("samp", "key", "biz", "sic", expression(alpha), expression(beta))
	F=get.fourier.matrix(dim(firms[[1]]$S)[1])
	#gaps=inv.gap()
	#o=load.mle.estimates(c(4,2))
	#mle=o$short; mle.full=o$full
}

###Counterfactual###
if(counterfactual){
	cf=as.data.frame(do.call(rbind, lapply(lapply(f.0, counter.factual, metric.fun=m.bw, pref.fun=function(x){c(x[1],0)}), unlist)))
	cf2=lapply(f.0, counter.factual, metric.fun=m.ASA, S.fun=improve.forecast)
	cf2.aviv=as.data.frame(t(vapply(seq(cf2), function(s){ 
		k=length(cf2[[1]]$o1)
		D=makeD(k, 1)
		S=f.0[[s]]$S
		c(diag(D%*%S%*%t(D)-S)%*%(2^(-1:-k)), (cf2[[s]]$o2-cf2[[s]]$o1)%*%(2^(-1:-k)))
	}, c(0,0))))
}

###Make plots###
if(make.plots){
	#plot.Impulse.Responses(A.0, c(expression(widehat(A)), expression(paste(A^"*", "(", widehat(alpha), ",", widehat(beta),",", widehat(Sigma),",", widehat(Lambda),",", widehat(Gamma), ")"))), "impulse_responses")
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
	#plot.newspapers()	
	#plot.F()
	#plot.Preference.Parameters()
	#plot.cf.bw(cf)
	#plot.cf.forc(cf2.aviv)
	
	#median plots
	#to.plot=melt(m.all[, c("id.samp", "sic", "s1")], c("id.samp", "sic"))
	#names(to.plot)=c("samp", "sic", "stat", "value")
	#plot.metric.medians(to.plot, bound=c(-1, 3), file="u_and_s/s_med.png", tick.length=.5)
	#to.plot=melt(m.all[, c("id.samp", "sic", "u1")], c("id.samp", "sic"))
	#names(to.plot)=c("samp", "sic", "stat", "value")
	#plot.metric.medians(to.plot, bound=c(-.5, 2), file="u_and_s/u_med.png", tick.length=.5)
	plot.theory.policies(12, c(1, 2), list(c(expression(alpha[0]==5^-2), expression(alpha[0]==5^-1), expression(alpha[0]==5^0), expression(alpha[0]==5^1), expression(alpha[0]==5^2)), c(expression(alpha[1]==5^-2), expression(alpha[1]==5^-1), expression(alpha[1]==5^0), expression(alpha[1]==5^1), expression(alpha[1]==5^2))))
	plot.theory.policies(12, c(3, 4), list(c(expression(alpha[2]==5^-2), expression(alpha[2]==5^-1), expression(alpha[2]==5^0), expression(alpha[2]==5^1), expression(alpha[2]==5^2)), c(expression(alpha[3]==5^-2), expression(alpha[3]==5^-1), expression(alpha[3]==5^0), expression(alpha[3]==5^1), expression(alpha[3]==5^2))))
	plot.theory.policies(12, c(1, 8), list(c(expression(alpha[0]==5^-2), expression(alpha[0]==5^-1), expression(alpha[0]==5^0), expression(alpha[0]==5^1), expression(alpha[0]==5^2)), c(expression(beta[3]==5^-2), expression(beta[3]==5^-1), expression(beta[3]==5^0), expression(beta[3]==5^1), expression(beta[3]==5^2))))	
	plot.theory.policies(12, c(4, 8), list(c(expression(alpha[3]==5^-2), expression(alpha[3]==5^-1), expression(alpha[3]==5^0), expression(alpha[3]==5^1), expression(alpha[3]==5^2)), c(expression(beta[3]==5^-2), expression(beta[3]==5^-1), expression(beta[3]==5^0), expression(beta[3]==5^1), expression(beta[3]==5^2))))		
}

###Make tables###
if(make.tables){
	#table.pref.param()
	#table.year95()
	#my.metric.table("u1", "R_tables/u_metric.tex")
	#my.metric.table("s1", "R_tables/s_metric.tex")
	#u.metric.table()
	#table.BW.Break()
	table.cf()
}
