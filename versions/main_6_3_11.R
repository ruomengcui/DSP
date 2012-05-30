load.data=0
get.firms=0
make.plots=0
test.h=0
source('modules/header.R')

###Calc Parameters###
if(get.firms){
	f=mclapply(key, calc.firm, mc.preschedule=FALSE)
	save(f, file=paste(varSave, 'f', batch, '.txt', sep=''))
}
if(0){#combine batches
	firms=list();
	for(i in 1:4){
		load(file=paste(varSave, "f", i, ".txt", sep=""))
		firms=c(firms, f)	
		rm(f)
	}
	s<-function(sub) length(sub)>1
	firms=firms[sapply(firms, s)]
	mtx=mclapply(1:length(firms), calc.metrics, firms, mc.preschedule=FALSE)
	save(firms, file=paste(varSave, 'firms.txt', sep=''))
	save(mtx, file=paste(varSave, 'mtx.txt', sep=''))
}

###Make plots###
if(make.plots){
	#load(paste(varSave, 'firms.txt', sep=''))
	#load(paste(varSave, 'mtx.txt', sep=''))
	get.0<-function(s) s$id[1]==0
	get.sics<-function(s) c(s$id[2], s$sic)
	pull.A<-function(s) s$A
	pull.theta<-function(s) unlist(c(s$id, s$biz, s$theta))
	sics=unique(as.data.frame(t(sapply(firms, get.sics))))
	names(sics)<-c("id.gvkey", "sic")
	f.0=firms[sapply(firms, get.0)]
	m.0=merge(as.data.frame(t(sapply(mtx[sapply(mtx, get.0)], unlist))), sics)
	A.0=vapply(f.0, pull.A, f.0[[1]]$A)
	theta.0=as.data.frame(t(vapply(f.0, pull.theta, c(0,0,0,f.0[[1]]$theta))))
	#plot.Impulse.Responses(A.0)
	#plot.metrics(m.0)
	#plot.theta(theta.0)
	#save(f.0, file=paste(varSave, 'f_0.txt', sep=''))
}

###Test Hypotheses###
if(test.h){
	#load(paste(varSave, 'f_0.txt', sep=''))
	#D1=data.frame(read.table(paste(dataSource, 'toReg1.txt', sep=""), header=TRUE))
	#D1=D1[D1$samp==0,]
	test.outputs=list()
	debug(gmm.chi.squared)
	diff.stats=gmm.chi.squared(f.0)
	test.outputs$diff.stats.05=sum(diff.stats>30.14)/length(diff.stats)
	test.outputs$diff.stats.10=sum(diff.stats>27.2)/length(diff.stats)
	#test.outputs$diff.stats.total=c(pchisq(sum(pols$diff.stats), 19*length(pols$diff.stats)), sum(pols$diff.stats), qchisq(0.5, 19*length(pols$diff.stats)))
}	