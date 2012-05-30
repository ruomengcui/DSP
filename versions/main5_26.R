dataSource <<- '/Volumes/d2/tables/getA/out/case_studies/'
dataDump <<- '../inputs/Routputs/'
library('reshape')
library('ggplot2')
library('mnormt') 
if(0){
	D.eps<-data.frame(read.table(paste(dataSource, 'toReg.txt', sep=""), header=TRUE))
	D.stats<-data.frame(read.table(paste(dataSource, 'sigmacov.txt', sep=""), header=TRUE))
	D.eps=D.eps[order(D.eps$gvkey), ]
	D.stats=D.stats[order(D.stats$gvkey), ]
	#D.eps=D.eps[D.eps$gvkey%in%unique(D.eps$gvkey)[1:5], ]
	#D.stats=D.stats[D.stats$gvkey%in%unique(D.eps$gvkey)[1:5], ]
}	
source('modules/getA.R')
source('modules/calc_parameters.R')
source('modules/calc_metrics.R')
source('modules/make_plots.R')
source('modules/h_tests.R')

#debug(calc.theta)

###Calc Parameters###
if(0){
	pols=calc.policies()
	save(pols, file='data/pols.Rdata')
	Lambda=calc.Lambda(pols$A)
	dat=calc.SigmaCov()
	dat$A=pols$A
	dat$gamma=pols$gamma
	dat$Lambda=Lambda
	dat$metrics=calc.metrics(dat)
	dat$theta=calc.theta(dat)
	dat$biz=get.biz(dat$sic)
	save(dat, file='data/dat.Rdata')
}


###Make plots###
if(1){
	#plot.metrics(dat)
	#plot.Impulse.Responses(pols$A)
	plot.theta(dat)
}

###Test Hypotheses###
if(0){
	test.outputs=list()
	pols$diff.stats=gmm.chi.squared(pols)
	test.outputs$diff.stats.05=sum(pols$diff.stats>30.14)/length(pols$diff.stats)
	test.outputs$diff.stats.10=sum(pols$diff.stats>27.2)/length(pols$diff.stats)
	test.outputs$diff.stats.total=c(pchisq(sum(pols$diff.stats), 19*length(pols$diff.stats)), sum(pols$diff.stats), qchisq(0.5, 19*length(pols$diff.stats)))
}	