dataLoc<<-'/Volumes/d2/tables/getA/R/'
dataSource<<-paste(dataLoc, 'inputs/', sep="") 
dataDump<<-paste(dataLoc, 'outputs/', sep="")
varSave<<-paste(dataLoc, 'variables/', sep="")
library('reshape')
library('ggplot2')
library('mnormt') 
if(1){
	gv.list=unique(data.frame(read.table(paste(dataSource, 'toReg0.txt', sep=""), header=TRUE))$gvkey)
}	
source('modules/getA.R')
source('modules/calc_parameters.R')
source('modules/calc_metrics.R')
source('modules/make_plots.R')
source('modules/h_tests.R')

#debug(calc.SigmaCov)
print("yo")
###Calc Parameters###
if(0){
	o=calc.policies()
	#save(dat, file='data/o_booted.txt')
	#dat$Lambda=calc.Lambda(dat$A)
	#o=c(o, calc.SigmaCov())
	#o$metrics=calc.metrics(o)
	#o$theta=calc.theta(o)
	#o$biz=get.biz(o$sic)
	#save(o, file='data/o.Rdata')
}


###Make plots###
if(0){
	#plot.metrics(o)
	#plot.Impulse.Responses(pols$A)
	plot.theta(o)
}

###Test Hypotheses###
if(0){
	test.outputs=list()
	pols$diff.stats=gmm.chi.squared(pols)
	test.outputs$diff.stats.05=sum(pols$diff.stats>30.14)/length(pols$diff.stats)
	test.outputs$diff.stats.10=sum(pols$diff.stats>27.2)/length(pols$diff.stats)
	test.outputs$diff.stats.total=c(pchisq(sum(pols$diff.stats), 19*length(pols$diff.stats)), sum(pols$diff.stats), qchisq(0.5, 19*length(pols$diff.stats)))
}	