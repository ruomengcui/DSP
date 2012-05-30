options(error=utils::recover)
if(file.exists("/Volumes/d1")){
	s.d1='/Volumes/d1/'
	s.d3='/Volumes/d3/getA/'
	} else{
	s.d1='/afs/ir.stanford.edu/users/r/l/rlbray/drives/d1/'
	s.d3='/data/gsb/rlbray/getA/'
	.libPaths(paste(s.d1, 'code/R/R_lib', sep=""))
}
dataSource<-paste(s.d3, 'inputs/', sep="") 
varSave<-paste(s.d3, 'variables/', sep="")
finalOuts<-paste(s.d3, 'outputs/', sep="")
global.sed.folder<-paste(s.d1, 'code/R/sed_files/', sep="")
local.sed.folder<-paste(s.d1, 'code/getA/R/sed_files/', sep="")
ex.funs<-paste(s.d1, 'code/R/functions/', sep="")
ex.mods<-paste(s.d1, 'code/getA/R/modules/', sep="")
library('multicore')
library('reshape')
library('ggplot2')
library('mnormt')
library('expm')
library('MASS')
library('numDeriv')
library('BB')
library('bayesm')

my.funs=dir(ex.funs)
for(i in seq(my.funs)) source(paste(ex.funs, my.funs[i], sep=""))
my.mods=dir(ex.mods)
for(i in seq(my.mods)) source(paste(ex.mods, my.mods[i], sep=""))
