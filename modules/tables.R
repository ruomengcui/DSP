 my.metric.table<-function(met.name, out.name){
 	m=m.all
	m=m[,which(names(m)%in%c("id.samp", "sic", met.name))]
	m$sic=1+(m$sic<5200)*(m$sic>=5000)+2*(m$sic<4000)
	names(m)<-c(met.name, "samp", "biz")
	m2=m; m2$biz=0; m=rbind(m,m2)
	m=cast(m, samp+biz~., summary, value=met.name)
	m=as.data.frame(m[,which(names(m)%in%c("samp", "biz", "Median", "Mean"))])
	m=melt(m, id=c("samp", "biz"))
 	reg.table(m, variable~biz, boot=T, value="value", out.file=paste(finalOuts, out.name, sep=""), stars=c(.99, 1, 1), h.0=0, sed.files=c(paste(global.sed.folder, "stars.sed", sep=""), paste(local.sed.folder, "theta_MLE.sed", sep="")))
 }
 
 table.s_a<-function(){
 	d=m.all; d$biz=(1:3)[1+(d$sic<5200)*(d$sic>=5000)+2*(d$sic<4000)]
 	names(d)[which(names(d)=="id.samp")]<-"samp"
	d=melt(d[, c("samp", "biz", "s1")], c("samp", "biz"))
 	d$biz=factor(c("Retail", "Wholesale", "Manufacturing")[d$biz], levels=c("Retail", "Wholesale", "Manufacturing"))
	d=d[,-which(names(d)=="variable")]
 	meds=cast(d, samp+biz~., quantile, probs=.5, value="value")
 	means=cast(d, samp+biz~., mean, value="value")
 	meds$stat="Median"; means$stat="Mean"; d=rbind(meds, means)
 	reg.table(d, stat~biz, value="(all)", boot=T, out.file=paste(finalOuts, "R_tables/s_a.tex", sep=""), stars=c(.975, 1, 1), h.0=0, sed.files=c(paste(global.sed.folder, "stars.sed", sep="")))
}

table.BW.Break<-function(){
	d=as.data.frame(do.call(rbind, lapply(seq(length(firms)), function(x){
		A=firms[[x]]$A[,,1]
		S=firms[[x]]$S
		Lambda=firms[[x]]$Lambda[,,1]
		dsp=diag(A%*%S%*%t(A)-S)
		noise=diag(Lambda)
		o=cbind(rep(firms[[x]]$id[1], length(dsp)), seq(length(dsp))-1, dsp, noise, dsp+noise)
		o[4,]=colSums(o[4:5,])
		o[5,]=colSums(o[1:4,])
		o[1:5,]
	})))
	names(d)<-c("samp", "l", "dsp", "noise", "total")
	d=melt(d, id.vars=c("samp", "l"))
	d=cast(d, samp+l+variable~., mean)
	reg.table(d, l~variable, value="(all)", boot=T, out.file=paste(finalOuts, "R_tables/bw_break.tex", sep=""), round.digits=3, stars=c(.995, 1, 1), h.0=0, sed.files=c(paste(global.sed.folder, "stars.sed", sep="")))
}

table.pref.param<-function(){
	d=melt(alpha.beta.all[,c("samp", "biz", "alpha", "beta")], id.vars=c("samp", "biz"))
	d2=d; d2$biz=0; d=rbind(d,d2)
	names(d)[3]="param"
	d=melt(cast(d, samp+biz+param~., quantile, probs=c(.25, .5, .75)), id=c("samp", "biz", "param"))
	names(d)[5]="stat"
	d$stat=factor(d$stat, levels=c("X50.", "X25.", "X75."))
	d$biz=factor(c("Entire Economy", "Retail", "Wholesale", "Manufacturing")[d$biz+1], levels=c("Entire Economy", "Retail", "Wholesale", "Manufacturing"))
	reg.table(d, param+stat~biz, boot=T, value="value", out.file=paste(finalOuts, "R_tables/param_table.tex", sep=""), stars=c(.995, 1, 1), h.0=0,round.digits=2, sed.files=c(paste(global.sed.folder, "stars.sed", sep=""), paste(local.sed.folder, "theta_MLE.sed", sep="")))
}

table.cf<-function(){
	load(paste(varSave, 'cf2', d.tag, '.txt', sep=''))
	d=cf2
	d$i=d$V5-d$V3; d$s=d$V6-d$V4
	d=merge(d, m.0[,c("id.key", "sic")], by.x="key", by.y="id.key")
	d$biz=1+(d$sic<5200)*(d$sic>=5000)+2*(d$sic<4000)
	d2=d; d2$biz=0; d=rbind(d2,d)
	d=d[,c("samp", "biz", "i", "s")]
	d=melt(d, id=c("samp", "biz"))
	d=cast(d, samp+biz+variable~., summary)
	names(d)[3]<-"v"
	d=as.data.frame(d[,which(names(d)%in%c("v", "biz", "samp", "Mean", "Median"))])
	d=melt(d, id=c("v","samp", "biz"))
	d$variable<-factor(d$variable, levels=c("Mean", "Median"))
	reg.table(d, v+variable~biz, value="value", boot=T, out.file=paste(finalOuts, "R_tables/cf_forc.tex", sep=""), stars=c(.975, 1, 1), h.0=0, round.digits=3,sed.files=c(paste(global.sed.folder, "stars.sed", sep="")))
}

###Old###
 theta.mle.table<-function(m){
 	q=c(.05, .1, .25, .5, 1)
 	m=cbind(m, vapply(seq(q), function(i) plnorm(q[i], m$m, m$s), m$m))
 	names(m)[6:10]<-c("p05", "p10", "p25", "p50", "p1")
 	m=melt(as.data.frame(m[,-which(names(m)%in%c("m", "s"))]), id=c("samp", "biz", "theta"))
 	names(m)[which(names(m)=="variable")]<-"stat"
	m$stat<-factor(m$stat, levels=c("p05", "p10", "p25", "p50", "p1"))
 	m$biz=paste("b", m$biz, sep=""); m$theta=paste("t", m$theta, sep="")
 	reg.table(m, biz+theta~stat, boot=T, value="value", out.file=paste(finalOuts, "R_tables/theta_mle.tex", sep=""), stars=c(.995, 1, 1), h.0=1, sed.files=c(paste(global.sed.folder, "stars.sed", sep=""), paste(local.sed.folder, "theta_MLE.sed", sep="")))
 }
 
 table.year95<-function(){
	d=as.data.frame(t(sapply(firms, function(x) c(x$id[1], x$id[2], x$theta[1:3]))))
	d=merge(d, sics, by.x="key", by.y="id.key")
	d$biz=factor(c("Retail", "Wholesale", "Manufacturing")[1+(d$sic<5200)*(d$sic>=5000)+2*(d$sic<4000)], levels=c("Retail", "Wholesale", "Manufacturing"))
	d$r=10*zapsmall(d$key-floor(d$key)); d$key=floor(d$key)
	d=d[,-which(names(d)%in%c("sic", "key"))]
	d=melt(d, id=c("samp", "biz", "r"))
	d=cast(d, samp+biz+r+variable~., quantile, probs=.5, value="value")
	d=merge(d[d$r==0,], d[d$r==1,], by=c("samp", "biz", "variable"))
	d$val=d$"(all).y"-d$"(all).x"
	d=d[,-which(names(d)%in%c("r.x", "r.y", "(all).x", "(all).y"))]
	d.sd=cast(d, biz+variable~., sd, value="val")
	d=merge(d[d$samp==0,], d.sd, by=c("biz", "variable"))
	d=d[,-which(names(d)=="samp")]
	reg.table(d, biz~variable, value="val", se="(all)", out.file=paste(finalOuts, "R_tables/year95.tex", sep=""), stars=c(.975, 1, 1), h.0=0, sed.files=c(paste(global.sed.folder, "stars.sed", sep="")))
}