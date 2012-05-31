plot.Impulse.Responses<-function(A, strip.text.labels, file.name, interquartile=1, tick=.5){
	stats=array(dim=c(4, dim(A)[1:3]))
	stats[1,,,]=apply(A, 1:3, mean)
	stats[2:4,,,]=apply(A, 1:3, quantile, probs=c(.5, .25, .75))
	d=data.frame(melt(stats, varnames=c("stat", "m", "n", "est")))
	d$stat<-as.factor(d$stat)
	d$est<-as.factor(d$est)
	d$n<-as.factor(d$n)
	levels(d$stat)<-c('Mean', 'Median', 'Q1', 'Q3')
	if(!interquartile) d=d[d$stat%in%c('Mean', 'Median')	,]
	levels(d$est)<-strip.text.labels
	levels(d$n)<-c('n==0', 'n==1', 'n==2', 'n==3', 'n==4', 'n==5', 'n==6', 'n==7')
	d$m=d$m-1
	d$stat<-factor(d$stat)
	p<-ggplot(d, aes(x=m, y=value, colour=stat)) +geom_line()+ geom_hline(aes(yintercept=0), size=.1) +facet_grid(n~est, labeller = label_parsed)+opts(legend.title = theme_text(size = 0))+ scale_y_continuous(expression(paste("Matrix's ", mn^"th", " element")), breaks=seq(-100,100,tick))+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90), legend.text=theme_text(size = 13), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.7,start=0)
	ggsave(paste(finalOuts, "impulse_responses/", file.name, ".png", sep=""),  dpi=300)
}

plot.Inventory.IRFs<-function(A){
	k=dim(A)[1]
	avg.IRF=apply(A, c(1,2,3), mean)[,,1]
	inv.IRFs<-function(x){ 
		IRFs=array(dim=c(k,k,2))
		IRFs[,,1]=makeC(k)%*%(A[,,1,x]-diag(k))
		IRFs[,,2]=makeC(k)%*%(avg.IRF-diag(k))%*%A[,,1,x]
		IRFs
	}		
	M=vapply(seq(dim(A)[4]), inv.IRFs, array(0,dim=c(k,k,2)))
	plot.Impulse.Responses(M,  c(expression(paste("C(", D^L, widehat(A), "-I)")), expression(paste("C(", D^L[s], bar(A), "-I)", widehat(A)))), "inv", interquartile=0, tick=.25)
}

plot.Cov.Matricies<-function(d, fourier=0){
	stats=array(dim=c(2, 2, dim(d[[1]]$S)))
	if(fourier){
		F=get.fourier.matrix(dim(d[[1]]$S)[1])
		Gamma=vapply(d, function(s) F%*%s$S%*%Conj(F), F)
		Xi=vapply(d, function(s) F%*%s$Lambda[,,2]%*%Conj(F), F)
		stats[1,1,,]=apply(Re(Gamma), c(1,2), quantile, probs=c(.5))
		stats[1,2,,]=apply(Im(Gamma), c(1,2), quantile, probs=c(.5))	
		stats[2,1,,]=apply(Re(Xi), c(1,2), quantile, probs=c(.5))
		stats[2,2,,]=apply(Im(Xi), c(1,2), quantile, probs=c(.5))
		my.labs=c(expression(widehat(Gamma)), expression(widehat(Xi)))
		my.labs2=c("Real","Imaginary")
		my.out="covmat_fourier"
	}else{
		S=vapply(d, function(s) s$S, d[[1]]$S)
		Lambda=vapply(d, function(s) s$Lambda[,,2], d[[1]]$Lambda[,,2])
		stats[1,1,,]=apply(S, c(1,2), mean)
		stats[1,2,,]=apply(S, c(1,2), quantile, probs=c(.5))	
		stats[2,1,,]=apply(Lambda, c(1,2), mean)
		stats[2,2,,]=apply(Lambda, c(1,2), quantile, probs=c(.5))
		my.labs=c(expression(widehat(Sigma)), expression(widehat(Lambda)))
		my.labs2=c('Mean', 'Median')
		my.out="covmats"
	}
	d=data.frame(melt(stats, varnames=c("est", "stat", "m", "n")))
	d$stat<-as.factor(d$stat)
	d$est<-as.factor(d$est)
	d$n<-as.factor(d$n)
	levels(d$est)<-my.labs
	levels(d$stat)<-my.labs2
	levels(d$n)<-c('n==0', 'n==1', 'n==2', 'n==3', 'n==4', 'n==5', 'n==6', 'n==7')
	d$m=d$m-1
	p<-ggplot(d, aes(x=m, y=value, colour=stat)) +geom_line() +facet_grid(n~est, labeller = label_parsed)+opts(legend.title = theme_text(size = 0))+ scale_y_continuous(expression(paste("Matrix's ", mn^"th", " element")), breaks=seq(0,50, 10))+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90), legend.text=theme_text(size = 13), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.7,start=0)
	ggsave(paste(finalOuts, "impulse_responses/", my.out, ".png", sep=""),  dpi=300)
}

plot.freq.IRF<-function(H=8){
	S=diag(H); Lambda=array(0, dim=c(H,H)); C=makeC(H); W=get.fourier.matrix(H)
	cases=t(matrix(c(0,0,  0.001,1,  0.001,2,  1,0,  2,0,  1,1,  2,2), nr=2))
	d=vapply(seq(dim(cases)[1]), function(x) W%*%solve.for.A(cases[x,], S, Lambda)$A%*%Conj(W), array(complex(r=0, i=0), dim=c(H,H)))
	cases=floor(cases)
	r=melt(Re(d), varnames=c("m", "n", "a.b"))
	i=melt(Im(d), varnames=c("m", "n", "a.b"))
	names(r)[which(names(r)=="value")]="Real"
	names(i)[which(names(i)=="value")]="Imaginary"
	d=merge(r, i)
	d=melt(d, id.vars=c("m", "n", "a.b"))
	d$m=d$m-1; d$n=d$n-1
	d=d[d$n<=H/2,]
	d$n=as.factor(d$n)
	levels(d$n)<-c('n==0', 'n==1', 'n==2', 'n==3', 'n==4', 'n==5', 'n==6', 'n==7')
	d$a.b=as.factor(d$a.b)
	levels(d$a.b)=vapply(seq(dim(cases)[1]), function(x) paste(expression(alpha),"==", cases[x,2], "~", expression(beta),"==", cases[x,1]), "")
	d$variable=factor(d$variable, levels=c("Real", "Imaginary"))
	p<-ggplot(d, aes(x=m, y=value, colour=variable))+geom_line() +facet_grid(n~a.b, labeller = label_parsed)+opts(legend.title = theme_text(size = 0),legend.text = theme_text(size = 9),legend.position=c(.87,.56))+ scale_y_continuous(expression(paste("Matrix's ", mn^"th", " element")), breaks=c(0,1))+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 10), strip.text.y=theme_text(size = 11, angle=-90), strip.text.x=theme_text(size = 10), axis.title.x=theme_text(size = 13), axis.title.y=theme_text(size = 13, angle=90), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.7,start=0)
	ggsave(paste(finalOuts, "impulse_responses/freqIRF.png", sep=""),  dpi=300)
}

plot.freq.vars<-function(H=20){
	get.freq.vars<-function(H, a.b){
		S=diag(H); Lambda=array(0, dim=c(H,H)); C=makeC(H); J=array(0, dim=c(H,H)); J[1,1]=1; W=get.fourier.matrix(H)
		A=solve.for.A(c(a.b[2], a.b[1]), S, Lambda)$A
		M=W%*%A%*%Conj(W)
		freq.l.vars<-function(l){
			gamma=array(0, dim=c(H,H)); gamma[l,l]=1
			gamma.o=M%*%gamma%*%Conj(t(M))
			gamma.i=W%*%C%*%Conj(W)%*%(M-diag(H))%*%gamma%*%t(Conj(M)-diag(H))%*%W%*%t(C)%*%Conj(W)
			gamma.j=W%*%J%*%Conj(W)%*%M%*%gamma%*%t(Conj(M))%*%W%*%J%*%Conj(W)
			c(a.b, l-1, Re(sum(diag(gamma.o))), Re(sum(diag(gamma.i))), Re(sum(diag(gamma.j))))
		}	
		t(vapply(seq(H), freq.l.vars, rep(0,6)))
	}
	cases=expand.grid(c(.1, .5, 1), c(.1, 1, 2))
	g=do.call(rbind, lapply(seq(dim(cases)[1]), function(x) get.freq.vars(H, as.numeric(cases[x,]))))
	colnames(g)<-c("a", "b", "frequency", "vo", "vi", "vm")
	g=as.data.frame(g)
	g$a=paste(expression(alpha), "==", g$a); g$b=paste(expression(beta), "==", g$b)
	g$vm=(H-.1)*g$vm
	g=g[g$frequency<=H/2,]
	g=melt(g, id.vars=c("a", "b", "frequency"))
	my.labs<-list(bquote(paste(o[t],":                                  ",MJ[l], M^"*'")), bquote(paste(i[t],": ", WCW^"*","(",M-I,")",J[l],"(",M^"*'",-I,")",WC^"'",W^"*")),bquote(paste(i[t]^m,":          ",HWJ[1],W^"*",MJ[l],M^"*'",WJ[1],W^"*")))
	levels(g$variable)=0:2
	p<-ggplot(g, aes(frequency, value, colour=variable))+geom_line() +facet_grid(b~a, labeller = label_parsed)+geom_point()+labs(x="Frequency l", y="")+opts(legend.position=c(.8,.8), legend.title = theme_text(size = 0),legend.text = theme_text(size = 8),panel.background=theme_rect(colour="lightgrey", size=1), panel.background = theme_rect(colour=NA))+scale_colour_manual(values=c("0"="black","1"="darkgrey","2"="lightgrey"), breaks=0:2, labels=my.labs)
	ggsave(paste(finalOuts, "freq/freqvars.png", sep=""),  dpi=300)
}

plot.Inventory.deviations<-function(A){
	k=dim(A)[1]
	inv.IRFs<-function(x){ 
		IRFs=makeC(k)%*%(A[,,1,x]-diag(k))
	}		
	invIRF=vapply(seq(dim(A)[4]), inv.IRFs, array(0,dim=c(k,k)))
	avg.IRF=apply(invIRF, c(1,2), mean)
	S=vapply(f.0, function(s) s$S, f.0[[1]]$S)
	avg.S=apply(S, c(1,2), mean)
	inv.diviations=abs(avg.IRF%*%mat.sqrt(avg.S)%*%matrix(rnorm(80000), 8))
	stats=array(dim=c(4, 8))
	stats[1,]=apply(inv.diviations, 1, mean)
	stats[2:4,]=apply(inv.diviations, 1, quantile, probs=c(.5, .25, .75))
	d=data.frame(melt(stats, varnames=c("stat", "x")))
	d$stat=as.factor(d$stat)
	levels(d$stat)<-c('Mean', 'Median', 'Q1', 'Q3')
	d$x=d$x-1
	p<-ggplot(d, aes(x=x, y=value, colour=stat)) +geom_line()+opts(legend.title = theme_text(size = 0))+scale_y_continuous(expression(paste("|", e[l],"'",epsilon[t]^i,"|",)))+scale_x_continuous("Lead time l")+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90), legend.text=theme_text(size = 13), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.7,start=0)
	ggsave(paste(finalOuts, "impulse_responses/deviations.png", sep=""),  dpi=300)
}

plot.BwBreakdown<-function(){
	S=vapply(f.0, function(s) s$S, f.0[[1]]$S)
	Lambda=vapply(f.0, function(s) s$Lambda[,,2], f.0[[1]]$Lambda[,,2])
	d=t(vapply(seq(length(f.0)), function(x){ dsp=sum(diag(A.0[,,1,x]%*%S[,,x]%*%t(A.0[,,1,x])-S[,,x])); noise=sum(diag(Lambda[,,x])); c(dsp, noise)}, c(0,0)))
	d=as.data.frame(pmax(pmin(d, 100), -50))
	d$noisy[d$V2>d$V1]<-1; d$noisy[d$V2<=d$V1]<-0
	d$noisy=as.factor(d$noisy)
	my.labs <- list(bquote(paste(widehat(B)["SA"],">",widehat(B)["NC"])),bquote(paste(widehat(B)["SA"]<=widehat(B)["NC"])))
	ggplot(d, aes(y=V2, x=V1, colour=noisy))+geom_point()+scale_colour_manual(values=c("0"="black","1"="lightgrey"),breaks=0:1, labels=my.labs)+opts(legend.position=c(.85,.9), legend.title = theme_text(size =0), legend.text = theme_text(size = 14), axis.title.x=theme_text(size = 20),axis.title.y=theme_text(size = 20, angle=90),axis.text.x=theme_text(size = 16),axis.text.y=theme_text(size = 16),panel.background=theme_rect(colour=NA))+coord_cartesian(ylim=c(-2,102), xlim=c(-52,52))+labs(y=expression(widehat(B)["NC"]), x=expression(widehat(B)["SA"]))+scale_x_continuous(breaks=seq(-50,50,25))+scale_y_continuous(breaks=seq(0,100,25))	
	ggsave(paste(finalOuts, "RcaseStudies/BwBreakdown.png", sep=""),  dpi=300)
}

plot.metrics<-function(m){	
	m=melt(m, id.var = c('id.key', 'sic'))
	m=m[m$variable%in%c("adh1", "prod"),]	
	m$variable<-factor(m$variable, levels=c("adh1", "prod"))
	m$sic=as.factor(floor(m$sic/100))
	m=wins(m, "value", c("variable"), q=c(0, .92))
	m=wins(m, "value", c("sic", "variable"), q=c(.1, .9))
	m=normalize(m, INDICIES=m$variable)
	m[m$variable=="adh1","value"]=1-m[m$variable=="adh1","value"]
	print(multicol.unique(m, c(3,5)))
	levels(m$variable)<-c('Plan Adherence', 'Production  Stability')
	f<-function(x) {
		r<-quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
		names(r)<-c("ymin", "lower", "middle", "upper", "ymax")
		r
	}
	m$temp=""
	m=merge(m, get.sic.labels())
	ggplot(m, aes(x=sic.label, y=value)) + stat_summary(fun.data = f, geom="boxplot")+ facet_grid(temp~variable, scales="free")+coord_flip()+opts(legend.position="none", axis.text.x=theme_text(colour="black"))+labs(x="", y="")+scale_y_continuous(breaks=c(-2,2))
	ggsave(paste(finalOuts, "RcaseStudies/box.png", sep=""),  dpi=300)
}

plot.gap.closures<-function(d){	
	d$delta=1-d$inv1
	d=d[,c("biz", "delta")]
	x=sort(unique(d$delta))
	CDFs=by(d, d$biz, function(s) ecdf(s$delta)(x))
	d=as.data.frame(do.call(rbind, lapply(seq(CDFs), function(i) cbind(x, CDFs[[i]], i))))
	names(d)<-c("x", "y", "biz")
	d$biz=factor(d$biz); levels(d$biz)<-c("Retail", "Wholesale", "Manufacturing")
	p<-ggplot(d, aes(colour=biz, y=y, x=x))+geom_line(size=1)+coord_cartesian(xlim=c(0,1.05))+scale_colour_grey(end=0,start=0.7)+opts(panel.background = theme_rect(colour=NA), axis.text.y=theme_text(size=18), axis.text.x=theme_text(size = 18), axis.title.x=theme_text(size = 24), axis.title.y=theme_text(size = 24, angle=90))+scale_y_continuous(breaks=seq(0,1,.25))+scale_x_continuous(breaks=seq(0,1,.25))+opts(legend.position=c(.85,.75), legend.text=theme_text(size = 18))+labs(colour="", y="CDF Across Firms", x="Fraction of Gap Closed Within One Quarter") 
	ggsave(paste(finalOuts, "invGap/gap.png", sep=""),  dpi=300)	
}

plot.frontier<-function(d){	
	d=d[,c("sic", "adh1", "prod")]
	d=d[as.logical((d$adh1<80)*(d$prod<80)),]
	d$adh1=-d$adh1
	d$biz=c("Retail", "Wholesale", "Manufacturing")[1+(d$sic<5200)*(d$sic>=5000)+2*(d$sic<4000)]
	d$biz=factor(d$biz, levels=c("Retail", "Wholesale", "Manufacturing"))
	p<-ggplot(d, aes(x=adh1, y=prod))+stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE)+geom_density2d(bins = 12, size = .5, colour = "black")+scale_fill_gradient(limit = c(0.0003, .0025), low = "white", high="darkgrey")+scale_x_continuous(breaks=seq(-120, 120, 20))+scale_y_continuous(breaks=seq(0, 95, 5))+opts(legend.position=c(.25,.87), legend.title=theme_text(size = 0), axis.text.y=theme_text(size=15), axis.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 22), axis.title.y=theme_text(size = 22, angle=90), legend.text=theme_text(size = 16), panel.background = theme_rect(colour=NA))+labs(x=expression(paste(U^o, '(', widehat(A), ",",widehat(Sigma), ",",widehat(Lambda), ")")), y=expression(paste('S(', widehat(A), ",",widehat(Sigma), ",",widehat(Lambda), ")")), biz="")+coord_cartesian(ylim=c(0,30), xlim=c(0,-83))
	ggsave(paste(finalOuts, "invGap/frontier.png", sep=""),  dpi=300)	
}

plot.Preference.Parameters<-function(){	
	d=alpha.beta.0; names(d)<-c("samp", "key", "biz", "sic", "alpha", "beta")
	my.jitter=.013; d$alpha=d$alpha+my.jitter; d$beta=d$beta+my.jitter; 
	d$alpha=pmin(d$alpha, 1.5-my.jitter); d$beta=pmin(d$beta, 1.5-my.jitter)
	p<-ggplot(d, aes(x=alpha, y=beta))+geom_point(size=1.25, alpha=1/3, position=position_jitter(w=my.jitter, h=my.jitter))+ geom_hline(aes(yintercept=quantile(beta, probs=c(.25,.5,.75))),size=.3)+ geom_vline(aes(xintercept=quantile(alpha, probs=c(.25,.5,.75))), size=.3)  +opts(legend.position=c(.25,.87), legend.title=theme_text(size = 0), axis.text.y=theme_text(size=15), axis.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 24), axis.title.y=theme_text(size = 24, angle=90), legend.text=theme_text(size = 16), panel.background = theme_rect(colour=NA))+labs(x=expression(widehat(alpha)), y=expression(widehat(beta)), biz="")+scale_fill_gradient(low = "white", high="darkgrey")+scale_x_continuous(breaks=c(0,.5,1,1.5))+scale_y_continuous(breaks=c(0,.5,1,1.5))
	ggsave(paste(finalOuts, "preferenceParameters/pref.png", sep=""),  dpi=300)	
}

plot.cf.bw<-function(d){
	d$bw=d$o21-d$o11; d$bw0=d$o22-d$o12; d=d[,c("id.key", "bw", "bw0")]; names(d)[1]="key"
	d=merge(d, alpha.beta.0[,c("key", "V6")])
	names(d)[4]="beta"
	d=d[d$beta<2,]
	d=melt(d, id=c("key", "beta"))
	d$beta=-d$beta
	levels(d$variable)<-c(expression(paste("B(", A[1], ",", Sigma[1],",", Lambda[1],",", Gamma[1], ")-B(", A[0],",", Sigma[0],",", Lambda[0],",", Gamma[0], ")")), expression(paste(B^0, "(", A[1],",", Sigma[1],",", Lambda[1],",", Gamma[1], ")-", B^0,"(", A[0], ",",Sigma[0],",", Lambda[0],",", Gamma[0], ")")))
	d$temp=""
	p<-ggplot(d, aes(x=beta, y=value))+facet_grid(temp~variable, labeller = label_parsed)+geom_point(size=1, alpha=1/2)+geom_smooth(colour="black")+ geom_hline(aes(yintercept=0),size=.3)+opts(panel.background=theme_rect(colour="lightgrey"), legend.title=theme_text(size = 0), axis.text.y=theme_text(size=13), axis.text.x=theme_text(size = 13), axis.title.x=theme_text(size = 18),strip.text.x=theme_text(size=15), axis.title.y=theme_text(size = 20, angle=90))+labs(x=expression(paste(beta[1]-beta[0])), y="")+scale_y_continuous(breaks=seq(-30, 40, 5))+coord_cartesian(ylim=c(-7,27), xlim=c(.05, -1.7))
	ggsave(paste(finalOuts, "cf/cf_bw.png", sep=""),  dpi=300)	
}

plot.cf.forc<-function(d){
	d=pmax(d, -22)
	p<-ggplot(d, aes(x=V1, y=V2))+geom_point(size=1) + stat_quantile(colour="black") +opts(panel.background=theme_rect(colour="lightgrey"), legend.title=theme_text(size = 0), axis.text.y=theme_text(size=13), axis.text.x=theme_text(size = 13), axis.title.x=theme_text(size = 20),strip.text.x=theme_text(size=17), axis.title.y=theme_text(size = 20, angle=90))+opts(panel.background=theme_rect(colour="lightgrey", size=1))+scale_x_continuous(breaks=seq(-20, 0, 5))+labs(x=expression(paste(U, "(",Sigma[1], ")-", U,"(",Sigma[0], ")")), y=expression(paste(U^o, "(", A[1],",", Sigma[1],",", Lambda[1],",", Gamma[1], ")-", U^o,"(", A[0], ",",Sigma[0],",", Lambda[0],",", Gamma[0], ")")))
	ggsave(paste(finalOuts, "cf/cf_forc.png", sep=""),  dpi=300)	
}

plot.metric.medians<-function(t, bound, file, h.0=0, tick.length=.25){
	t$sic=as.factor(floor(t$sic/100))
	t=cast(t, ...~., median)
	names(t)[which(names(t)=="(all)")]="median"
	s=cast(t, sic+stat~., value="median", fun.aggregate=sd)
	names(s)[which(names(s)=="(all)")]="se"
	t=t[t$samp==0,]
	t=merge(t, s)
	t=merge(t, get.sic.labels())
	t$sic.label<-factor(t$sic.label)
	t$u=pmin(pmax(t$median+1.64*t$se, bound[1]), bound[2])
	t$l=pmin(pmax(t$median-1.64*t$se, bound[1]), bound[2])
	t$median=pmin(pmax(t$median, bound[1]), bound[2])
	t$sig[t$l>h.0]<-"Significant"
	t$sig[t$l<=h.0]<-"Insignificant"
	t$sig=factor(t$sig, levels=c("Significant", "Insignificant"))
	t$temp=""
	p<-ggplot(t, aes(colour=sig, y=median, x=sic.label))
	limits<-aes(ymax=u, ymin=l)
	p + geom_point() + geom_errorbar(limits, width=0.7)+ geom_hline(aes(yintercept=0),size=.3)+ facet_grid(temp~stat, labeller = label_parsed, scales="free_y")+scale_y_continuous(breaks=seq(from=bound[1], to = bound[2], by = tick.length))+labs(x="", y="")+opts(panel.background=theme_rect(colour="lightgrey", size=1), axis.line=theme_segment(colour="black"), legend.title = theme_text(colour = "white"), legend.position = "none", axis.text.x=theme_text(colour="black"), strip.text.x=theme_text(size=13), panel.background = theme_rect(colour=NA))+scale_colour_grey(end=0.7,start=0)+coord_flip(ylim=c(0,bound))
	ggsave(paste(finalOuts, file, sep=""),  dpi=300)	
}

plot.hausman<-function(d){	
	d=as.data.frame(d)
	d$R2=apply(A.0, 4, pull.R2)
	d$p=pchisq(22, d[,2])
	d$p=pmax(10^(-5), d$p);
	ggplot(d, aes(x=R2, y=p))+geom_point(size=I(.8))+scale_y_log10(breaks=10^(0:-4)) +opts(axis.text.x=theme_text(colour="black", size=11), axis.text.y=theme_text(colour="black", size=11), axis.title.x=theme_text(size=15), axis.title.y=theme_text(size=15, angle=90), panel.background = theme_rect(colour=NA))+scale_x_continuous(breaks=seq(from=0, to = 1, by = .2))+labs(x=expression(R^2), y="p-value")
	ggsave(paste(finalOuts, "hausman/haus.png", sep=""),  dpi=300)
}

plot.F<-function(){
	N=6
	i=complex(real=0, imaginary=1)
	delta=.01
	w=exp(-2*pi*i/(N/delta))
	sines=N^(-1/2)*sapply(seq(N)-1, function(x) (w^(seq((N/delta))-1))^x)
	d=melt(sines, varnames=c("m", "n"))
	d$m=(d$m-1)*delta
	d$int=d$m==floor(d$m)
	r=as.data.frame(d); r$value=as.real(Re(r$value)); r$type="Real"; 
	i=as.data.frame(d); i$value=as.real(Im(i$value)); i$type="Imaginary"; 
	d=rbind(r,i); d$m=as.real(d$m); d$n=as.real(d$n)
	d$type=factor(d$type, levels=c("Real", "Imaginary"))
	d$n<-as.factor(d$n)
	levels(d$n)<-c('n==0', 'n==1', 'n==2', 'n==3', 'n==4', 'n==5', 'n==6', 'n==7')
	d$temp=""
	p<-ggplot(d, aes(x=m, y=value, colour=type, size=int)) +geom_point() +facet_grid(n~temp, labeller = label_parsed)+opts(legend.title = theme_text(size = 0))+ scale_y_continuous(expression(paste(mn^"th", " element of W")), breaks=seq(-1,1,.2))+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90), legend.text=theme_text(size = 13), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.8,start=0)
	ggsave(paste(finalOuts, "impulse_responses/fourier.png", sep=""),  dpi=300)
}

get.sic.labels<-function(){
	levels=c('Hardware and garden (R)',
		'General merchandise (R)',
		'Food (R)',
		'Automotive dealers (R)',
		'Apparel and accessory (R)',
		'Homefurnishings (R)',
		'Restaurants (R)',
		'Miscellaneous (R)',
		'Durable goods (W)',
		'Nondurable goods (W)'	,
		'Tobacco (M)',
		'Textile mill (M)',
		'Apparel (M)',
		'Lumber and wood (M)',
		'Furniture and fixtures (M)',
		'Paper (M)',
		'Printing and publishing (M)',
		'Chemicals (M)',
		'Petroleum and coal (M)',
		'Rubber and plastics (M)',
		'Leather goods (M)',
		'Stone and glass (M)',
		'Primary metal (M)',
		'Fabricated metal (M)',
		'Industrial machinery (M)',
		'Electronic  equipment (M)',
		'Transportation equipment (M)',
		'Instruments and related (M)',
		'Miscellaneous (M)')
	o=data.frame(sic=c(52:59, 50:51, 21:39), sic.label=levels)
	levels(o$sic.label)=levels[length(levels):1]
	o
}

####Old###
plot.counterfactual<-function(m, append=""){	
	names(m)[1]<-"sic"
	m=melt(m, id.var='sic')
	m=m[m$variable%in%c("adh", "stock", "lead0", "inv", "prod", "sinv"),]
	m$variable<-factor(m$variable, levels=c("adh", "stock", "lead0", "inv", "prod", "sinv"))
	m$sic=as.factor(floor(m$sic/100))
	m=wins(m, "value", c("sic", "variable"), q=c(.1, .9))
	m=wins(m, "value", c("variable"), q=c(.05, .95))
	m=normalize(m, INDICIES=m$variable)
	m$value=m$value+(m$tick<0)	
	print(multicol.unique(m, c(2,4)))
	levels(m$variable)<-c('Adherence', 'Stock Var', 'Lead-0 BW', 'Inv. Exceed', 'Order Exceed', "Supplier Inv.")
	f<-function(x) {
		r<-quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
		names(r)<-c("ymin", "lower", "middle", "upper", "ymax")
		r
	}
	m$temp=""
	m=merge(m, get.sic.labels())
	ggplot(m, aes(x=sic.label, y=value)) + stat_summary(fun.data = f, geom="boxplot")+ facet_grid(temp~variable, scales="free") +coord_flip()+opts(legend.position="none", axis.text.x=theme_text(colour="black"))+labs(x="", y="")+scale_y_continuous(breaks=c(0,1))+opts(strip.text.x=theme_text(size=8.5))
	ggsave(paste(finalOuts, "RcaseStudies/counterfactual", append, ".png", sep=""),  dpi=300)
}

plot.theta<-function(mle, f.0){
	x=seq(.0005, .95, .001)
	u=1; 
	d1=mle[mle$samp==0,-which(names(mle)=="samp")]
	d1=do.call(rbind, lapply(seq(dim(d1)[1]), function(i) cbind(d1[i, 1:2], x, plnorm(x, d1[i, "m"], d1[i, "s"]))))
	names(d1)<-c("sic", "theta", "x", "y")
	d2=do.call(rbind, lapply(f.0, function(x) c(x$sic, x$theta[1:3])))
	findCDF<-function(col) do.call(rbind, lapply(1:3, function(i) cbind(i, col-2, x, ecdf(d2[d2[,1]==i, col])(x))))
	d2=as.data.frame(do.call(rbind, lapply(2:4, findCDF)))
	names(d2)<-names(d1)
	d1$est<-"MLE"; d2$est<-"GMM"
	d=rbind(d1, d2)
	d$theta=factor(d$theta)
	levels(d$theta)<-c(expression(theta), expression(theta[1]), expression(theta[2]))
	d$sic=factor(d$sic)
	levels(d$sic)<-c("Retail", "Wholesale", "Manufacturing")
	p<-ggplot(d, aes(colour=est, y=y, x=x))
	limits<-aes(ymax=u)
	p + geom_line()+facet_grid(theta~sic, labeller = label_parsed, scales="free_y")+coord_cartesian(xlim=c(0,.5), ylim=c(0,1))+labs(x="", y="")+scale_colour_grey(end=0,start=0.7)+opts(panel.background = theme_rect(colour=NA), axis.text.y=theme_text(size=10, face="bold"), axis.text.x=theme_text(size = 11, face="bold"), strip.text.y=theme_text(size = 13), strip.text.x=theme_text(size = 13), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90))+scale_y_continuous(breaks=seq(.2,1,.2))+scale_x_continuous(breaks=seq(.1,.4,.1))+labs(colour="Estimator", y=expression(CDF), x="Newsvendor Ratio Value") 
	ggsave(paste(finalOuts, "theta/theta_MLE_plot.png", sep=""),  dpi=300)
}

plot.newspapers<-function(){
	my.jitter=.04
	ind=c(3674, 3823:3841, 2711:2721, 5411)
	d=alpha.beta.0[alpha.beta.0$V4%in%ind, c("V4", "V5", "V6")]; 
	d$V5=d$V5+my.jitter; d$V6=d$V6+my.jitter;
	d$Industry=factor(c("Newspapers/Magazine Publishing", "Grocery Retailing", "Semiconductor Manufacturing", "Instruments Manufacturing")[1+(d$V4%in%5411)+2*(d$V4%in%3674)+3*(d$V4%in%3823:3841)], levels=c("Semiconductor Manufacturing", "Instruments Manufacturing", "Newspapers/Magazine Publishing", "Grocery Retailing"))
	d$V5=pmin(1.05, d$V5); d$V6=pmin(.8, d$V6);
	p<-ggplot(d, aes(y=V6, x=V5))
	p+geom_point(size=1.25, position=position_jitter(w=my.jitter, h=my.jitter))+facet_wrap(~Industry, nrow=2)+scale_colour_grey()+labs(x=expression(widehat(alpha)),y=expression(widehat(beta)))+scale_x_continuous(breaks=seq(0,1,.2))+scale_y_continuous(breaks=seq(0,.8,.2))+opts(panel.background=theme_rect(colour="lightgrey"), panel.background = theme_rect(colour=NA), axis.text.y=theme_text(size=10, face="bold"), axis.text.x=theme_text(size = 11, face="bold"), strip.text.y=theme_text(size = 13), strip.text.x=theme_text(size = 13), axis.title.x=theme_text(size = 18), axis.title.y=theme_text(size = 18))
	ggsave(paste(finalOuts, "theta/newspapers.png", sep=""),  dpi=300)
}

plot.theory.policies<-function(H, margins, strip.labels, range = 5^c(-2:2)){
#margins is two numbers telling which two variables are "turned on". The first four numbers relate to a's and the rest to b's.
	S=diag(H+1); Lambda=matrix(0, H+1, H+1); Gamma=Lambda; L=0
	x=expand.grid(range, range)
	cases=matrix(0, dim(x)[1], max(c(margins, 5)))
	cases[,margins[1]]=x[,1]
	cases[,margins[2]]=x[,2]
	calc.A<-function(x){
		A=solve.for.A(list(alpha=x[1:4], beta=x[5:length(x)], Sigma=S, Lambda=Lambda, Gamma=Gamma, L=L))$A
		cbind(rbind(x), melt(A))
	}
	d=adply(cases, 1, calc.A)
	d=d[,c(margins, (dim(d)[2]-2):dim(d)[2])]
	names(d)[which(names(d)%in%c("X1", "X2"))]<-c("m", "n")
	d$n=factor(d$n); 
	names(d)[1:2]=c("m1", "m2")		
	d$m1<-as.factor(d$m1); d$m2<-as.factor(d$m2); 		
	levels(d$m1)<-strip.labels[[1]]
	levels(d$m2)<-strip.labels[[2]]
	p<-ggplot(d, aes(x=m, y=value, colour=n)) +geom_line(size=.25)+facet_grid(m1~m2, labeller = label_parsed, scales="free")+ labs(y=expression(paste(e[m], "'", A, e[n])), x="m")+opts(axis.text.y=theme_text(size=11), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 18), axis.title.y=theme_text(size = 18, angle=90), legend.text=theme_text(size = 15),legend.title=theme_text(size = 18), panel.background=theme_rect(colour="lightgrey"))+scale_colour_grey(end=0.9,start=0)+scale_x_continuous(breaks=seq(0,15,2))
	ggsave(paste(finalOuts, "theoryPics/theory", paste(margins, collapse = ''), ".png", sep=""),  dpi=300)
}

plot.mc<-function(trials, graph.name, which.alpha, which.beta){
  d=ldply(trials, function(x){
    load(paste(varSave, 'simPanel', x, '.txt', sep=""))
    load(paste(varSave, 'simEst', x, '.txt', sep=""))
    p=unique(d[,names(d)=="cid"| substr(names(d), 1, 5)=="alpha"|substr(names(d), 1, 4)=="beta"])
    e=ldply(o, function(x){
      if(attributes(x)[1]=="try-error") return() #exit if Sigma isn't singular
      o=cbind(x$firm, rbind(x$alpha), rbind(x$beta))
      colnames(o)<-c("samp", "cid", paste("alpha_", which.alpha, sep=""), paste("beta_", which.beta, sep=""))
      o
    })
    s=ddply(e[,-which(names(e)=="samp")], "cid", function(x) sd(x)[2:dim(x)[2]])
    e=e[e$samp==0,-which(names(e)=="samp")]
    names(p)[2:dim(p)[2]]<-paste(names(p)[2:dim(p)[2]], "_p", sep="")
    names(e)[2:dim(e)[2]]<-paste(names(e)[2:dim(e)[2]], "_e", sep="")
    names(s)[2:dim(s)[2]]<-paste(names(s)[2:dim(s)[2]], "_s", sep="")
    a=merge(p, merge(e, s))
    a=melt(a, id="cid")
    a$var=as.factor(paste(substr(a$variable, 1,1), substr(a$variable, regexpr("_", a$variable)+1, regexpr("_", a$variable)+1), sep=""))
    a$var<-factor(a$var, levels=c(paste("a", which.alpha, sep=""), paste("b", which.beta, sep="")))
    levels(a$var)<-c(laply(which.alpha, function(x) paste("alpha[", x, "]", sep="")), laply(which.beta, function(x) paste("beta[", x, "]", sep="")))
    a$t=substr(a$variable, nchar(as.character(a$variable)), nchar(as.character(a$variable)))
    a=a[,-which(names(a)==c("variable"))]
    a=cast(a, cid+var~t)
    a$u=a$e+1.96*a$s
    a$l=a$e-1.96*a$s
    a$c=a$p<a$u & a$p>a$l
    a$num.obs=x
    a
  })
  d$N=paste("T==", d$num.obs, sep="")
  d$N<-factor(d$N, levels=paste("T==", unique(d$num.obs)[order(unique(d$num.obs))], sep=""))
  ggplot(d,aes(y=e,x=p,colour=c))+geom_point()+geom_errorbar(aes(max=u,ymin=l))+geom_abline(intercept=0,slope=1)+facet_grid(var~N,labeller=label_parsed)+scale_colour_grey(end=0.7,start=0)+labs(x="Parameter", y="Estimate")+coord_cartesian(ylim=c(-.025,1.32), xlim=c(-.05, 1.05))+scale_x_continuous(breaks=seq(0, 1, .25))+scale_y_continuous(breaks=seq(0, 1.25, .25))+opts(panel.background=theme_rect(colour="lightgrey"), panel.background = theme_rect(colour=NA), legend.position = "none", strip.text.y=theme_text(size=12), strip.text.x=theme_text(size=12), axis.title.x=theme_text(size = 14), axis.title.y=theme_text(size = 14, angle=90), axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 10))
  ggsave(paste(finalOuts, "mc/", graph.name, ".png", sep=""),  dpi=300)	
  print(ddply(d, c("var", "N"), function(x) c(length(x$e), sum(!x$c), sum(x$e>x$p)/length(x$e), sum(x$e-x$p)/length(x$e))))
  print(c(length(d$e), sum(!d$c), sum(d$e>d$p)/length(d$e), sum(d$e-d$p)/length(d$e)))
}