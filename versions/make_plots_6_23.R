plot.Impulse.Responses<-function(A){
	A=A[,,2:3,]
	stats=array(dim=c(6, dim(A)[1:3]))
	stats[1,,,]=apply(A, c(1,2,3), mean)
	stats[2:6,,,]=apply(A, c(1,2,3), quantile, probs=c(.1, .25, .5, .75, .9))
	d=data.frame(melt(stats, varnames=c("stat", "m", "n", "est")))
	d$stat<-as.factor(d$stat)
	d$est<-as.factor(d$est)
	d$n<-as.factor(d$n)
	levels(d$stat)<-c('Mean', '10%', '25%', '50%', '75%', '90%')
	levels(d$est)<-c('Unrestricted', 'Restricted', 'Residual')
	levels(d$n)<-c('n=0', 'n=1', 'n=2', 'n=3', 'n=4', 'n=5')
	d$m=d$m-1
	p<-ggplot(d, aes(x=m, y=value, colour=stat)) +geom_line() +facet_grid(n~est)+opts(legend.title = theme_text(size = 0))+ scale_y_continuous(expression("A"["m,n"]))+opts(axis.text.y=theme_text(size=10), axis.text.x=theme_text(size = 11), strip.text.y=theme_text(size = 15, angle=-90), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90), legend.text=theme_text(size = 13), panel.background = theme_rect(colour=NA))+scale_colour_grey(end=0.9,start=0)
	ggsave(paste(finalOuts, "impulse_responses/impulse_responses.png", sep=""),  dpi=300)
}

plot.theta<-function(theta){
	t=theta[,c(3, 5:dim(theta)[2])]
	names(t)<-c("Sector", expression(theta), expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]))
	t$Sector<-factor(t$Sector)
	levels(t$Sector)<-c("Retail", "Wholesale", "Manufacturing")
	t=melt(t, id.var = "Sector")
	t$value=pmin(t$value, 3.75)
	p<-ggplot(t, aes(x=value))+geom_density(adjust=.34, fill="white")+facet_grid(variable~Sector, scales="free", labeller = label_parsed)+ xlab("Parameter Value")+ylab("Density")+coord_cartesian(xlim=c(0,.55))+scale_x_continuous(breaks=seq(from=0, to = 1, by = .2))+opts(legend.position="none", axis.text.y=theme_text(size=10, face="bold"), axis.text.x=theme_text(size = 11, face="bold"), strip.text.y=theme_text(size = 15), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90))
	ggsave(paste(finalOuts, "theta/theta_density.png", sep=""),  dpi=300)
}

plot.metrics<-function(m){	
	m$stock2=m$stock2/10
	m$inv=m$inv
	m=melt(m, id.var = c('id.gvkey', 'sic'))
	m=m[m$variable%in%c("adh2", "stock2", "lead02", "inv", "prod"),]	
	m$variable<-factor(m$variable, levels=c("adh2", "stock2", "lead02", "inv", "prod"))
	m$sic=as.factor(floor(m$sic/100))
	m=wins(m, "value", c("variable"), q=c(.1, .9))
	m=wins(m, "value", c("sic", "variable"), q=c(.1, .9))
	levels(m$variable)<-c('Adherence', 'Stock Var', 'Lead-0 BW', 'Inv. Exceed', 'Prod. Exceed')
	f<-function(x) {
		r<-quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
		names(r)<-c("ymin", "lower", "middle", "upper", "ymax")
		r
	}
	m$temp=""
	m=merge(m, get.sic.labels())
	ggplot(m, aes(x=sic.label, y=value)) + stat_summary(fun.data = f, geom="boxplot")+ facet_grid(temp~variable, scales="free")+coord_flip()+opts(legend.position="none", axis.text.x=theme_text(colour="black"))+labs(x="", y="")+scale_y_continuous(breaks=seq(from=0, to = 40, by = 20))
	ggsave(paste(finalOuts, "RcaseStudies/box.png", sep=""),  dpi=300)
}

plot.medians<-function(){
	to.plot=melt(theta.all[,-(2:3)], c("samp", "sic"), variable_name="stat")
	plot.metric.medians(to.plot, bound=c(0, .7), file="theta/theta_medians.png", tick.length=.5)
	to.plot=melt(m.all[, c("id.samp", "sic", "u2")], c("id.samp", "sic"))
	names(to.plot)=c("samp", "sic", "stat", "value")
	plot.metric.medians(to.plot, bound=c(-0, 1.5), file="u_and_s/u_median.png", tick.length=.25)
	to.plot=melt(m.all[, c("id.samp", "sic", "s2")], c("id.samp", "sic"))
	names(to.plot)=c("samp", "sic", "stat", "value")
	plot.metric.medians(to.plot, bound=c(-.5, .5), file="u_and_s/s_median.png", tick.length=.25)
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
	t$u=pmin(pmax(t$median+1.64*t$se, bound[1]), bound[2])
	t$l=pmin(pmax(t$median-1.64*t$se, bound[1]), bound[2])
	t$median=pmin(pmax(t$median, bound[1]), bound[2])
	t$sig[t$l>h.0]<-"Significant"
	t$sig[t$l<=h.0]<-"Insignificant"
	t$sig=factor(t$sig, levels=c("Significant", "Insignificant"))
	t$temp=""
	p<-ggplot(t, aes(colour=sig, y=median, x=sic.label))
	limits<-aes(ymax=u, ymin=l)
	p + geom_point() + geom_errorbar(limits, width=0.7)+ facet_grid(temp~stat, labeller = label_parsed, scales="free_y")+scale_y_continuous(breaks=seq(from=bound[1], to = bound[2], by = tick.length))+labs(x="", y="")+opts(panel.background=theme_rect(colour="lightgrey", size=1), axis.line=theme_segment(colour="black"), legend.title = theme_text(colour = "white"), legend.position = "none", axis.text.x=theme_text(colour="black"), strip.text.x=theme_text(size=13), panel.background = theme_rect(colour=NA))+scale_colour_grey(end=0.7,start=0)+coord_flip(ylim=c(0,bound))
	ggsave(paste(finalOuts, file, sep=""),  dpi=300)	
}

plot.hausman<-function(d){	
	d=as.data.frame(d)
	d$R2=apply(A.0, 4, pull.R2)
	d$p=pchisq(19, d[,2])
	d$p=pmax(10^(-5), d$p);
	ggplot(d, aes(x=R2, y=p))+geom_point(size=I(.8))+scale_y_log10(breaks=10^(0:-4)) +opts(axis.text.x=theme_text(colour="black", size=11), axis.text.y=theme_text(colour="black", size=11), axis.title.x=theme_text(size=15), axis.title.y=theme_text(size=15, angle=90), panel.background = theme_rect(colour=NA))+scale_x_continuous(breaks=seq(from=0, to = 1, by = .2))+labs(x=expression(R^2), y="p-value")
	ggsave(paste(finalOuts, "hausman/haus.png", sep=""),  dpi=300)
}

plot.theta.mle<-function(d){
	x=seq(.0005, .95, .001)
	u=1; 
	d=d[d$samp==0,-which(names(d)=="samp")]
	d=do.call(rbind, lapply(seq(dim(d)[1]), function(i) cbind(d[i, 1:2], 1-plnorm(.2, d[i, "m"], d[i, "s"]), x, plnorm(x, d[i, "m"], d[i, "s"]))))
	names(d)<-c("biz", "theta", "cdf20", "x", "y")
	d$theta=factor(d$theta)
	levels(d$theta)<-c(expression(theta[0]), expression(theta[1]), expression(theta[2]))
	d$biz=factor(d$biz)
	levels(d$biz)<-c("Retail", "Wholesale", "Manufacturing")
	p<-ggplot(d, aes(fill=cdf20, y=y, x=x))
	limits<-aes(ymax=u)
	p + geom_area()+ facet_grid(theta~biz, labeller = label_parsed, scales="free_y")+coord_cartesian(xlim=c(0,.5), ylim=c(0,1))+labs(x="", y="")+opts(panel.background = theme_rect(colour=NA), axis.text.y=theme_text(size=10, face="bold"), axis.text.x=theme_text(size = 11, face="bold"), strip.text.y=theme_text(size = 15), strip.text.x=theme_text(size = 15), axis.title.x=theme_text(size = 16), axis.title.y=theme_text(size = 16, angle=90))+scale_y_continuous(breaks=seq(.2,1,.2))+scale_x_continuous(breaks=seq(.1,.4,.1))+scale_fill_gradient(high="black", low="lightgrey")+labs(fill="Fraction > .2") 
	ggsave(paste(finalOuts, "theta/theta_MLE_plot.png", sep=""),  dpi=300)
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