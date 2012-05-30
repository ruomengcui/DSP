plot.Impulse.Responses<-function(A){
	#plot IRF means and quantiles, faceting by estimator and lead time.
	stats=array(dim=c(4, dim(A)[1:3]))
	stats[1,,,]=apply(A, c(1,2,3), mean)
	stats[2:4,,,]=apply(A, c(1,2,3), quantile, probs=c(.25, .5, .75))
	d=data.frame(melt(stats, varnames=c("stat", "m", "n", "est")))
	d$stat<-as.factor(d$stat)
	d$est<-as.factor(d$est)
	d$n<-as.factor(d$n)
	levels(d$stat)<-c('Mean', '1st Quartile', 'Median', '2nd Quartile')
	levels(d$est)<-c('All Feasible', 'Economically Sound', 'Triple Newsvendor')
	levels(d$n)<-c('n=0', 'n=1', 'n=2', 'n=3', 'n=4')
	d$m=d$m-1
	p <- ggplot(d, aes(x=m, y=value, colour=stat)) +geom_line() +facet_grid(n~est)+opts(legend.title = theme_text(size = 0))+ scale_y_continuous(expression("A"["m,n"]))
	ggsave(paste(finalOuts, "impulse_responses/impulse_responses.png", sep=""),  dpi=400)
}

plot.theta<-function(theta){
	t=theta[,3:dim(theta)[2]]
	names(t)<-c("Sector", expression(theta), expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]))
	t$Sector<-factor(t$Sector)
	levels(t$Sector)<-c("Wholesale", "Retail", "Manufacturing")
	t=melt(t, id.var = "Sector")
	t$value=pmin(t$value, 3.75)
	p <- ggplot(t, aes(x=value, fill=variable))+geom_density(adjust=.34)+facet_grid(variable~Sector, scales="free", labeller = label_parsed)+opts(axis.text.x = theme_text(size = 7))+opts(legend.position="none")+ xlab("Parameter Value")+ylab("Density")+coord_cartesian(xlim=c(0,.55))+scale_x_continuous(breaks=seq(from=0, to = 1, by = .1)) 
	ggsave(paste(finalOuts, "theta/theta_density.png", sep=""),  dpi=400)
}

plot.UandSmetrics<-function(A){
	#plot IRF means and quantiles, faceting by estimator and lead time.
}

plot.metrics<-function(m){
	m=melt(m, id.var = c('id.gvkey', 'sic'))
	m=m[m$variable%in%c("adh2", "stock2", "lead02", "inv", "prod"),]	
	m$variable<-factor(m$variable, levels=c("adh2", "stock2", "lead02", "inv", "prod"))
	m$sic=as.factor(floor(m$sic/100))
	q=do.call(rbind, tapply(m$value, m$variable, quantile, probs=c(.1, .9)))
	m$value=pmin(pmax(m$value, q[match(m$variable, rownames(q)), 1]), q[match(m$variable, rownames(q)), 2])
	levels(m$variable)<-c('Adherence to production plans', 'Variability of stock-up-to level', 'Lead-0 bullwhip', 'Expected time until inventory exceeds [a,b]', 'Expected time until production exceeds y')
	f <- function(x) {
		r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
		names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
		r
	}
	ggplot(m, aes(sic, value, fill = sic)) + stat_summary(fun.data = f, geom="boxplot")+ facet_wrap(~variable, ncol = 1, scales="free")+opts(legend.position="none")+ xlab("")+ ylab("") 
	ggsave(paste(finalOuts, "RcaseStudies/box.png", sep=""),  dpi=400)
}