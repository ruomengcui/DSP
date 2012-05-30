my.scale = 1
gn=paste('sim_graph_', my.scale, '.txt', sep="")
graph.name=gn


a$ab=as.factor(a$ab)
levels(a$ab)<-c(expression(alpha), expression(beta))
ggplot(a,aes(y=e,x=p,colour=c))+geom_point()+geom_errorbar(aes(max=u,ymin=l))+geom_abline(intercept=0,slope=1)+facet_grid(ab~n,labeller=label_parsed)+scale_colour_grey(end=0.7,start=0)+labs(x="Parameter", y="Estimate")+coord_cartesian(ylim=c(-.05,6.3), xlim=c(-.05, 3.05))+scale_x_continuous(breaks=seq(0, 3, .5))+scale_y_continuous(breaks=seq(0, 6, 1))+opts(panel.background = theme_rect(colour=NA), legend.position = "none", strip.text.y=theme_text(size=16), axis.title.x=theme_text(size = 18), axis.title.y=theme_text(size = 18, angle=90), axis.text.y=theme_text(size=13), axis.text.x=theme_text(size = 13))
ggsave(paste(finalOuts, "mc/", graph.name, ".png", sep=""),  dpi=300)	
