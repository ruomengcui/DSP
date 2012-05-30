run.monte.carlo<-function(){
	alpha.max=3
	beta.max=3
	l_ply(c(1, 10, 100), function(my.scale){
		pn=paste('sim_panel_', my.scale, '.txt', sep="") 
		en=paste('sim_est_', my.scale, '.txt', sep="")
		gn=paste('sim_graph_', my.scale, sep="")
		create.sim.panel(data.name=pn, alpha.max=alpha.max, beta.max=beta.max, N=1000, S=100, num.boot=30, positive.Gamma=F, z.scale=my.scale)
		o=calc.parameters(data.name=pn, L=0, mc=T, param.count=c(length(alpha.max),length(beta.max)))
		save(o, file=paste(varSave, en, sep=''))
		plot.mc(pn, en, gn)
	})	
}

