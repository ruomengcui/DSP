run.monte.carlo<-function(){
  which.alpha=c(0, 1); which.beta=c(0, 3)
  H=16; z.scale=100; S=100; num.boot=30
  trials=c(100, 1000)
	l_ply(trials, function(my.N){
	  pn=data.name=paste('simPanel', my.N, '.txt', sep="")
	  en=paste('simEst', my.N, '.txt', sep="")
    create.sim.panel(data.name=pn, which.alpha=which.alpha, which.beta=which.beta, N=my.N, S=S, num.boot=num.boot, z.scale=z.scale, H=H)
	  calc.parameters(data.in=pn, data.out=en, which.alpha=which.alpha, which.beta=which.beta, mc=T)
	})	
  plot.mc(trials, 'simGraph', which.alpha, which.beta)
}