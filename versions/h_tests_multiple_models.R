hausman<-function(){
	pull.gv<-function(s) s$id[2]
	gvs=unique(sapply(f.0, pull.gv))
	hausman.est<-function(est){
		haus<-function(my.firm){
			get.firms<-function(s) (s$id[2]==my.firm)
			boots=firms[sapply(firms, get.firms)]
			real=f.0[sapply(f.0, get.firms)]
			get.est.diff<-function(t) c(t$A[,,2]-t$A[,,est])
			boot.diff=t(vapply(boots, get.est.diff, c(firms[[1]]$A[,,1])))
			real.diff=vapply(real, get.est.diff, c(firms[[1]]$A[,,1])) 
			t<-try(t(real.diff)%*%solve(cov(boot.diff), real.diff))
			if("try-error" %in% class(t)) t=-1
			t
		}
		do.call(rbind, mclapply(gvs, haus, mc.preschedule=FALSE))
	}
	vapply(3:5, hausman.est, as.double(gvs)) 
}

