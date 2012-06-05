create.compustat.panel<-function(){
	if(load.data){
		D1=data.frame(read.table(paste(dataSource, 'toReg1', d.tag, '.txt', sep=""), header=TRUE))
		D2=data.frame(read.table(paste(dataSource, 'sigmacov1', d.tag, '.txt', sep=""), header=TRUE))
		key=lapply(unique(apply(D1[,c(1,2)], 1, list)), unlist)
	}
	
	x=2*(1:H)+1; y=x+1; 
	A=array(dim=c(H, H, 2))	
	Lambda=A
		
	r1=(D1$samp==id[1]) & (D1$key==id[2])
	r2=(D2$samp==id[1]) & (D2$key==id[2])
	S.cols=substring(names(D2), 2, 2)%in%as.character(0:9)
	E=data.matrix(D1[r1, x])
	Eo=data.matrix(D1[r1, y]) 
	S=matrix(as.real(D2[r2,S.cols]), sqrt(sum(S.cols)));
}

create.wards.panel<-function(){
	#base panel
	d=llply(c("sales", "prod"), function(var){	
		x=rbind(melt(read.csv(paste(dataSource, "Wards/", var, "1.csv", sep=""), header = TRUE), 1:4), melt(read.csv(paste(dataSource, "Wards/", var, "2.csv", sep=""), header = TRUE), 1:4))
		x$value=as.numeric(gsub(",","", as.character(x$value)))
		x$year=laply(as.character(x$variable), function(y){as.integer(substr(y,nchar(y)-3,nchar(y)))})
		x$month=as.integer(substr(x$variable, 2,3))
		x$quart=x$year+(x$month-1)/12
		names(x)[which(names(x)=="value")]<-var
		x[c(1:4, 9, 6)]
	})
	d=merge(d[[1]], d[[2]])
	d=d[(d$Company != "")&(d$sales>30)&(d$prod>30),]
	unique.cars=unique(d[1:3])
	unique.cars$cid=1:dim(unique.cars)[1]
	d=merge(d, unique.cars)	
	d=d[order(d$cid, d$quart),]
	d=block.boot(d, 40, 36)
	d=d[,c(8, 10, 9, 5,  1:4, 6, 7)]
	
	#create signals
	d<-ddply(d, c("samp", "cid"), function(s){
		#Remove firms that violate regularity conditions
		if((length(unique(s$quart))<num.periods)||sum(s$prod)>3/2*sum(s$sales)||sum(s$sales)>3/2*sum(s$prod)) return()
			
		#Remove quadratic trend, seperate out seasonal shifts, and add detrended and seasonal series back to into s
		s=s[order(s$nboot),]
		detrend=aperm(aaply(s[,c("sales", "prod")], 2, function(v){
			quatratic.trend=lm(as.matrix(v)~seq(dim(s)[1])+seq(dim(s)[1])^2)$residuals
			seasonality=lm(as.matrix(quatratic.trend)~as.factor(ceiling(12*(s$quart-floor(s$quart))+1)))
			cbind(seasonality$residuals, seasonality$fitted.values)
		}), c(2,1,3))
		append=c(".dt", ".season")
		for(i in seq(append)){
			x=as.data.frame(detrend[,,i])
			names(x)<-paste(names(x), append[i], sep="")
			s=cbind(s, x)
		}
		
		#Run regressions and get signals; Use X*beta, rather than the fitted values of the regression, because this way we get more forecasts
		t.vars=ts(cbind(detrend[,,1], detrend[,,1]^2))
		reg.panel=ts(s[,c("nboot", "sales.dt", "prod.dt")])
		for(i in seq(H+p-1)) reg.panel=ts.union(reg.panel, lag(t.vars, -i))
		Y.hat=aperm(aaply(reg.panel[,2:3], 2, function(Y){
			laply(seq(H), function(h){
				X=reg.panel[,3+(h-1)*dim(t.vars)[2]+seq(p*dim(t.vars)[2])]
				beta=lm(as.matrix(Y, length(Y))~X-1)$coefficients
				X%*%beta
			})
		}), 3:1)
		
		#Difference forecasts into signals; align signals into epsilon vectors; then add these vectors back to s.		
		for(i in 1:2){
			forecasts=cbind(reg.panel[,1+i], Y.hat[,,i], c(rep(NA, H+p), rep(0, dim(Y.hat)[1]-H-p)))
			signals=-t(diff(t(forecasts)))
			epsilon=as.data.frame(rbind(matrix(rep(NA, p*(H+1)), p), matrix(signals[!is.na(signals)], sum(!is.na(signals[,1])))))
			names(epsilon)<-paste(c("eps", "eps^o")[i], 0:H, sep="_")
			s=cbind(s, epsilon)
		}
		s
	})
	save(d, file=paste(varSave, 'wards_signal_panel.txt', sep=''))
}

create.sim.panel<-function(H=6, N=2000, S=20, which.alpha=c(0, 1), which.beta=c(0, 3), data.name='sim_signal_panel.txt', num.boot=30, positive.Gamma=T, z.scale=1){
	alpha.max=1
  beta.max=1
  d=ldply(seq(S), function(cid){
		alpha=alpha.max*runif(length(which.alpha))
		beta=beta.max*runif(length(which.beta))
		block.mat=rwishart(2*(H+1), diag(2*(H+1)))$W
    Omega=rwishart(H+1, z.scale*diag(H+1))$W
		signals=mvrnorm(N, rep(0, 2*(H+1)), block.mat)
		z=mvrnorm(N, rep(0, H+1), Omega)
		eps=positive.Gamma*signals[,seq(H+1)]+z
		eta=signals[,H+1+seq(H+1)]
    Sigma=block.mat[seq(H+1),seq(H+1)]+Omega
    Gamma=block.mat[seq(H+1),H+1+seq(H+1)]
    Lambda=block.mat[H+1+seq(H+1),H+1+seq(H+1)]
    A=solve.for.A(list(alpha=alpha, which.alpha=which.alpha, beta=beta, which.beta=which.beta, Sigma=Sigma, Lambda=Lambda, Gamma=Gamma, L=0))$A
		epso=eps%*%t(A)+eta
		colnames(eps)<-paste("e_", seq(dim(eps)[2])-1, sep="")
		colnames(epso)<-paste("o_", seq(dim(epso)[2])-1, sep="")
		colnames(z)<-paste("z_", seq(dim(z)[2])-1, sep="")
		names(alpha)<-paste("alpha_", which.alpha, sep="")
		names(beta)<-paste("beta_", which.beta, sep="")
		quart=seq(N)
		cbind(cid, quart, merge(merge(rbind(alpha),rbind(beta)), as.data.frame(cbind(eps, epso, z))))
	})
  d=block.boot(d, num.boot, 20)
	save(d, file=paste(varSave, data.name, sep=''))
}

