calc.metrics<-function(iter, f){
	attach(f[[iter]])
	met.funs=define.m(A, Lambda, S, Omega, Xi)
	applier<-function(x) x()
	o=lapply(met.funs[[1]], applier)
	names(o)<-met.funs[[2]]
	o$id=f[[iter]]$id
	detach(f[[iter]])
	o
}

#metric definitions
define.m<-function(A, Lambda, S, Omega, Xi){
	if(length(dim(A))==2) A=array(A, c(dim(A), 1))
	if(length(dim(Lambda))==2) Lambda=array(Lambda, c(dim(Lambda), 1))
	m.stockVar<-function(){
		k=dim(A)[1]
		C=matrix(0, k, k)
		C[lower.tri(C, diag = TRUE)] <-1
		to.calc<-function(i){
			sum(diag(C%*%(A[,,i]-diag(k))%*%S%*% t(A[,,i]-diag(k))%*%t(C)+C%*%Lambda[,,i]%*%t(C)))
		}
		lapply(seq(dim(A)[3]), to.calc)
	}
	m.adherence<-function(){
		to.calc<-function(i){
			diag(A[,,i]%*%S%*%t(A[,,i])+Lambda[,,i])%*%(1/2^(1:dim(A)[1]))
		}
		lapply(seq(dim(A)[3]), to.calc)
	}	
	m.lead0<-function(){
		to.calc<-function(i){
			(A[,,i]%*%S%*%t(A[,,i])+Lambda[,,i]-S)[1,1]
		}
		lapply(seq(dim(A)[3]), to.calc)
	}
	m.prodExceed<-function(){
		a=-99999; b=1.645*Psi[1,1]^(1/2)
		draws=rmnorm(50000, c(0, 0), Xi)
		(1/2*(pnorm((b)* Xi[1,1]^(-1/2))- pnorm((a)*Xi[1,1]^(-1/2)))^2) *(mean(((draws[,1]>a)*(draws[,1]<b)) * ((draws[,2]<a)+(draws[,2]>b)))^(-1))		
	}
	m.invExceed<-function(){
		a=m-Psi[1,1]^(1/2); b=m+Psi[1,1]^(1/2)
		draws=rmnorm(50000, c(m, m), Omega)
		(1/2*(pnorm((b-m)*Omega[1,1]^(-1/2))- pnorm((a-m)*Omega[1,1]^(-1/2)))^2)*(mean(((draws[,1]>a)*(draws[,1]<b)) * ((draws[,2]<a)+(draws[,2]>b)))^(-1))
	}
	m.supplier.Inv<-function(){
		k=dim(A)[1]
		C=makeC(k)
		D=makeD(k, 1)
		R=makeR(k, 2)
		to.calc<-function(i){
			sqrt(sum(diag(C%*%(D%*%R-diag(k))%*%(A[,,i]%*%S%*%t(A[,,i])+Lambda[,,i])%*%t(D%*%R-diag(k))%*%t(C))))
		}
		lapply(seq(dim(A)[3]), to.calc)	
	}
	m.u<-function(){
		to.calc<-function(i) sum(A[,,i][upper.tri(A[,,i], diag = FALSE)])
		lapply(seq(dim(A)[3]), to.calc)
	}
	m.s<-function(){
		to.calc<-function(i) sum(A[,,i][lower.tri(A[,,i], diag = FALSE)]-A[,,i][upper.tri(A[,,i], diag = FALSE)])
		lapply(seq(dim(A)[3]), to.calc)		
	}
	list(
		list(m.stockVar, m.adherence, m.lead0, m.prodExceed, m.invExceed, m.u, m.s, m.supplier.Inv),
		c("stock", "adh", "lead0", "prod", "inv", "u", "s", "sinv")
	)
}
