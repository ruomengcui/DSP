counter.factual<-function(my.firm, metric.fun, pref.fun=identity, S.fun=identity, Lambda.fun=identity){
#takes four functions as arguments: first is metric function---which takes A, S, and Lambda as arguments---and last three tell how each parameter is updated
	attach(my.firm)
	print(id)
	o1=metric.fun(A[,,2], S, Lambda[,,2])
	o2=metric.fun(solve.for.A(pref.fun(alpha.beta), S.fun(S), Lambda.fun(Lambda[,,2]))$A, S.fun(S), Lambda.fun(Lambda[,,2]))
	o=list(id=id, o1=o1, o2=o2)
	detach(my.firm)
	o
}

m.bw<-function(A, S, Lambda){
	x=A%*%S%*%t(A)+Lambda-S
	c(sum(diag(x)), x[1,1])
}

m.ASA<-function(A, S, Lambda) diag(A%*%S%*%t(A)) 

improve.forecast<-function(S){
	D=makeD(dim(S)[1], 1)
	D%*%S%*%t(D)
}