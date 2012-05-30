get.fourier.matrix<-function(N){
	i=complex(real=0, imaginary=1)
	w=exp(-2*pi*i/N)
	N^(-1/2)*sapply(seq(N)-1, function(x) (w^(seq(N)-1))^x)
} 

get.Gamma<-function(my.firm){
	attach(my.firm)
	o=list(id=id, Gamma=diag(F%*%S%*%Conj(F)), Gammao=diag(F%*%(A[,,1]%*%S%*%t(A[,,1])+Lambda[,,1])%*%Conj(F)))
	detach(my.firm)
	o
}