makeC <- function(N){
	C=matrix(0, N, N)
	C[lower.tri(C, diag = TRUE)] <-1
	C
}

makeD <- function(N, L){
	x=rbind(matrix(0, abs(L), N), diag(1, N-abs(L), N))
	if(L>=0) x
	else t(x)
}

makeK<-function(H, L){
	o=matrix(0, H+1, H+1)
	o[H+1-L,]<-1
	o
}

makeG<-function(H, L){
	rbind(diag(H-L), rep(-1, H-L), matrix(0, L, H-L))
}

make.S<-function(H, j){
	((diag(H+1+j)-makeD(H+1+j,1))%^%j)[,seq(H+1)]
}

makeIj <- function(H, L){diag(c(array(1,L),array(0,H-L)))}

makeR <- function(H, L){R=t(makeD(H, L)); R[1,seq(L)]=1; R}

make.el<-function(H, l){ x=rep(0, H); x[l]=1; x}