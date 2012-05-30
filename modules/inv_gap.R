inv.gap<-function(){
	D.gap=data.frame(read.table(paste(dataSource, 'invDat.txt', sep=""), header=TRUE))
	coeffs=by(D.gap, D.gap$gvkey, lm, formula=inv+e0~c_f1+inv1)
	coeffs=as.data.frame(t(vapply(coeffs, function(x) x$coefficients, rep(0,3))))
	coeffs$gvkey=row.names(coeffs)
	bizlist=t(vapply(unique(apply(D.gap[,c('gvkey', 'biz')], 1, list)), unlist, c(0,0)))
	coeffs=merge(coeffs, bizlist)
	coeffs
}
