#Simulate some data
E=rbind(matrix(rnorm(500), 5), matrix(0, 4, 100))
ot=laply(0:100, function(l){sum(diag(t(makeD(100,l))%*%t(E[nrow(E):1, ])))})

E.Delta=(diag(9)-makeD(9,1))%*%E
Delta.ot=laply(0:100, function(l){sum(diag(t(makeD(100,l))%*%t(E.Delta[nrow(E.Delta):1, ])))})

E.Delta2=(diag(9)-makeD(9,1))%*%(diag(9)-makeD(9,1))%*%E
Delta2.ot=laply(0:100, function(l){sum(diag(t(makeD(100,l))%*%t(E.Delta2[nrow(E.Delta2):1, ])))})

direct.Delta.ot=ot-makeD(101,1)%*%ot
direct.Delta2.ot=Delta.ot-makeD(101,1)%*%Delta.ot
