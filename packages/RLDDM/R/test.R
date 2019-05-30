# # 
vv <- matrix(1:8, nrow=4, ncol=2)
vv[2,1] <- NA
vv[2,2] <- 3.2
# 
# .C('tester', VV=as.double(vv), nrow=nrow(vv), ncol=ncol(vv), NAOK=TRUE)
#