#!/bin/R
library('data.table'); library('mvmeta')

# load data
# use this line if not using chunks -- otherwise see below
# X=fread('data/combined.21.genPC.nominal.maf05.txt')

# control which chunk
i=commandArgs(trailingOnly=TRUE)
if (length(i) == 0) {
	i=1
} else {
	i=i[1]
}
X=fread(paste0('temp/in/af.eqtl.hg38.',i,'.txt'))
print(head(X))
colnames(X)[1:2]=c("Variant","Transcript")

# lil helper guy
meta=function(row) {
	Y=as.matrix(as.numeric(row[3:8])) # 3:8 if AFR, 3:6 if EUR
	SS=as.matrix(as.numeric(row[9:14]))^2 # 9:14 if AFR, 7:10 if EUR
	if (sum(!is.na(Y))==1) {
		ix=which(!is.na(Y))
		SE=sqrt(SS[ix])
		Zz=Y[ix]/SE
		Pp=pnorm(-abs(Zz))
		Ca=qnorm(0.975)*SE
		return(c(Y[ix],SE,Zz,Pp,Y[ix]-Ca,Y[ix]+Ca,'NA',1,'NA'))
	}
	try({ # healthy yolo
		meta=mvmeta(Y~1, S=SS, method='mm') # fixed effects meta-analysis, with q-test
		test=qtest(meta)
		coef=as.matrix(cbind(summary(meta)$coefficients, as.data.frame(unclass(test)[1:3])))
		return(coef)
	})
	return(rep(NA, 9))
}

# munge, write to file -- use bottom line if not using chunks
Z=t(apply(X, 1, meta))
colnames(Z)=c("BETA","SE","ZSTAT","P","95PCT_CI_LOWER", "95PCT_CI_UPPER","Q","DF","P_HET")
write.table(cbind(X[,1:2], Z), file=paste0('temp/out/af.eQTL.mvmeta.hg38.',i,'.txt'), 
             sep="\t", row.names=F, quote=F)
