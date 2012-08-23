boxcoxnc <-
function(data,method="all",lam=seq(-2,2,0.01),plotit = TRUE,rep=30)
{


library(nortest)
library(tseries)
library(MASS)

data<-as.numeric(data)

if (is.na(min(data))==TRUE) stop("Data include NA")
if (min(data)<=0) stop("Data must include positive values")


if (method=="all") {


# Shapiro Wilk - in nortest
sw<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sw<-rbind(sw,c(lam[i],shapiro.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sw<-rbind(sw,c(lam[i],shapiro.test(log(data))$statistic))
}
swlam<-sw[which.max(sw[,2]),1]
data.sw<-((data^swlam)-1)/swlam
sw.pvalue_data.sw<-shapiro.test(data.sw)$p.value
ad.pvalue_data.sw<-ad.test(data.sw)$p.value
cvm.pvalue_data.sw<-cvm.test(data.sw)$p.value
pt.pvalue_data.sw<-pearson.test(data.sw)$p.value
sf.pvalue_data.sw<-sf.test(data.sw)$p.value
lt.pvalue_data.sw<-lillie.test(data.sw)$p.value
jb.pvalue_data.sw<-jarque.bera.test(data.sw)$p.value

# Anderson Darling - in nortest
ad<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) ad<-rbind(ad,c(lam[i],ad.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) ad<-rbind(ad,c(lam[i],ad.test(log(data))$statistic))
}
adlam<-ad[which.min(ad[,2]),1]

data.ad<-((data^adlam)-1)/adlam
sw.pvalue_data.ad<-shapiro.test(data.ad)$p.value
ad.pvalue_data.ad<-ad.test(data.ad)$p.value
cvm.pvalue_data.ad<-cvm.test(data.ad)$p.value
pt.pvalue_data.ad<-pearson.test(data.ad)$p.value
sf.pvalue_data.ad<-sf.test(data.ad)$p.value
lt.pvalue_data.ad<-lillie.test(data.ad)$p.value
jb.pvalue_data.ad<-jarque.bera.test(data.ad)$p.value

# Cramer Von-Mises - in nortest
cvm<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) cvm<-rbind(cvm,c(lam[i],cvm.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) cvm<-rbind(cvm,c(lam[i],cvm.test(log(data))$statistic))
}
cvmlam<-cvm[which.min(cvm[,2]),1]

data.cvm<-((data^cvmlam)-1)/cvmlam
sw.pvalue_data.cvm<-shapiro.test(data.cvm)$p.value
ad.pvalue_data.cvm<-ad.test(data.cvm)$p.value
cvm.pvalue_data.cvm<-cvm.test(data.cvm)$p.value
pt.pvalue_data.cvm<-pearson.test(data.cvm)$p.value
sf.pvalue_data.cvm<-sf.test(data.cvm)$p.value
lt.pvalue_data.cvm<-lillie.test(data.cvm)$p.value
jb.pvalue_data.cvm<-jarque.bera.test(data.cvm)$p.value

# Pearson Test - in nortest
pt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) pt<-rbind(pt,c(lam[i],pearson.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) pt<-rbind(pt,c(lam[i],pearson.test(log(data))$statistic))
}
ptlam<-pt[which.min(pt[,2]),1]

data.pt<-((data^ptlam)-1)/ptlam
sw.pvalue_data.pt<-shapiro.test(data.pt)$p.value
ad.pvalue_data.pt<-ad.test(data.pt)$p.value
cvm.pvalue_data.pt<-cvm.test(data.pt)$p.value
pt.pvalue_data.pt<-pearson.test(data.pt)$p.value
sf.pvalue_data.pt<-sf.test(data.pt)$p.value
lt.pvalue_data.pt<-lillie.test(data.pt)$p.value
jb.pvalue_data.pt<-jarque.bera.test(data.pt)$p.value

#Shapiro Francia - in nortest
sf<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sf<-rbind(sf,c(lam[i],sf.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sf<-rbind(sf,c(lam[i],sf.test(log(data))$statistic))
}
sflam<-sf[which.max(sf[,2]),1]

data.sf<-((data^sflam)-1)/sflam
sw.pvalue_data.sf<-shapiro.test(data.sf)$p.value
ad.pvalue_data.sf<-ad.test(data.sf)$p.value
cvm.pvalue_data.sf<-cvm.test(data.sf)$p.value
pt.pvalue_data.sf<-pearson.test(data.sf)$p.value
sf.pvalue_data.sf<-sf.test(data.sf)$p.value
lt.pvalue_data.sf<-lillie.test(data.sf)$p.value
jb.pvalue_data.sf<-jarque.bera.test(data.sf)$p.value

# Lilliefors - Kolmogorov Smirnov - in nortest
lt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) lt<-rbind(lt,c(lam[i],lillie.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) lt<-rbind(lt,c(lam[i],lillie.test(log(data))$statistic))
}
ltlam<-lt[which.min(lt[,2]),1]

data.lt<-((data^ltlam)-1)/ltlam
sw.pvalue_data.lt<-shapiro.test(data.lt)$p.value
ad.pvalue_data.lt<-ad.test(data.lt)$p.value
cvm.pvalue_data.lt<-cvm.test(data.lt)$p.value
pt.pvalue_data.lt<-pearson.test(data.lt)$p.value
sf.pvalue_data.lt<-sf.test(data.lt)$p.value
lt.pvalue_data.lt<-lillie.test(data.lt)$p.value
jb.pvalue_data.lt<-jarque.bera.test(data.lt)$p.value

# Jarque-Bera in tseries
jb<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) jb<-rbind(jb,c(lam[i],jarque.bera.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) jb<-rbind(jb,c(lam[i],jarque.bera.test(log(data))$statistic))
}
jblam<-jb[which.min(jb[,2]),1]

data.jb<-((data^jblam)-1)/jblam
sw.pvalue_data.jb<-shapiro.test(data.jb)$p.value
ad.pvalue_data.jb<-ad.test(data.jb)$p.value
cvm.pvalue_data.jb<-cvm.test(data.jb)$p.value
pt.pvalue_data.jb<-pearson.test(data.jb)$p.value
sf.pvalue_data.jb<-sf.test(data.jb)$p.value
lt.pvalue_data.jb<-lillie.test(data.jb)$p.value
jb.pvalue_data.jb<-jarque.bera.test(data.jb)$p.value

aclam1<-NULL
for (q in 1:rep) {
ac<-rnorm(length(data),0,100)
lm1<-glm(data~ac,family=gaussian)

bc1<-boxcox(lm1,lam,plotit=FALSE)
aclam<-bc1$x[which.max(bc1$y)]
aclam1<-cbind(aclam1,aclam)
}



if (plotit==TRUE){
par(mfrow=c(2,4))
plot(sw[,1],sw[,2],ylab="test statistic",xlab=expression(lambda), main="Shapiro-Wilk")
abline(v=swlam,lty=2)
plot(ad[,1],ad[,2],ylab="test statistic",xlab=expression(lambda),main="Anderson-Darling")
abline(v=adlam,lty=2)
plot(cvm[,1],cvm[,2],ylab="test statistic",xlab=expression(lambda),main="Cramer-von Mises")
abline(v=cvmlam,lty=2)
plot(pt[,1],pt[,2],ylab="test statistic",xlab=expression(lambda),main="Pearson Chi-Square")
abline(v=ptlam,lty=2)
plot(sf[,1],sf[,2],ylab="test statistic",xlab=expression(lambda),main="Shapiro-Francia")
abline(v=sflam,lty=2)
plot(lt[,1],lt[,2],ylab="test statistic",xlab=expression(lambda),main="Lilliefors")
abline(v=ltlam,lty=2)
plot(jb[,1],jb[,2],ylab="test statistic",xlab=expression(lambda),main="Jarque-Bera")
abline(v=jblam,lty=2)
boxcox(lm1,lam,plotit)

if (plotit==TRUE){
title("Artificial Covariate")
}

}



aclam1<-as.numeric(aclam1)
aclam2<-mean(aclam1)

data.ac<-((data^aclam2)-1)/aclam2
sw.pvalue_data.ac<-shapiro.test(data.ac)$p.value
ad.pvalue_data.ac<-ad.test(data.ac)$p.value
cvm.pvalue_data.ac<-cvm.test(data.ac)$p.value
pt.pvalue_data.ac<-pearson.test(data.ac)$p.value
sf.pvalue_data.ac<-sf.test(data.ac)$p.value
lt.pvalue_data.ac<-lillie.test(data.ac)$p.value
jb.pvalue_data.ac<-jarque.bera.test(data.ac)$p.value


result<-matrix(0,8,8)
colnames(result)<-c("sw","ad","cvm","pt","sf","lt","jb","ac")
rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")

result[1,]<-c(swlam,adlam,cvmlam,ptlam,sflam,ltlam,jblam,aclam2)
result[2,]<-c(sw.pvalue_data.sw,sw.pvalue_data.ad,sw.pvalue_data.cvm,sw.pvalue_data.pt,sw.pvalue_data.sf,sw.pvalue_data.lt,sw.pvalue_data.jb,sw.pvalue_data.ac)
result[3,]<-c(ad.pvalue_data.sw,ad.pvalue_data.ad,ad.pvalue_data.cvm,ad.pvalue_data.pt,ad.pvalue_data.sf,ad.pvalue_data.lt,ad.pvalue_data.jb,ad.pvalue_data.ac)
result[4,]<-c(cvm.pvalue_data.sw,cvm.pvalue_data.ad,cvm.pvalue_data.cvm,cvm.pvalue_data.pt,cvm.pvalue_data.sf,cvm.pvalue_data.lt,cvm.pvalue_data.jb,cvm.pvalue_data.ac)
result[5,]<-c(pt.pvalue_data.sw,pt.pvalue_data.ad,pt.pvalue_data.cvm,pt.pvalue_data.pt,pt.pvalue_data.sf,pt.pvalue_data.lt,pt.pvalue_data.jb,pt.pvalue_data.ac)
result[6,]<-c(sf.pvalue_data.sw,sf.pvalue_data.ad,sf.pvalue_data.cvm,sf.pvalue_data.pt,sf.pvalue_data.sf,sf.pvalue_data.lt,sf.pvalue_data.jb,sf.pvalue_data.ac)
result[7,]<-c(lt.pvalue_data.sw,lt.pvalue_data.ad,lt.pvalue_data.cvm,lt.pvalue_data.pt,lt.pvalue_data.sf,lt.pvalue_data.lt,lt.pvalue_data.jb,lt.pvalue_data.ac)
result[8,]<-c(jb.pvalue_data.sw,jb.pvalue_data.ad,jb.pvalue_data.cvm,jb.pvalue_data.pt,jb.pvalue_data.sf,jb.pvalue_data.lt,jb.pvalue_data.jb,jb.pvalue_data.ac)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="All"
out$date<-date()
out$result<-result
out
}

else if (method=="sw") {
sw<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sw<-rbind(sw,c(lam[i],shapiro.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sw<-rbind(sw,c(lam[i],shapiro.test(log(data))$statistic))
}

swlam<-sw[which.max(sw[,2]),1]

if (plotit==TRUE){
plot(sw[,1],sw[,2],ylab="test statistic",xlab=expression(lambda), main="Shapiro-Wilk")
abline(v=swlam,lty=2)
}


data.sw<-((data^swlam)-1)/swlam
sw.pvalue_data.sw<-shapiro.test(data.sw)$p.value
ad.pvalue_data.sw<-ad.test(data.sw)$p.value
cvm.pvalue_data.sw<-cvm.test(data.sw)$p.value
pt.pvalue_data.sw<-pearson.test(data.sw)$p.value
sf.pvalue_data.sw<-sf.test(data.sw)$p.value
lt.pvalue_data.sw<-lillie.test(data.sw)$p.value
jb.pvalue_data.sw<-jarque.bera.test(data.sw)$p.value


result<-matrix(0,8,1)
colnames(result)<-c("sw");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(swlam,sw.pvalue_data.sw,ad.pvalue_data.sw,cvm.pvalue_data.sw,pt.pvalue_data.sw,sf.pvalue_data.sw,lt.pvalue_data.sw,jb.pvalue_data.sw)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Shapiro-Wilk"
out$date<-date()
out$result<-result
out
}

else if (method=="ad") {
library(nortest)
ad<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) ad<-rbind(ad,c(lam[i],ad.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) ad<-rbind(ad,c(lam[i],ad.test(log(data))$statistic))
}

adlam<-ad[which.min(ad[,2]),1]
adstat<-ad[which.min(ad[,2]),2]

if (plotit==TRUE){
plot(ad[,1],ad[,2],ylab="test statistic",xlab=expression(lambda),main="Anderson-Darling")
abline(v=adlam,lty=2)
}



data.ad<-((data^adlam)-1)/adlam
sw.pvalue_data.ad<-shapiro.test(data.ad)$p.value
ad.pvalue_data.ad<-ad.test(data.ad)$p.value
cvm.pvalue_data.ad<-cvm.test(data.ad)$p.value
pt.pvalue_data.ad<-pearson.test(data.ad)$p.value
sf.pvalue_data.ad<-sf.test(data.ad)$p.value
lt.pvalue_data.ad<-lillie.test(data.ad)$p.value
jb.pvalue_data.ad<-jarque.bera.test(data.ad)$p.value


result<-matrix(0,8,1)
colnames(result)<-c("ad");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(adlam,sw.pvalue_data.ad,ad.pvalue_data.ad,cvm.pvalue_data.ad,pt.pvalue_data.ad,sf.pvalue_data.ad,lt.pvalue_data.ad,jb.pvalue_data.ad)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Anderson-Darling"
out$date<-date()
out$result<-result
out
}

else if (method=="cvm") {
library(nortest)
# Cramer Von-Mises - in nortest
cvm<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) cvm<-rbind(cvm,c(lam[i],cvm.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) cvm<-rbind(cvm,c(lam[i],cvm.test(log(data))$statistic))
}

cvmlam<-cvm[which.min(cvm[,2]),1]
cvmstat<-cvm[which.min(cvm[,2]),2]


if (plotit==TRUE){
plot(cvm[,1],cvm[,2],ylab="test statistic",xlab=expression(lambda),main="Cramer-von Mises")
abline(v=cvmlam,lty=2)
}


data.cvm<-((data^cvmlam)-1)/cvmlam
sw.pvalue_data.cvm<-shapiro.test(data.cvm)$p.value
ad.pvalue_data.cvm<-ad.test(data.cvm)$p.value
cvm.pvalue_data.cvm<-cvm.test(data.cvm)$p.value
pt.pvalue_data.cvm<-pearson.test(data.cvm)$p.value
sf.pvalue_data.cvm<-sf.test(data.cvm)$p.value
lt.pvalue_data.cvm<-lillie.test(data.cvm)$p.value
jb.pvalue_data.cvm<-jarque.bera.test(data.cvm)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("cvm");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(cvmlam,sw.pvalue_data.cvm,ad.pvalue_data.cvm,cvm.pvalue_data.cvm,pt.pvalue_data.cvm,sf.pvalue_data.cvm,lt.pvalue_data.cvm,jb.pvalue_data.cvm)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Cramer-von Mises"
out$date<-date()
out$result<-result
out
}

else if (method=="pt") {
library(nortest)
pt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) pt<-rbind(pt,c(lam[i],pearson.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) pt<-rbind(pt,c(lam[i],pearson.test(log(data))$statistic))
}

ptlam<-pt[which.min(pt[,2]),1]
ptstat<-pt[which.min(pt[,2]),2]

if (plotit==TRUE){
plot(pt[,1],pt[,2],ylab="test statistic",xlab=expression(lambda),main="Pearson Chi-Square")
abline(v=ptlam,lty=2)
}



data.pt<-((data^ptlam)-1)/ptlam
sw.pvalue_data.pt<-shapiro.test(data.pt)$p.value
ad.pvalue_data.pt<-ad.test(data.pt)$p.value
cvm.pvalue_data.pt<-cvm.test(data.pt)$p.value
pt.pvalue_data.pt<-pearson.test(data.pt)$p.value
sf.pvalue_data.pt<-sf.test(data.pt)$p.value
lt.pvalue_data.pt<-lillie.test(data.pt)$p.value
jb.pvalue_data.pt<-jarque.bera.test(data.pt)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("pt");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(ptlam,sw.pvalue_data.pt,ad.pvalue_data.pt,cvm.pvalue_data.pt,pt.pvalue_data.pt,sf.pvalue_data.pt,lt.pvalue_data.pt,jb.pvalue_data.pt)



out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Pearson Chi-Square"
out$date<-date()
out$result<-result
out
}

else if (method=="sf") {
library(nortest)
sf<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) sf<-rbind(sf,c(lam[i],sf.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) sf<-rbind(sf,c(lam[i],sf.test(log(data))$statistic))
}

sflam<-sf[which.max(sf[,2]),1]
sfstat<-sf[which.max(sf[,2]),2]

if (plotit==TRUE){
plot(sf[,1],sf[,2],ylab="test statistic",xlab=expression(lambda),main="Shapiro-Francia")
abline(v=sflam,lty=2)
}


data.sf<-((data^sflam)-1)/sflam
sw.pvalue_data.sf<-shapiro.test(data.sf)$p.value
ad.pvalue_data.sf<-ad.test(data.sf)$p.value
cvm.pvalue_data.sf<-cvm.test(data.sf)$p.value
pt.pvalue_data.sf<-pearson.test(data.sf)$p.value
sf.pvalue_data.sf<-sf.test(data.sf)$p.value
lt.pvalue_data.sf<-lillie.test(data.sf)$p.value
jb.pvalue_data.sf<-jarque.bera.test(data.sf)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("sf");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(sflam,sw.pvalue_data.sf,ad.pvalue_data.sf,cvm.pvalue_data.sf,pt.pvalue_data.sf,sf.pvalue_data.sf,lt.pvalue_data.sf,jb.pvalue_data.sf)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Shapiro-Francia"
out$date<-date()
out$result<-result
out
}

else if (method=="lt") {
library(nortest)
lt<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) lt<-rbind(lt,c(lam[i],lillie.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) lt<-rbind(lt,c(lam[i],lillie.test(log(data))$statistic))
}

ltlam<-lt[which.min(lt[,2]),1]
ltstat<-lt[which.min(lt[,2]),2]

if (plotit==TRUE){
plot(lt[,1],lt[,2],ylab="test statistic",xlab=expression(lambda),main="Lilliefors")
abline(v=ltlam,lty=2)
}


data.lt<-((data^ltlam)-1)/ltlam
sw.pvalue_data.lt<-shapiro.test(data.lt)$p.value
ad.pvalue_data.lt<-ad.test(data.lt)$p.value
cvm.pvalue_data.lt<-cvm.test(data.lt)$p.value
pt.pvalue_data.lt<-pearson.test(data.lt)$p.value
sf.pvalue_data.lt<-sf.test(data.lt)$p.value
lt.pvalue_data.lt<-lillie.test(data.lt)$p.value
jb.pvalue_data.lt<-jarque.bera.test(data.lt)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("lt");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(ltlam,sw.pvalue_data.lt,ad.pvalue_data.lt,cvm.pvalue_data.lt,pt.pvalue_data.lt,sf.pvalue_data.lt,lt.pvalue_data.lt,jb.pvalue_data.lt)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Lilliefors"
out$date<-date()
out$result<-result
out
}

else if (method=="jb") {
library(tseries)
jb<-NULL
for (i in 1:length(lam)){
if (round(lam[i],2)!=0) jb<-rbind(jb,c(lam[i],jarque.bera.test((data**(lam[i])-1)/(lam[i]))$statistic))
if (round(lam[i],2)==0) jb<-rbind(jb,c(lam[i],jarque.bera.test(log(data))$statistic))
}


jblam<-jb[which.min(jb[,2]),1]
jbstat<-jb[which.min(jb[,2]),2]

if (plotit==TRUE){
plot(jb[,1],jb[,2],ylab="test statistic",xlab=expression(lambda),main="Jarque-Bera")
abline(v=jblam,lty=2)
}



data.jb<-((data^jblam)-1)/jblam
sw.pvalue_data.jb<-shapiro.test(data.jb)$p.value
ad.pvalue_data.jb<-ad.test(data.jb)$p.value
cvm.pvalue_data.jb<-cvm.test(data.jb)$p.value
pt.pvalue_data.jb<-pearson.test(data.jb)$p.value
sf.pvalue_data.jb<-sf.test(data.jb)$p.value
lt.pvalue_data.jb<-lillie.test(data.jb)$p.value
jb.pvalue_data.jb<-jarque.bera.test(data.jb)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("jb");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(jblam,sw.pvalue_data.jb,ad.pvalue_data.jb,cvm.pvalue_data.jb,pt.pvalue_data.jb,sf.pvalue_data.jb,lt.pvalue_data.jb,jb.pvalue_data.jb)


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Jarque-Bera"
out$date<-date()
out$result<-result
out
}

else if (method=="ac") {
aclam1<-NULL

for (t in 1:rep) {
library(MASS)
ac<-rnorm(length(data),0,100)
lm1<-glm(data~ac,family=gaussian)
bc1<-boxcox(lm1,lam,plotit=FALSE)
aclam<-bc1$x[which.max(bc1$y)]
aclam1<-cbind(aclam1,aclam)
}

boxcox(lm1,lam,plotit)
if (plotit==TRUE){
title("Artificial Covariate")
}
aclam1<-as.numeric(aclam1)
aclam2<-mean(aclam1)

data.ac<-((data^aclam2)-1)/aclam2
sw.pvalue_data.ac<-shapiro.test(data.ac)$p.value
ad.pvalue_data.ac<-ad.test(data.ac)$p.value
cvm.pvalue_data.ac<-cvm.test(data.ac)$p.value
pt.pvalue_data.ac<-pearson.test(data.ac)$p.value
sf.pvalue_data.ac<-sf.test(data.ac)$p.value
lt.pvalue_data.ac<-lillie.test(data.ac)$p.value
jb.pvalue_data.ac<-jarque.bera.test(data.ac)$p.value

result<-matrix(0,8,1)
colnames(result)<-c("ac");rownames(result)<-c("lambda.hat","sw.pvalue","ad.pvalue","cvm.pvalue","pt.pvalue","sf.pvalue","lt.pvalue","jb.pvalue")
result[,1]<-c(aclam2,sw.pvalue_data.ac,ad.pvalue_data.ac,cvm.pvalue_data.ac,pt.pvalue_data.ac,sf.pvalue_data.ac,lt.pvalue_data.ac,jb.pvalue_data.ac)

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$version<-"Version 1.0"
out$method="Artifical Covariate"
out$date<-date()
out$result<-result
out
}

}
