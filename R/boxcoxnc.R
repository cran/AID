boxcoxnc <-
function(data,method="sw",lam=seq(-3,3,0.01),plotit = TRUE, alpha=0.05, verbose = TRUE)
{
dname<-deparse(substitute(data))

data<-as.numeric(data)



if (is.na(min(data))==TRUE) stop("Data include NA")
if (min(data)<=0) stop("Data must include positive values")




if (method=="sw") {
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


if (swlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (swlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")



if (swlam!=0) data.sw<-((data^swlam)-1)/swlam
if (swlam==0) data.sw<-log(data)


sw.pvalue<-shapiro.test(data.sw)$p.value


result<-matrix(0,2,1)
colnames(result)<-c("sw");rownames(result)<-c("lambda.hat","sw.pvalue")
result[,1]<-c(swlam,sw.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", swlam, "\n", sep = " ")
cat("  p.value    :", sw.pvalue, "\n\n", sep = " ")
cat(if(sw.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Shapiro-Wilk"
out$lambda.hat<-as.numeric(swlam)
out$p.value<-sw.pvalue
out$tf.data<-data.sw
out$date<-date()
invisible(out)
}

else if (method=="ad") {

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


if (adlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (adlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (adlam!=0) data.ad<-((data^adlam)-1)/adlam
if (adlam==0) data.ad<-log(data)

ad.pvalue<-ad.test(data.ad)$p.value



result<-matrix(0,2,1)
colnames(result)<-c("ad");rownames(result)<-c("lambda.hat","ad.pvalue")
result[,1]<-c(adlam,ad.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", adlam, "\n", sep = " ")
cat("  p.value    :", ad.pvalue, "\n\n", sep = " ")
cat(if(ad.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Anderson-Darling"
out$lambda.hat<-as.numeric(adlam)
out$p.value<-ad.pvalue
out$tf.data<-data.ad
out$date<-date()
invisible(out)
}

else if (method=="cvm") {

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


if (cvmlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (cvmlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (cvmlam!=0) data.cvm<-((data^cvmlam)-1)/cvmlam
if (cvmlam==0) data.cvm<-log(data)

cvm.pvalue<-cvm.test(data.cvm)$p.value


result<-matrix(0,2,1)
colnames(result)<-c("cvm");rownames(result)<-c("lambda.hat","cvm.pvalue")
result[,1]<-c(cvmlam,cvm.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", cvmlam, "\n", sep = " ")
cat("  p.value    :", cvm.pvalue, "\n\n", sep = " ")
cat(if(cvm.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Cramer-von Mises"
out$lambda.hat<-as.numeric(cvmlam)
out$p.value<-cvm.pvalue
out$tf.data<-data.cvm
out$date<-date()
invisible(out)
}

else if (method=="pt") {

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


if (ptlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ptlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (ptlam!=0) data.pt<-((data^ptlam)-1)/ptlam
if (ptlam==0) data.pt<-log(data)

pt.pvalue<-pearson.test(data.pt)$p.value


result<-matrix(0,2,1)
colnames(result)<-c("pt");rownames(result)<-c("lambda.hat","pt.pvalue")
result[,1]<-c(ptlam,pt.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", ptlam, "\n", sep = " ")
cat("  p.value    :", pt.pvalue, "\n\n", sep = " ")
cat(if(pt.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Pearson Chi-Square"
out$lambda.hat<-as.numeric(ptlam)
out$p.value<-pt.pvalue
out$tf.data<-data.pt
out$date<-date()
invisible(out)
}

else if (method=="sf") {

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


if (sflam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (sflam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (sflam!=0) data.sf<-((data^sflam)-1)/sflam
if (sflam==0) data.sf<-log(data)


sf.pvalue<-sf.test(data.sf)$p.value


result<-matrix(0,2,1)
colnames(result)<-c("sf");rownames(result)<-c("lambda.hat","sf.pvalue")
result[,1]<-c(sflam,sf.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", sflam, "\n", sep = " ")
cat("  p.value    :", sf.pvalue, "\n\n", sep = " ")
cat(if(sf.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Shapiro-Francia"
out$lambda.hat<-as.numeric(sflam)
out$p.value<-sf.pvalue
out$tf.data<-data.sf
out$date<-date()
invisible(out)
}

else if (method=="lt") {

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

if (ltlam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (ltlam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (ltlam!=0) data.lt<-((data^ltlam)-1)/ltlam
if (ltlam==0) data.lt<-log(data)

lt.pvalue<-lillie.test(data.lt)$p.value

result<-matrix(0,2,1)
colnames(result)<-c("lt");rownames(result)<-c("lambda.hat","lt.pvalue")
result[,1]<-c(ltlam,lt.pvalue)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", ltlam, "\n", sep = " ")
cat("  p.value    :", lt.pvalue, "\n\n", sep = " ")
cat(if(lt.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Lilliefors"
out$lambda.hat<-as.numeric(ltlam)
out$p.value<-lt.pvalue
out$tf.data<-data.lt
out$date<-date()
invisible(out)
}

else if (method=="jb") {

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

if (jblam==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (jblam==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (jblam!=0) data.jb<-((data^jblam)-1)/jblam
if (jblam==0) data.jb<-log(data)

jb.pvalue<-jarque.bera.test(data.jb)$p.value

result<-matrix(0,2,1)
colnames(result)<-c("jb");rownames(result)<-c("lambda.hat","jb.pvalue")
result[,1]<-c(jblam,jb.pvalue)

if (verbose){

cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", jblam, "\n", sep = " ")
cat("  p.value    :", jb.pvalue, "\n\n", sep = " ")
cat(if(jb.pvalue > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}

out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Jarque-Bera"
out$lambda.hat<-as.numeric(jblam)
out$p.value<-as.numeric(jb.pvalue)
out$tf.data<-data.jb
out$date<-date()

invisible(out)
}

else if (method=="ac") {
aclam1<-NULL

set.seed(100)
for (t in 1:30) {

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

if (aclam2==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
if (aclam2==min(lam)) stop("Enlarge the range of the lambda in a negative direction")


if (aclam2!=0) data.ac<-((data^aclam2)-1)/aclam2
if (aclam2==0) data.ac<-log(data)


result<-matrix(0,1,1)
colnames(result)<-c("ac");rownames(result)<-c("lambda.hat")
result[,1]<-c(aclam2)

if (verbose){
cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("-----------------------------------------------", "\n", sep = " ")
cat("  data :", dname, "\n\n", sep = " ")
cat("  lambda.hat :", aclam2, "\n", sep = " ")
cat("  p.value    :", shapiro.test(data.ac)$p, "\n\n", sep = " ")
cat(if(shapiro.test(data.ac)$p > alpha){"  Result     : Transformed data are normal."}
     else {"  Result     : Transformed data are not normal."},"\n")
    cat("-----------------------------------------------", "\n\n", sep = " ")
}


out<-list()
out$title<-"Implementation of Box-Cox Power Transformation when No Covariate Is Available"
out$method="Artificial Covariate"
out$lambda.hat<-as.numeric(aclam2)
out$tf.data<-data.ac
out$date<-date()
invisible(out)
}

}
