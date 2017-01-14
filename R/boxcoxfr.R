boxcoxfr <- function(y, x, option="both",lam = seq(-3,3,0.01), tau = 0.05, alpha = 0.05, verbose = TRUE){

dname1<-deparse(substitute(y))
dname2<-deparse(substitute(x))


x=factor(x)

k=length(levels(x))


if (length(y) != length(x)) {stop("The lengths of x and y must be equal")}

if (is.na(min(y)) == TRUE) {stop("Data include NA")}

if (min(y) <= 0) {stop("Data must include positive values")}

if(((option=="both")|(option=="nor")|(option=="var"))==F) {stop("Correct option argument")}


####################################
if ((option=="both")|(option=="nor")){

stor_w=NULL
for (i in 1:k){

for (j in 1:length(lam)) {

if (lam[j]!=0){
y1=y[which(x==(levels(x)[i]))]
w=(shapiro.test((y1^(lam[j]) - 1)/(lam[j])))
stor_w=rbind(stor_w,c(lam[j],w$statistic,w$p))
}

if (lam[j]==0){
y1=y[which(x==(levels(x)[i]))]
w=shapiro.test(log(y1))
stor_w=rbind(stor_w,c(lam[j],w$statistic,w$p))
}

} #closing for loop

lam=stor_w[which(stor_w[,3]>=tau),1]
if (length(lam)==0) {stop("Feasible region is null set. No solution. \n  Try to enlarge the range of feasible lambda values, lam. \n  Try to decrease feasible region parameter, tau.")}
stor_w=NULL

} #closing for loop
}
################################


##########

if ((option=="both")|(option=="var")){
stor_w=NULL
for (j in 1:length(lam)) {

if (lam[j]!=0){
lt=bartlett.test((y^(lam[j]) - 1)/(lam[j]),x)
stor_w=rbind(stor_w,c(lam[j],lt$statistic,lt$p.value))
}

if (lam[j]==0){
lt=bartlett.test(log(y),x)
stor_w=rbind(stor_w,c(lam[j],lt$statistic,lt$p.value))
}
}
lam=stor_w[which(stor_w[,3]>=tau),1]
if (length(lam)==0) {stop("Feasible region is null set. No solution.")}
}

##########



####

van=boxcox(y~x, lam, plotit = FALSE)
lam=van$x[which.max(van$y)]

####


################################


stor=NULL
for(i in 1:k){

if(lam!=0){

kk=shapiro.test((y[which(x==(levels(x)[i]))]^lam-1)/lam)

}else{
kk=shapiro.test(log(y[which(x==(levels(x)[i]))]))
}

stor=c(stor,kk$p)

} 

store = data.frame(matrix(NA, nrow = k, ncol = 3))
colnames(store) = c("Level", "p.value", "Normality")

store$p.value=stor
store$Normality = ifelse(store$p.value > alpha, "YES", "NO")
store$Level=levels(x)




if(lam!=0){

kk2=bartlett.test((y^lam-1)/lam,x)

}else{

kk2=bartlett.test(log(y),x)
}

store2 = data.frame(matrix(NA, nrow = 1, ncol = 3))
colnames(store2) = c("Level","p.value", "Homogenity")

store2$p.value=kk2$p.value
store2$Homogenity= ifelse(store2$p.value > alpha, "YES", "NO")
store2$Level="All"


if(lam!=0){
tf.data=(y^lam-1)/lam
}else{
tf.data=log(y)
}



if(tau==0){
method="MLE"
}else{
method="MLEFR"
}

if (verbose){

cat("\n"," Box-Cox Power Transformation", "\n", sep = " ")
cat("------------------------------------------------------------", "\n", sep = " ")
cat("  data :", dname1, "vs",dname2, "\n\n", sep = " ")
cat("  lambda.hat :", lam, "\n\n", sep = " ")

cat("\n"," Shapiro-Wilk Normality Test for Transformed Data", "\n", sep = " ")
cat("----------------------------------------------------", "\n", sep = " ")
print(store)

cat("\n\n"," Bartlett's Homogenity Test for Transformed Data", "\n", sep = " ")
cat("----------------------------------------------------", "\n", sep = " ")
print(store2)

cat("------------------------------------------------------------", "\n\n", sep = " ")
}


out <- list()
out$method <-method
out$lambda.hat <-lam
out$shapiro <- store
out$bartlett <- store2
out$tf.data <- tf.data
out$x <- x
out$date <- date()

invisible(out)

}


 