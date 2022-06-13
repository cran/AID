boxcoxfr <- function(y, x, option="both",lambda = seq(-3,3,0.01), lambda2 = NULL, tau = 0.05, alpha = 0.05, verbose = TRUE){

dname1<-deparse(substitute(y))
dname2<-deparse(substitute(x))


x=factor(x)

k=length(levels(x))


if (length(y) != length(x)) {stop("The lengths of x and y must be equal")}


if(is.null(lambda2)) lambda2<-0
    
y <- y+lambda2

if (is.na(min(y)) == TRUE) {stop("Data include NA")}

if (min(y) <= 0) {stop("Data must include positive values. Specify shifting parameter, lambda2")}


if(((option=="both")|(option=="nor")|(option=="var"))==F) {stop("Correct option argument")}


####################################
if ((option=="both")|(option=="nor")){

stor_w=NULL
for (i in 1:k){

for (j in 1:length(lambda)) {

if (lambda[j]!=0){
y1=y[which(x==(levels(x)[i]))]
w=(shapiro.test((y1^(lambda[j]) - 1)/(lambda[j])))
stor_w=rbind(stor_w,c(lambda[j],w$statistic,w$p))
}

if (lambda[j]==0){
y1=y[which(x==(levels(x)[i]))]
w=shapiro.test(log(y1))
stor_w=rbind(stor_w,c(lambda[j],w$statistic,w$p))
}

} #closing for loop

lambda=stor_w[which(stor_w[,3]>=tau),1]
if (length(lambda)==0) {stop("Feasible region is null set. No solution. \n  Try to enlarge the range of feasible lambda values, lambda. \n  Try to decrease feasible region parameter, tau.")}
stor_w=NULL

} #closing for loop
}
################################


##########

if ((option=="both")|(option=="var")){
stor_w=NULL
for (j in 1:length(lambda)) {

if (lambda[j]!=0){
lt=bartlett.test((y^(lambda[j]) - 1)/(lambda[j]),x)
stor_w=rbind(stor_w,c(lambda[j],lt$statistic,lt$p.value))
}

if (lambda[j]==0){
lt=bartlett.test(log(y),x)
stor_w=rbind(stor_w,c(lambda[j],lt$statistic,lt$p.value))
}
}
lambda=stor_w[which(stor_w[,3]>=tau),1]
if (length(lambda)==0) {stop("Feasible region is null set. No solution. \n  Try to enlarge the range of feasible lambda values, lambda. \n  Try to decrease feasible region parameter, tau.")}
}

##########



####

van=boxcox(y~x, lambda, plotit = FALSE)
lambda=van$x[which.max(van$y)]

####


################################


stor1=stor2=NULL
for(i in 1:k){

if(lambda!=0){

kk=shapiro.test((y[which(x==(levels(x)[i]))]^lambda-1)/lambda)

}else{
kk=shapiro.test(log(y[which(x==(levels(x)[i]))]))
}

stor1=c(stor1,kk$statistic)
stor2=c(stor2,kk$p)

} 

store = data.frame(matrix(NA, nrow = k, ncol = 4))
colnames(store) = c("Level", "statistic", "p.value", "Normality")
store$statistic=stor1
store$p.value=stor2
store$Normality = ifelse(store$p.value > alpha, "YES", "NO")
store$Level=levels(x)




if(lambda!=0){

kk2=bartlett.test((y^lambda-1)/lambda,x)

}else{

kk2=bartlett.test(log(y),x)
}

store2 = data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(store2) = c("Level","statistic", "p.value", "Homogeneity")
store2$statistic=kk2$statistic
store2$p.value=kk2$p.value
store2$Homogeneity= ifelse(store2$p.value > alpha, "YES", "NO")
store2$Level="All"


if(lambda!=0){
tf.data=(y^lambda-1)/lambda
}else{
tf.data=log(y)
}



if(tau==0){
method="MLE"
}else{
method="MLEFR"
}

if (verbose){

cat("\n"," Box-Cox power transformation", "\n", sep = " ")
cat("---------------------------------------------------------------------", "\n\n", sep = " ")
cat("  lambda.hat :", lambda, "\n\n", sep = " ")

cat("\n","  Shapiro-Wilk normality test for transformed data ","(alpha = ",alpha,")", "\n", sep = "")
cat("-------------------------------------------------------------------", "\n", sep = " ")
print(store)

cat("\n\n","  Bartlett's homogeneity test for transformed data ","(alpha = ",alpha,")", "\n", sep = "")
cat("-------------------------------------------------------------------", "\n", sep = " ")
print(store2)

cat("---------------------------------------------------------------------", "\n\n", sep = " ")
}


out <- list()
out$method <-method
out$lambda.hat <-lambda
out$lambda2 <-lambda2
out$shapiro <- store
out$bartlett <- store2
out$alpha<-as.numeric(alpha)
out$tf.data <- tf.data
out$x <- x
out$y.name <- dname1
out$x.name <- dname2

attr(out, "class") <- "boxcoxfr"
invisible(out)

}


 