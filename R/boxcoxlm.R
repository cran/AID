boxcoxlm <-function(x, y, method="lse", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, alpha = 0.05, verbose = TRUE)
  {
    
    dname1<-deparse(substitute(y))
    dname2<-deparse(substitute(x))
    if(typeof(y)=="list"){
      dname1 <-colnames(as.data.frame(y))
      
    }
    if(typeof(x)=="list"){
      dname2 <-colnames(as.data.frame(x))
      
    }else{
      dname2 <-deparse(substitute(x))
    }
    x<-data.matrix(x)
    y<-as.numeric(data.matrix(y))
    
    
    if(is.null(lambda2)) lambda2<-0
    
    y <- y + lambda2
    
    if (!any(class(x)=="matrix")) stop("x must be a matrix")
    
    if (is.na(min(y))==TRUE) stop("x matrix includes NA")
    if (is.na(min(x))==TRUE) stop("response y includes NA")
    
    if (min(y)<=0) stop("response y must include positive values. Specify shifting parameter, lambda2")
    
    x<-cbind(1,x)
    
    
    
    
    if (method=="sw") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) shapiro.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.max(store6)]]
      
      method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Wilk test statistic"
      
      
    }
    
    else if (method=="ad") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) ad.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.min(store6)]]
      
      method.name<-"Estimating Box-Cox transformation parameter via Anderson-Darling test statistic"
      
      
    }
    
    else if (method=="cvm") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) cvm.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.min(store6)]]
      method.name<-"Estimating Box-Cox transformation parameter via Cramer-von Mises test statistic"
      
      
    }
    
    else if (method=="pt") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) pearson.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.min(store6)]]
      method.name<-"Estimating Box-Cox transformation parameter via Pearson Chi-Square test statistic"
      
    }
    
    else if (method=="sf") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) sf.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.max(store6)]]
      method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Francia test statistic"
      
      
    }
    
    else if (method=="lt") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) lillie.test(store5[[i]])$statistic)
      
      pred.lamb<-store1[[which.min(store6)]]
      method.name<-"Estimating Box-Cox transformation parameter via Lilliefors test statistic"
      
      
    }
    
    else if (method=="jb") {
      
      store1<-lapply(1:length(lambda), function(i) lambda[i])
      store2<-lapply(1:length(lambda), function(i) if (store1[[i]] != 0) (y^store1[[i]]-1)/store1[[i]] else log(y))
      store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
      store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
      store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
      
      store6<-lapply(1:length(lambda), function(i) jarque.bera.test(store5[[i]][1,])$statistic)
      
      pred.lamb<-store1[[which.min(store6)]]
      method.name<-"Estimating Box-Cox transformation parameter via Jarque-Bera test statistic"
      
      
    }
    
    else if (method=="mle") {
      
store1<-lapply(1:length(lambda), function(i) lambda[i])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (y^store1[[x]]-1)/(store1[[x]]*(geometric.mean(y)^(store1[[x]]-1))) else geometric.mean(y)*log(y))
store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
store6<-lapply(1:length(lambda), function(i) sum(log(dnorm(store5[[i]], mean = mean(store5[[i]]), sd = sd(store5[[i]])))))

pred.lamb<-store1[[which.max(store6)]]         
method.name<-"Estimating Box-Cox transformation parameter via maximum likelihood estimation"
      
      
    }
    
    else if (method=="lse") {
      
  store1<-lapply(1:length(lambda), function(i) lambda[i])
  store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (y^store1[[x]]-1)/(store1[[x]]*(geometric.mean(y)^(store1[[x]]-1))) else geometric.mean(y)*log(y))
  store3<-lapply(1:length(lambda), function(i) ginv(t(x) %*% x) %*% t(x) %*% store2[[i]])
  store4<-lapply(1:length(lambda), function(i) t(store3[[i]]) %*% t(x))
  store5<-lapply(1:length(lambda), function(i) store4[[i]]-store2[[i]])
  store6<-lapply(1:length(lambda), function(i) sum(store5[[i]]^2))
  
  pred.lamb<-store1[[which.min(store6)]]
  method.name<-"Estimating Box-Cox transformation parameter via least square estimation"
      
    }
    
  
    if (pred.lamb==max(lambda)) stop("Enlarge the range of the lambda")
    if (pred.lamb==min(lambda)) stop("Enlarge the range of the lambda")
    
    
    
    
    coef = ginv(t(x) %*% x) %*% t(x) %*% y
    ypred = t(coef) %*% t(x)
    residual = ypred - y
    
    if (pred.lamb!=0) y.transformed<-((y^pred.lamb)-1)/pred.lamb
    if (pred.lamb==0) y.transformed<-log(y)
    
    
    coef.transformed = ginv(t(x) %*% x) %*% t(x) %*% y.transformed
    ypred.transformed = t(coef.transformed) %*% t(x)
    residual.transformed = ypred.transformed - y.transformed
    
    
    if(plot){
      
      
      par(mfrow=c(2,2))
      hist(residual, xlab = "Residuals", prob=TRUE, main = "Histogram of residuals")
      lines(density(residual))
      hist(residual.transformed, xlab = "Residuals after transformation", prob=TRUE, main = paste("Histogram of residuals after transformation"))
      lines(density(residual.transformed))
      qqnorm(residual, main = "Q-Q plot of residuals")
      qqline(residual)
      qqnorm(residual.transformed, main = "Q-Q plot of residuals after transformation")
      qqline(residual.transformed)
      
      
    }
    
    
    
    
    if (method=="sw") {
      statistic<-shapiro.test(residual.transformed)$statistic
      pvalue<-shapiro.test(residual.transformed)$p.value
      nortest.name<-"Shapiro-Wilk normality test"
      
    }
    if (method=="ad") {
      statistic<-ad.test(residual.transformed)$statistic
      pvalue<-ad.test(residual.transformed)$p.value
      nortest.name<-"Anderson-Darling normality test"
      
    }
    if (method=="cvm") {
      statistic<-cvm.test(residual.transformed)$statistic
      pvalue<-cvm.test(residual.transformed)$p.value
      nortest.name<-"Cramer-von Mises normality test"
      
    }
    if (method=="pt") {
      statistic<-pearson.test(residual.transformed)$statistic
      pvalue<-pearson.test(residual.transformed)$p.value
      nortest.name<-"Pearson Chi-square normality test"   
      
    }
    if (method=="sf") {
      statistic<-sf.test(residual.transformed)$statistic
      pvalue<-sf.test(residual.transformed)$p.value
      nortest.name<-"Shapiro-Francia normality test"   
      
    }
    if (method=="lt") {
      statistic<-lillie.test(residual.transformed)$statistic
      pvalue<-lillie.test(residual.transformed)$p.value
      nortest.name<-"Lilliefors normality test"   
      
    }
    if (method=="jb") {
      statistic<-jarque.bera.test(residual.transformed[1,])$statistic
      pvalue<-jarque.bera.test(residual.transformed[1,])$p.value
      nortest.name<-"Jarque-Bera normality test"   
      
    }
    if ((method=="mle")|(method=="lse")) {
      statistic<-shapiro.test(residual.transformed)$statistic
      pvalue<-shapiro.test(residual.transformed)$p.value
      nortest.name<-"Shapiro-Wilk normality test"   
  
    }
    
    
    if (verbose){
      cat("\n"," Box-Cox power transformation", "\n", sep = " ")
      cat("-----------------------------------------------------------", "\n", sep = " ")
      cat("  data :", "Residuals", "\n\n", sep = " ")
      cat("  lambda.hat :", pred.lamb, "\n\n", sep = " ")
      cat("\n", "  ",nortest.name," (alpha = ",alpha,")", "\n", sep = "")
      cat("----------------------------------------------------", "\n\n", sep = " ")
      cat("  statistic  :", statistic, "\n", sep = " ")
      cat("  p.value    :", pvalue, "\n\n", sep = " ")
      cat(if(pvalue > alpha){"  Result     : Residuals are normal after transformation."}
          else {"  Result     : Residuals are not normal after transformation."},"\n")
      cat("-----------------------------------------------------------", "\n\n", sep = " ")
    }
    
    out<-list()
    out$method<-method.name
    out$lambda.hat<-as.numeric(pred.lamb)
    out$lambda2<-as.numeric(lambda2)
    out$statistic<-as.numeric(statistic)
    out$p.value<-as.numeric(pvalue)
    out$alpha<-as.numeric(alpha)
    out$tf.data<-as.numeric(y.transformed)
    out$tf.residuals<-as.numeric(residual.transformed[1,])
    out$y.name<-dname1
    out$x.name<-dname2
    attr(out, "class") <- "boxcoxlm"
    invisible(out)
    
    
}
