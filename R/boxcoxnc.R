boxcoxnc <-
  function(data, method="sw", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, alpha = 0.05, verbose = TRUE)
  {
    dname<-deparse(substitute(data))
    
    data<-as.numeric(data)
   
    
    if(is.null(lambda2)) lambda2<-0
    
    data <- data+lambda2

    if (is.na(min(data))==TRUE) stop("Data include NA")
    if (min(data)<=0) stop("Data must include positive values. Specify shifting parameter, lambda2")



    if (method=="sw") {
  
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) shapiro.test(store2[[x]])$statistic)
     
pred.lamb<-store1[[which.max(store3)]]
method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Wilk test statistic"
      

    }
    
    else if (method=="ad") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) ad.test(store2[[x]])$statistic)
      
      pred.lamb<-store1[[which.min(store3)]]

      method.name<-"Estimating Box-Cox transformation parameter via Anderson-Darling test statistic"
    
      
    }
    
    else if (method=="cvm") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) cvm.test(store2[[x]])$statistic)
      
      pred.lamb<-store1[[which.min(store3)]]
      method.name<-"Estimating Box-Cox transformation parameter via Cramer-von Mises test statistic"
 
      
    }
    
    else if (method=="pt") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) pearson.test(store2[[x]])$statistic)
      
      pred.lamb<-store1[[which.min(store3)]]
      method.name<-"Estimating Box-Cox transformation parameter via Pearson Chi-Square test statistic"
      
    }
    
    else if (method=="sf") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) sf.test(store2[[x]])$statistic)
      
      pred.lamb<-store1[[which.max(store3)]]
      method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Francia test statistic"
      

    }
    
    else if (method=="lt") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) lillie.test(store2[[x]])$statistic)
      
      pred.lamb<-store1[[which.min(store3)]]
      method.name<-"Estimating Box-Cox transformation parameter via Lilliefors test statistic"

 
    }
    
    else if (method=="jb") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/store1[[x]] else log(data))
store3<-lapply(1:length(lambda), function(x) jarque.bera.test(store2[[x]])$statistic)
         
      pred.lamb<-store1[[which.min(store3)]]
      method.name<-"Estimating Box-Cox transformation parameter via Jarque-Bera test statistic"


    }
    
    else if (method=="ac") {
   
set.seed(100)
stor1<-lapply(1:30, function(x) rnorm(length(data),0,100))
stor2<-lapply(1:30, function(x) glm(data~stor1[[x]],family=gaussian))
stor3<-lapply(1:30, function(x) boxcox(stor2[[x]],lambda,plotit=FALSE))
stor4<-sapply(1:30, function(x) stor3[[x]]$x[which.max(stor3[[x]]$y)])

pred.lamb<-mean(stor4)

      method.name<-"Estimating Box-Cox transformation parameter via artificial covariate method"

    }

    else if (method=="mle") {
      
store1<-lapply(1:length(lambda), function(x) lambda[x])
store2<-lapply(1:length(lambda), function(x) if (store1[[x]] != 0) (data^store1[[x]]-1)/(store1[[x]]*(geometric.mean(data)^(store1[[x]]-1))) else geometric.mean(data)*log(data))
store3<-lapply(1:length(lambda), function(x) sum(log(dnorm(store2[[x]], mean = mean(store2[[x]]), sd = sd(store2[[x]])))))

pred.lamb<-store1[[which.max(store3)]]         
method.name<-"Estimating Box-Cox transformation parameter via maximum likelihood estimation"

    }





    if (pred.lamb==max(lambda)) stop("Enlarge the range of the lambda")
    if (pred.lamb==min(lambda)) stop("Enlarge the range of the lambda")
      


    if (pred.lamb!=0) data.transformed<-((data^pred.lamb)-1)/pred.lamb
    if (pred.lamb==0) data.transformed<-log(data)


    if(plot){
           

      par(mfrow=c(2,2))
      hist(data, xlab = dname, prob=TRUE, main = paste("Histogram of", dname))
      lines(density(data))
      hist(data.transformed, xlab = paste("Transformed", dname), prob=TRUE, main = paste("Histogram of tf", dname))
      lines(density(data.transformed))
      qqnorm(data, main = paste("Q-Q plot of", dname))
      qqline(data)
      qqnorm(data.transformed, main = paste("Q-Q plot of tf", dname))
      qqline(data.transformed)
           
           
    }




      if (method=="sw") {
         statistic<-shapiro.test(data.transformed)$statistic
         pvalue<-shapiro.test(data.transformed)$p.value
         nortest.name<-"Shapiro-Wilk normality test"

      }
      if (method=="ad") {
         statistic<-ad.test(data.transformed)$statistic
         pvalue<-ad.test(data.transformed)$p.value
         nortest.name<-"Anderson-Darling normality test"

      }
      if (method=="cvm") {
         statistic<-cvm.test(data.transformed)$statistic
         pvalue<-cvm.test(data.transformed)$p.value
         nortest.name<-"Cramer-von Mises normality test"

      }
      if (method=="pt") {
         statistic<-pearson.test(data.transformed)$statistic
         pvalue<-pearson.test(data.transformed)$p.value
         nortest.name<-"Pearson Chi-square normality test"   

      }
      if (method=="sf") {
         statistic<-sf.test(data.transformed)$statistic
         pvalue<-sf.test(data.transformed)$p.value
         nortest.name<-"Shapiro-Francia normality test"   

      }
      if (method=="lt") {
         statistic<-lillie.test(data.transformed)$statistic
         pvalue<-lillie.test(data.transformed)$p.value
         nortest.name<-"Lilliefors normality test"   

      }
      if (method=="jb") {
         statistic<-jarque.bera.test(data.transformed)$statistic
         pvalue<-jarque.bera.test(data.transformed)$p.value
         nortest.name<-"Jarque-Bera normality test"   

      }
      if ((method=="ac")|(method=="mle")) {
         statistic<-shapiro.test(data.transformed)$statistic
         pvalue<-shapiro.test(data.transformed)$p.value
         nortest.name<-"Shapiro-Wilk normality test"   

      }
      
      
      if (verbose){
        cat("\n"," Box-Cox power transformation", "\n", sep = " ")
        cat("-------------------------------------------------------------------", "\n\n", sep = " ")
        cat("  lambda.hat :", pred.lamb, "\n\n", sep = " ")
	cat("\n", "  ",nortest.name," for transformed data ", "(alpha = ",alpha,")", "\n", sep = "")
        cat("-------------------------------------------------------------------", "\n\n", sep = " ")
        cat("  statistic  :", statistic, "\n", sep = " ")
        cat("  p.value    :", pvalue, "\n\n", sep = " ")
        cat(if(pvalue > alpha){"  Result     : Transformed data are normal."}
            else {"  Result     : Transformed data are not normal."},"\n")
        cat("-------------------------------------------------------------------", "\n\n", sep = " ")
      }
      
      out<-list()
      out$method <- method.name
      out$lambda.hat <- as.numeric(pred.lamb)
      out$lambda2 <- as.numeric(lambda2)
      out$statistic <- as.numeric(statistic)
      out$p.value <- as.numeric(pvalue)
      out$alpha <- as.numeric(alpha)
      out$tf.data <- data.transformed
      out$var.name <- dname
      attr(out, "class") <- "boxcoxnc"
      invisible(out)


  }
