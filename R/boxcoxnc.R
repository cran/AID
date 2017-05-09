boxcoxnc <-
  function(data, method="sw", lam = seq(-3,3,0.01), plot = TRUE, alpha = 0.05, verbose = TRUE)
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
      
      pred.lamb<-sw[which.max(sw[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Wilk test statistic"
      nortest.name<-"Shapiro-Wilk normality test"
      

    }
    
    else if (method=="ad") {
      
      ad<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) ad<-rbind(ad,c(lam[i],ad.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) ad<-rbind(ad,c(lam[i],ad.test(log(data))$statistic))
      }
      
      pred.lamb<-ad[which.min(ad[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Anderson-Darling test statistic"
      nortest.name<-"Anderson-Darling normality test"
    
      
    }
    
    else if (method=="cvm") {
      
      # Cramer Von-Mises - in nortest
      cvm<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) cvm<-rbind(cvm,c(lam[i],cvm.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) cvm<-rbind(cvm,c(lam[i],cvm.test(log(data))$statistic))
      }
      
      pred.lamb<-cvm[which.min(cvm[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Cramer-von Mises test statistic"
      nortest.name<-"Cramer-von Mises normality test"
 
      
    }
    
    else if (method=="pt") {
      
      pt<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) pt<-rbind(pt,c(lam[i],pearson.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) pt<-rbind(pt,c(lam[i],pearson.test(log(data))$statistic))
      }
      
      pred.lamb<-pt[which.min(pt[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Pearson Chi-Square test statistic"
      nortest.name<-"Pearson Chi-square normality test"   
      
    }
    
    else if (method=="sf") {
      
      sf<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) sf<-rbind(sf,c(lam[i],sf.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) sf<-rbind(sf,c(lam[i],sf.test(log(data))$statistic))
      }
      
      pred.lamb<-sf[which.max(sf[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Shapiro-Francia test statistic"
      nortest.name<-"Shapiro-Francia normality test"   
      

    }
    
    else if (method=="lt") {
      
      lt<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) lt<-rbind(lt,c(lam[i],lillie.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) lt<-rbind(lt,c(lam[i],lillie.test(log(data))$statistic))
      }
      
      pred.lamb<-lt[which.min(lt[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Lilliefors test statistic"
      nortest.name<-"Lilliefors normality test"   

 
    }
    
    else if (method=="jb") {
      
      jb<-NULL
      for (i in 1:length(lam)){
        if (round(lam[i],2)!=0) jb<-rbind(jb,c(lam[i],jarque.bera.test((data**(lam[i])-1)/(lam[i]))$statistic))
        if (round(lam[i],2)==0) jb<-rbind(jb,c(lam[i],jarque.bera.test(log(data))$statistic))
      }
      
      
      pred.lamb<-jb[which.min(jb[,2]),1]
      method.name<-"Estimating Box-Cox transformation parameter via Jarque-Bera test statistic"
      nortest.name<-"Jarque-Bera normality test"   


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
      
   
      aclam1<-as.numeric(aclam1)
      pred.lamb<-mean(aclam1)
      method.name<-"Estimating Box-Cox transformation parameter via artificial covariate method"
      nortest.name<-"Shapiro-Wilk normality test"   


    }
    

    if (pred.lamb==max(lam)) stop("Enlarge the range of the lambda in a positive direction")
    if (pred.lamb==min(lam)) stop("Enlarge the range of the lambda in a negative direction")
      


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
      }
      if (method=="ad") {
         statistic<-ad.test(data.transformed)$statistic
         pvalue<-ad.test(data.transformed)$p.value
      }
      if (method=="cvm") {
         statistic<-cvm.test(data.transformed)$statistic
         pvalue<-cvm.test(data.transformed)$p.value
      }
      if (method=="pt") {
         statistic<-pearson.test(data.transformed)$statistic
         pvalue<-pearson.test(data.transformed)$p.value
      }
      if (method=="sf") {
         statistic<-sf.test(data.transformed)$statistic
         pvalue<-sf.test(data.transformed)$p.value
      }
      if (method=="lt") {
         statistic<-lillie.test(data.transformed)$statistic
         pvalue<-lillie.test(data.transformed)$p.value
      }
      if (method=="jb") {
         statistic<-jarque.bera.test(data.transformed)$statistic
         pvalue<-jarque.bera.test(data.transformed)$p.value
      }
      if (method=="ac") {
         statistic<-shapiro.test(data.transformed)$statistic
         pvalue<-shapiro.test(data.transformed)$p.value
      }
   

      
      if (verbose){
        cat("\n"," Box-Cox power transformation", "\n", sep = " ")
        cat("--------------------------------------------------------", "\n", sep = " ")
        cat("  data :", dname, "\n\n", sep = " ")
        cat("  lambda.hat :", pred.lamb, "\n\n", sep = " ")
	cat("\n", "",nortest.name,"for transformed data", "\n", sep = " ")
        cat("--------------------------------------------------------", "\n\n", sep = " ")
        cat("  statistic  :", statistic, "\n", sep = " ")
        cat("  p.value    :", pvalue, "\n\n", sep = " ")
        cat(if(pvalue > alpha){"  Result     : Transformed data are normal."}
            else {"  Result     : Transformed data are not normal."},"\n")
        cat("--------------------------------------------------------", "\n\n", sep = " ")
      }
      
      out<-list()
      out$method=method.name
      out$lambda.hat<-as.numeric(pred.lamb)
      out$statistic<-as.numeric(statistic)
      out$p.value<-as.numeric(pvalue)
      out$tf.data<-data.transformed
      invisible(out)


  }
