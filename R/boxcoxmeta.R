boxcoxmeta<-function(data, 
                     lambda = seq(-3, 3, 0.01), nboot = 100,
                     lambda2 = NULL, plot = TRUE, alpha = 0.05, verbose = TRUE){
  method = c("sw","ad","jb")
  if(is.null(lambda2)) lambda2<-0
  
  if(length(method)<2){
    boxcoxnc(data,method=method,alpha=alpha,lambda = lambda,lambda2 = lambda2)
  }else{
    
    lambdas<-data.frame(matrix(nrow=1, ncol=length(method)))
    boost_lambda<-data.frame(matrix(nrow = nboot, ncol = length(method)))
    colnames(boost_lambda)<-method
    colnames(lambdas)<-method
    
    for(i in method){
      lambdas[,i]<-boxcoxnc(data, method = i, plot =FALSE, verbose = FALSE,lambda = lambda,lambda2 = lambda2,alpha = alpha)$lambda.hat
      for (j in c(1:nboot)) {
        sample<-sample(1:length(data),length(data),replace = TRUE)
        boost_lambda[j,i]<-boxcoxnc(data[sample],method = i, plot = FALSE,verbose =FALSE,lambda = lambda,lambda2 = lambda2,alpha = alpha)$lambda.hat
      }
    }
    
    sd <- apply(boost_lambda,2,sd)
    pred.lamb<-metamean(n=rep(length(data),length(method)),mean =as.double(lambdas[method]),sd=as.double(sd[method]) )$TE.random
    
    if (pred.lamb == max(lambda)) stop("Enlarge the range of the lambda")
    if (pred.lamb == min(lambda)) stop("Enlarge the range of the lambda")
    if (pred.lamb != 0) data.transformed <- ((data^pred.lamb) - 1)/pred.lamb
    if (pred.lamb == 0) data.transformed <- log(data)
    
    if(typeof(data)=="list"){
      dname <-colnames(as.data.frame(data))
      if(length(dname)>1)dname<-"list"
      data<-unlist(data)
    }else{
      dname <-deparse(substitute(data))
      data <- as.numeric(data)
    }
    nortest.name <- str_replace_all(paste(method,collapse = " "),pattern = " ",replacement = ",")
    
    results<-data.frame(matrix(nrow=length(method),ncol=4))
    colnames(results)<-c("Test","Statistic","P.Value","Normality")
    row.names(results)<-method
    for (i in method) {
      if(i=="sw"){
        results[i,"Test"]<-"Shapiro-Wilk"
        results[i,"Statistic"]<-shapiro.test(data.transformed)$statistic
        results[i,"P.Value"]<-shapiro.test(data.transformed)$p.value
        results[i,"Normality"]<-ifelse(results[i,"P.Value"]<0.05,"Reject","Not reject")
      }else if(i=="ad"){
        results[i,"Test"]<-"Anderson Darling"
        results[i,"Statistic"]<-ad.test(data.transformed)$statistic
        results[i,"P.Value"]<-ad.test(data.transformed)$p.value
        results[i,"Normality"]<-ifelse(results[i,"P.Value"]<0.05,"Reject","Not reject")
      }else if(i=="jb"){
        results[i,"Test"]<-"Jarque-Bera"
        results[i,"Statistic"]<-jarque.bera.test(data.transformed)$statistic
        results[i,"P.Value"]<-jarque.bera.test(data.transformed)$p.value
        results[i,"Normality"]<-ifelse(results[i,"P.Value"]<0.05,"Reject","Not reject")
      }
    }
    row.names(results)<-NULL
    
    if (verbose){
      maxentry <- 60
      if (maxentry < 25) maxentry <- 25
      line<- paste(c(rep("-", round((maxentry + 
                                       10 - 32)/2, 0)), rep("-", 20), rep("-", 
                                                                          round((maxentry + 10 - 32)/2, 0))), sep = "")
      cat("\n"," Box-Cox power transformation via meta analysis", "\n", sep = " ")
      cat( line, sep = "")
      cat("\n","  data :", dname, "\n\n", sep = " ")
      cat("\n","  lambda.hat :", pred.lamb, "\n\n", sep = " ")
      cat("\n", "  ","Normality tests for transformed data ",
          "(alpha = ", alpha, ")", "\n",
          sep = "")
      cat(line[1:(length(line)-5)],"\n", sep = "")
      print(results)
      cat(line, sep = "")
    }
    if (plot) {
      par(mfrow = c(2, 2))
      hist(data, xlab = dname, prob = TRUE, main = paste("Histogram of", dname))
      lines(density(data))
      
      hist(data.transformed, xlab = paste("Transformed", dname), 
           prob = TRUE, main = paste("Histogram of tf", dname))
      lines(density(data.transformed))
      
      qqnorm(data, main = paste("Q-Q plot of", dname))
      qqline(data)
      
      qqnorm(data.transformed, main = paste("Q-Q plot of tf", dname))
      qqline(data.transformed)
    }
    out <- list()
    
    out$method <- "Ensemble Based Box-Cox Transformation via Meta Analysis"
    out$lambda.hat <- as.numeric(pred.lamb)
    out$lambda2 <- as.numeric(lambda2)
    out$result <- results
    out$alpha <- as.numeric(alpha)
    out$tf.data <- data.transformed
    out$var.name <- dname
    attr(out, "class") <- "boxcoxmeta"
    invisible(out)
    
  }
  
}
