
confInt<- function(x,...) UseMethod("confInt")

confInt.default <- function(x,...) confInt.boxcoxnc(x,...)

confInt.boxcoxnc<- function(x, level = 0.95, verbose = TRUE,...){


if ((level<=0)|(level>=1)) stop("Confidence level must be between 0 and 1")
if (x$p.value<=x$alpha) stop(paste("Transformed data must be normally distributed at alpha = ",x$alpha,sep = ""))



if (x$p.value>x$alpha){


meantf <- mean(x$tf.data)
lowertf <- mean(x$tf.data)-qt((1-level)/2,df = (length(x$tf.data)-1),lower.tail = FALSE)*sd(x$tf.data)/sqrt(length(x$tf.data))
uppertf <- mean(x$tf.data)+qt((1-level)/2,df = (length(x$tf.data)-1),lower.tail = FALSE)*sd(x$tf.data)/sqrt(length(x$tf.data))
vectf <- c(meantf, lowertf, uppertf)
if (x$lambda.hat != 0) vecbt <- (vectf*x$lambda.hat+1)^(1/x$lambda.hat)
if (x$lambda.hat == 0) vecbt <- exp(vectf)

}

vecbt<- vecbt-x$lambda2
vecbt<- matrix(vecbt,1,3)
colnames(vecbt)<-c("Mean", paste((1-level)/2*100, "%",sep = ""), paste((1-(1-level)/2)*100, "%",sep = ""))
rownames(vecbt)<-x$var.name

     if (verbose){
        cat("\n"," Back transformed data", "\n", sep = " ")
        cat("---------------------------------------------", "\n", sep = " ")
        print(vecbt)
        cat("---------------------------------------------", "\n\n", sep = " ")
      }

invisible(vecbt)
}




