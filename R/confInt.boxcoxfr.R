
confInt.boxcoxfr<- function(x, level = 0.95, plot = TRUE, xlab = NULL, ylab = NULL, title = NULL, width = NULL, verbose = TRUE,...){


if ((level<=0)|(level>=1)) stop("Confidence level must be between 0 and 1")
if (max((x$shapiro$p.value<x$alpha))==1) stop(paste("Transformed data in each group must be normally distributed at alpha = ",x$alpha,sep = ""))

if (max((x$shapiro$p.value<x$alpha))!=1){

k = length(levels(x$x))


stor1 = stor2 = stor3 = NULL
    for (i in 1:k) {
       datasub <- x$tf.data[which(x$x == (levels(x$x)[i]))]
       meantf<-mean(datasub)       
       lowertf <- meantf-qt((1-level)/2,df = (length(datasub)-1),lower.tail = FALSE)*sd(datasub)/sqrt(length(datasub))
       uppertf <- meantf+qt((1-level)/2,df = (length(datasub)-1),lower.tail = FALSE)*sd(datasub)/sqrt(length(datasub))


        stor1 = c(stor1, meantf)
        stor2 = c(stor2, lowertf)
        stor3 = c(stor3, uppertf)

    }

mattf <- cbind(stor1, stor2, stor3)


if (x$lambda.hat != 0) matbt <- (mattf*x$lambda.hat+1)^(1/x$lambda.hat)-x$lambda2
if (x$lambda.hat == 0) matbt <- exp(mattf)-x$lambda2

}

colnames(matbt)<-c("Mean", paste((1-level)/2*100, "%",sep = ""), paste((1-(1-level)/2)*100, "%",sep = ""))
rownames(matbt)<-levels(x$x)

     if (verbose){
        cat("\n"," Back transformed data", "\n", sep = " ")
        cat("-----------------------------------------", "\n", sep = " ")
        print(matbt)
        cat("-----------------------------------------", "\n\n", sep = " ")
      }

resp <- trt <- NULL
if (plot == TRUE){
df <- data.frame(trt = levels(x$x), resp = matbt[,1])
limits <- aes(ymax = matbt[,3], ymin = matbt[,2])
out <- ggplot(df, aes(y = resp, x = trt))

if (is.null(width)) width <- 0.15 else width <- width

out <- out + geom_point() + geom_errorbar(limits, width = width, size = 0.8)

if (is.null(ylab)) out <- out + ylab(x$y.name) else out <- out + ylab(ylab)

if (is.null(xlab)) out <- out + xlab(x$x.name) else out <- out + xlab(xlab)

if (is.null(title)) out <- out + ggtitle("") else out <- out + ggtitle(title)
plot(out)
}

invisible(matbt)

}
