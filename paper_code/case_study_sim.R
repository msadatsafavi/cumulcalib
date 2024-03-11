#March 10, 2024

#You need to run case_study.R first.

library(cumulcalib)
library(sqldf)
library(generalhoslem)

val_sample_sizes <- c(500,2000,8000,Inf)
n_sim <- 100

val_data <- data_us

val_data$pi1  <- predict(model1, type="response", newdata=val_data)
val_data$pi2 <- predict(model2, type="response", newdata=val_data)

out <- data.frame(i=integer(), sample_size=integer(), dev_type=character(), method=character(), pval=double())

index <- 0
for(i in 1:n_sim)
{
  cat('.')
  for(ss in val_sample_sizes)
  {
    if(is.infinite(ss)) ss<-nrow(val_data)
    this_data <- val_data[sample(nrow(val_data), size=ss, replace=T), ]

    res <- cumulcalib::cumulcalib(y=this_data$Y, p=this_data$pi1)
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "S"
    out[index,"method"] <- "BM"
    out[index,"pval"] <- res$by_method$BM$pval
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "S"
    out[index,"method"] <- "BB"
    out[index,"pval"] <- res$by_method$BB$pval

    res <- logitgof(this_data$Y, this_data$pi1)
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "S"
    out[index,"method"] <- "HL"
    out[index,"pval"] <- 1-pchisq(res$statistic,10)



    res <- cumulcalib::cumulcalib(y=this_data$Y, p=this_data$pi2)
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "L"
    out[index,"method"] <- "BM"
    out[index,"pval"] <- res$by_method$BM$pval
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "L"
    out[index,"method"] <- "BB"
    out[index,"pval"] <- res$by_method$BB$pval

    res <- logitgof(this_data$Y, this_data$pi2)
    index <- index+1
    out[index,"i"] <- i
    out[index,"sample_size"] <- ss
    out[index,"dev_type"] <- "L"
    out[index,"method"] <- "HL"
    out[index,"pval"] <- 1-pchisq(res$statistic,10)


  }
}


sqldf("SELECT sample_size, dev_type, method, round(AVG(pval),5) AS epval, round(AVG(pval<0.05),5) AS psig FROM out GROUP BY sample_size, dev_type, method")
