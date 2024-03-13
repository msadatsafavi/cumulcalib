library(rms)
library(pROC)
library(predtools)
library(cumulcalib)
library(predtools) #For  data
library(generalhoslem)
library(ggplot2)

set.seed(111111)

settings <- list(
  dev_sample_size1 = 500,
  dev_sample_size2 = Inf,
  val_sample_size1 = 500,
  val_sample_size2 = Inf
)


data(gusto)
gusto$kill <- (as.numeric(gusto$Killip)>1)*1
gusto$Y <- gusto$day30
data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]
data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15),]

if(is.infinite(settings$dev_sample_size1))
{
  dev_data1<- data_other
}else
{
  dev_data1 <- data_other[sample(1:(dim(data_other)[1]),settings$dev_sample_size1,F),]
}
if(is.infinite(settings$dev_sample_size2))
{
  dev_data2<- data_other
}else
{
  dev_data2 <- data_other[sample(1:(dim(data_other)[1]),settings$dev_sample_size2,F),]
}


if(is.infinite(settings$val_sample_size1))
{
  val_data1 <- data_us
}else
{
  val_data1 <- data_us[sample(1:(dim(data_us)[1]),settings$val_sample_size1,F),]  #data is for external validation
}
if(is.infinite(settings$val_sample_size2))
{
  val_data2 <- data_us
}else
{
  val_data2 <- data_us[sample(1:(dim(data_us)[1]),settings$val_sample_size2,F),]  #data is for external validation
}


#model <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + lsp(pulse,50), data=dev_data, family=binomial(link="logit"))
model1 <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data1, family=binomial(link="logit"))
model2 <- glm(Y ~ age + miloc + pmi + kill + pmin(sysbp,100) + pulse, data=dev_data2, family=binomial(link="logit"))

val_data1$pi1  <- predict(model1, type="response", newdata=val_data1)
val_data1$pi2 <- predict(model2, type="response", newdata=val_data1)

val_data2$pi1  <- predict(model1, type="response", newdata=val_data2)
val_data2$pi2 <- predict(model2, type="response", newdata=val_data2)

# coefficients(model1)
# pROC::auc(val_data$Y,val_data$pi1)
# coefficients(model2)
# pROC::auc(val_data$Y,val_data$pi2)

#2024.03.12
plot_loess <- function(Y, pi)
{
  legendPosition <- "right"
  lc <-loess(Y~pi)
  res <- predict(lc, se=T)
  df <- data.frame(x=pi, y=res$fit, yl=res$fit-1.96*res$se.fit, yh=res$fit+1.96*res$se.fit)

  plt <- ggplot(data=df, aes(x=x, y=y)) +
          geom_line(color="#619cff") +
          #geom_errorbar(aes(ymax = yh, ymin = yl)) +
          geom_line(aes(y=yl), color="#619cff", linetype=2) +
          geom_line(aes(y=yh), color="#619cff", linetype=2) +
          geom_line(aes(y=x), color="black", linetype=1) +
          labs(x="Prediction",  y="Observation") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   axis.text = element_text(colour = "black", size = 12),
                                   legend.position = legendPosition)

  plt
}

res11 <- res12 <- res21 <- res22 <-list()

res12$coeffs <- coefficients(model1)
res22$coeffs <- coefficients(model2)

res12$auc <- roc(val_data2$Y, val_data2$pi1)$auc
res22$auc <- roc(val_data2$Y, val_data2$pi2)$auc

res12$calib <- predtools::calibration_plot(val_data2, obs="Y", pred="pi1")
res12$scalib <- plot_loess(val_data2$Y, val_data2$pi1)
res12$cumul_calib <- cumulcalib::cumulcalib(val_data2$Y, val_data2$pi1)
res12$HL <- generalhoslem::logitgof(val_data2$Y, val_data2$pi1)

res22$calib <- predtools::calibration_plot(val_data2, obs="Y", pred="pi2")
res22$scalib <- plot_loess(val_data2$Y, val_data2$pi2)
res22$cumul_calib <- cumulcalib::cumulcalib(val_data2$Y, val_data2$pi2)
res22$HL <- generalhoslem::logitgof(val_data2$Y, val_data2$pi2)


saveRDS(list(res_small=res12,res_full=res22),"case_study.RDS")






