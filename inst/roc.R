
library(plotROC)

getROC4 <- function(m0, m1, m2, m3, dat, endpoint) {
  
  dat <- m0$data
  p0 <- predict(m0, type="response")
  p1 <- predict(m1, type="response", newdata=m1$data)
  p2 <- predict(m2, type="response", newdata=m2$data)
  p3 <- predict(m3, type="response", newdata=m3$data)
  
  df <- rbind(data.frame(predictor = p0,
                         known.truth = as.numeric(as.character(dat[[endpoint]])),
                         Model="AvgRO% alone"),
              data.frame(predictor = p1,
                         known.truth = as.numeric(as.character(dat[[endpoint]])),
                         Model="AUC alone"),
              data.frame(predictor = p2,
                         known.truth = as.numeric(as.character(dat[[endpoint]])),
                         Model="Full"),
              data.frame(predictor = p3,
                         known.truth = as.numeric(as.character(dat[[endpoint]])),
                         Model="Reduced (p<0.05)"))
  
  df$Model <- ordered(df$Model, c("AvgRO% alone","AUC alone","Full","Reduced (p<0.05)"))
  roc <- ggplot(df, aes(d= known.truth, m = predictor, col=Model)) + 
    geom_roc(n.cuts=0) +
    style_roc(theme = theme_light) +
    scale_color_manual(values=pal_futurama()(4)) +
    theme(legend.position = "bottom")#+
  #labs(subtitle="Better AUC indicates better predictive performance. 0.9-1 = excellent;  0.8-0.9 = good; 0.7-0.8 = fair; 0.6-0.7 = poor; 0.5-0.6 = fail")
  
  auc <- calc_auc(roc)
  auc$Model <- c("AvgRO% alone","AUC alone","Full","Reduced (p<0.05)")
  auc$Model <- ordered(auc$Model, c("AvgRO% alone","AUC alone","Full","Reduced (p<0.05)"))
  auc$predictor <- c(0.75, 0.75, 0.75, 0.75)
  auc$known.truth <- c(0.35, 0.3, 0.25, 0.2)
  
  roc + geom_text(data=auc, aes(x=predictor, y=known.truth, label = paste("AUC[ROC] = ", fmt_signif(AUC,3)), col=Model)) + 
    plot_annotation(caption="Better AUC indicates better predictive performance.\n0.9-1 = excellent;  0.8-0.9 = good; 0.7-0.8 = fair; 0.6-0.7 = poor; 0.5-0.6 = fail")
} 
