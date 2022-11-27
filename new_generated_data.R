ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("mvtnorm","ggplot2","caret","ipred","devtools","sm","hrbrthemes","tidyr","dplyr","nnet","viridis","scatterplot3d", "patchwork", "plot3D")
ipak(packages)

library(MASS)
library(mvtnorm)
library(ggplot2)
library(caret)
library(ipred)
library(devtools)
library(sm)
library(hrbrthemes)
library(tidyr)
library(dplyr)
library(viridis)
library(scatterplot3d)
library(patchwork)
library(plotly)
library(e1071)
library(ggpubr)
library(rgl)
library(plot3D)

set.seed(2000)

d <- 2 # dimensions
N <- 2000 # number of generated samples
num_of_trials <- 100
num_of_classes <- 4
alphas <- seq(0.01,0.2,0.01)

labels.generated <- NULL

for(l in 1:N)
  labels.generated[l] <- sample(1:num_of_classes,prob = rep(1/num_of_classes,num_of_classes))
df <- matrix(0,ncol = d, nrow = N)

sigma <- matrix(c(1,0.1,0.1,1),ncol = d)
for (index in 1:N) {
  
  if(labels.generated[index] == 1)
     df[index,] <- mvrnorm(1,c(-1,0),sigma)
  else if(labels.generated[index] == 2)
     df[index,] <- mvrnorm(1,c(-1,-1),sigma)
  else if(labels.generated[index] == 3)
     df[index,] <- mvrnorm(1,c(0,-1),sigma)
  else
     df[index,] <- mvrnorm(1,c(1,-1),sigma)
}

df <- data.frame(cbind(df,labels.generated))
names(df)[d+1] <- "label"
View(df)

#scatter3D(df$V1,df$V2,df$V3,colvar=df$label,box=TRUE,pch=1,cex=1,bty="u",axes=TRUE,label=TRUE, nticks=2, ticktype="detailed",theta=40, phi=40, xlab="X-val", ylab="Y-val", zlab="Z-val", main="3D scatter plot")

# scatter3D(df$V1, df$V2, df$V3, phi = 0, colvar = df$label, bty = "g",pch = 20, cex = 1)


ggplot(df,aes(x=V1,y=V2)) + xlab(bquote(x[1])) + ylab(bquote(x[2])) + geom_point(aes(color=label)) + scale_color_viridis(option = "D")

marginal_coverage_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

set_size_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

ordinal_CDF_score_function <- function(calibration.scores,calibration.labels)
{
  cusum_cal <- list()
  cusum_cal <- vector("list",length = nrow(calibration.scores))
  argmaxes <- NULL
  conformity_scores <- NULL
  for (i in 1:nrow(calibration.scores)) {
    
    cusum_cal[[i]] <- cumsum(calibration.scores[i,])
    argmaxes[i] <- which.max(calibration.scores[i,])
    t <- cusum_cal[[i]]
    conformity_scores[i] <- t[argmaxes[i]] - t[calibration.labels[i]]
  }
  return(conformity_scores)
}

ordinal_CDF_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  cusum_val <- list()
  cusum_val <- vector("list",length = nrow(validation.scores))
  argmaxes <- NULL
  prediction_interval <- matrix(0,ncol = 2,nrow = nrow(validation.scores))
  for (j in 1:nrow(validation.scores)) 
  {
    cusum_val[[j]] <- cumsum(validation.scores[j,])
    argmaxes[j] <- which.max(validation.scores[j,])
    t <- cusum_val[[j]]
    prediction_interval[j,1] <- min(which(t >= t[argmaxes[j]]-qhat)) 
    prediction_interval[j,2] <- max(which(t <= t[argmaxes[j]]+qhat))
  }
  return(prediction_interval)
}

lac_score_function <- function(calibration.scores,calibration.labels)
{
  conformity_scores <- NULL
  for(i in 1:nrow(calibration.scores))
    conformity_scores[i] <- 1 - calibration.scores[i,calibration.labels[i]]
  return(conformity_scores)
}

naive_lac_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  lac.prediction.set <- matrix(0,nrow = nrow(validation.scores), ncol = num_of_classes)
  for(j in 1:nrow(validation.scores))
  {
    for(k in 1:num_of_classes)
    {
      if(validation.scores[j,k] >= 1 - qhat)
        lac.prediction.set[j,k] <- k
      else
        lac.prediction.set[j,k] <- NA
    }
  }
  return(lac.prediction.set)
}

opi_prediction_intervals <- function(p_values,level.of.significance)
{
  lower_edge <- NULL
  upper_edge <- NULL
  
  for(i in 1:nrow(p_values))
  {
    lower_edge[i] <- min(which(p_values[i,] >= level.of.significance))
    upper_edge[i] <- max(which(p_values[i,] >= level.of.significance))
  }
  return(cbind(lower_edge,upper_edge))
}


opi_prediction_sets <- function(p_values,level.of.significance)
{
  opi.prediction.set = list()
  opi.prediction.set = vector("list",length = nrow(p_values))
  for(i in 1:nrow(p.values))
  {
    idx <- which(p.values[i,] >= level.of.significance)
    opi.prediction.set[[i]] <- idx
  }
  return(opi.prediction.set)
}

evaluate_set_size_for_lac <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:nrow(prediction_set))
    set_size[i] <- sum(!is.na(prediction_set[i,]))
  return(mean(set_size))
}

evaluate_set_size_ops <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:length(prediction_set))
    set_size[i] <- length(as.vector(prediction_set[[i]]))
  return(mean(set_size))
}

estimated_qhat <- function(conformity_scores,a)
{
  q <- floor((cal_pct*(N-500) + 1)*(1 - a))/(cal_pct*(N-500) + 1)
  qhat <- quantile(conformity_scores, probs = c(0.25,q))[2]
  return(qhat)
}

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials)
  {
    training.index <- sample(1:N,500)
    df_training <- df[training.index,]
    df_test <- df[-training.index,]
    
    # logistic regression
    
    multinomial_model <- multinom(label~.,data = df_training)
    
    cal_pct <- 0.35
    
    calibration.indexes <- createDataPartition(y=df_test$label,p=cal_pct,list=FALSE)
    df_calibration <- df[calibration.indexes,]
    df_validation <- df[-calibration.indexes,]
    
    # Calibration data
    
    cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(d+1)],"probs")
    cal_labels <- df_calibration$label
    
    cal_scores <- NULL
    
    for(i in 1:nrow(cal_probs))
      cal_scores[i] <- cal_probs[i,cal_labels[i]]
    
    # Validation data
    
    val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
    val_labels <- df_validation$label
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] >= cal_scores)+1)/(length(cal_scores) + 1)
      }
    }
    
    ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
    lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
    
    Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
    Q_hat_lac <- estimated_qhat(lac_conformity_scores, alphas[alpha])
    
    # ordinal cdf 
    
    cdf_prediction_interval <- ordinal_CDF_prediction_sets(val_probs,val_labels,Q_hat_cdf)
    cdf_coverage <- 0
    
    for(m in 1:nrow(cdf_prediction_interval))
    {
      if(cdf_prediction_interval[m,1] <= val_labels[m] && val_labels[m] <= cdf_prediction_interval[m,2] )
        cdf_coverage <- cdf_coverage + 1
    }
    marginal_coverage_cdf[alpha,trial] <- cdf_coverage/length(val_labels)
    set_size_cdf[alpha,trial] <- mean(cdf_prediction_interval[,2] - cdf_prediction_interval[,1] + 1)
    
    # naive lac
    
    lac_prediction_set <- naive_lac_prediction_sets(val_probs,val_df,Q_hat_lac)
    lac_coverage <- 0
    
    for(m in 1:nrow(lac_prediction_set))
    {
      if(val_labels[m] %in% lac_prediction_set[m,])
        lac_coverage <- lac_coverage + 1
    }
    
    marginal_coverage_lac[alpha,trial] <- lac_coverage/length(val_labels)
    set_size_lac[alpha,trial] <- evaluate_set_size_for_lac(lac_prediction_set)
    
    # opi
    
    ordinal_prediction_intervals <- opi_prediction_intervals(p.values,alphas[alpha])
    opi_coverage <- 0
    
    for(m in 1:nrow(ordinal_prediction_intervals))
    {
      if(ordinal_prediction_intervals[m,1] <= val_labels[m] && val_labels[m] <= ordinal_prediction_intervals[m,2] )
        opi_coverage <- opi_coverage + 1
    }
    marginal_coverage_opi[alpha,trial] <- opi_coverage/length(val_labels)  
    set_size_opi[alpha,trial] <- mean(ordinal_prediction_intervals[,2] - ordinal_prediction_intervals[,1] + 1)
    
    # ops
    
    ordinal_prediction_sets <- opi_prediction_sets(p.values,alphas[alpha])
    ops_coverage <- 0
    
    for(m in 1:length(ordinal_prediction_sets))
    {
      if(val_labels[m] %in% as.vector(ordinal_prediction_sets[[m]]))
        ops_coverage <- ops_coverage + 1
    }
    marginal_coverage_ops[alpha,trial] <- ops_coverage/length(val_labels)
    set_size_ops[alpha,trial] <- evaluate_set_size_ops(ordinal_prediction_sets)
  }
}

cal_df <- cbind(cal_probs,cal_labels)
val_df <- cbind(val_probs,val_labels)
mydat = data.frame(rbind(cal_df,val_df))
write.csv(data.frame(cbind(cal_probs,cal_labels)),"/Users/Owner/Desktop/mydata.csv")

average_cvg_cdf <- NULL
average_cvg_lac <- NULL
average_cvg_opi <- NULL
average_cvg_ops <- NULL

avg_set_size_cdf <- NULL
avg_set_size_lac <- NULL
avg_set_size_opi <- NULL
avg_set_size_ops <- NULL

for(alpha in 1:length(alphas))
{
  average_cvg_cdf[alpha] <- mean(marginal_coverage_cdf[alpha,])
  average_cvg_lac[alpha] <- mean(marginal_coverage_lac[alpha,])
  average_cvg_opi[alpha] <- mean(marginal_coverage_opi[alpha,])
  average_cvg_ops[alpha] <- mean(marginal_coverage_ops[alpha,])
  
  avg_set_size_cdf[alpha] <- mean(set_size_cdf[alpha,][is.finite(set_size_cdf[alpha,])])
  avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,][is.finite(set_size_opi[alpha,])])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

average_cvg_aps <- c(1.0,
                     0.987719298245614,
                     0.976315789473684,
                     0.9681286549707602,
                     0.9523391812865498,
                     0.9444444444444444,
                     0.9374269005847955,
                     0.9345029239766083,
                     0.9157894736842106,
                     0.9076023391812866,
                     0.9029239766081872,
                     0.8812865497076023,
                     0.881578947368421,
                     0.87046783625731,
                     0.8616959064327485,
                     0.8421052631578947,
                     0.8494152046783625,
                     0.8128654970760234,
                     0.8,
                     0.7970760233918128)
avg_set_size_aps <- c(3.9239766081871337,
                      3.3432748538011694,
                      3.2374269005847958,
                      3.127192982456141,
                      2.9982456140350875,
                      2.864327485380117,
                      2.7918128654970764,
                      2.7491228070175437,
                      2.639473684210526,
                      2.588011695906433,
                      2.526608187134503,
                      2.3812865497076023,
                      2.369883040935673,
                      2.3020467836257312,
                      2.2567251461988307,
                      2.1953216374269005,
                      2.18859649122807,
                      2.083040935672515,
                      2.0105263157894737,
                      2.0125730994152047)

cond_coverage_lac <- matrix(0,nrow = num_of_classes,ncol = length(alphas))
cond_coverage_cdf <- matrix(0,nrow = num_of_classes,ncol = length(alphas))
cond_coverage_opi <- matrix(0,nrow = num_of_classes,ncol = length(alphas))
cond_coverage_ops <- matrix(0,nrow = num_of_classes,ncol = length(alphas))

for(alpha in 1:length(alphas))
{
  training.index <- sample(1:N,500)
  df_training <- df[training.index,]
  df_test <- df[-training.index,]
  
  multinomial_model <- multinom(label~.,data = df_training)
  
  cal_pct <- 0.35
  
  calibration.indexes <- createDataPartition(y=df_test$label,p=cal_pct,list=FALSE)
  df_calibration <- df[calibration.indexes,]
  df_validation <- df[-calibration.indexes,]
  
  # Calibration data
  
  cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(d+1)],"probs")
  cal_labels <- df_calibration$label
  
  cal_scores <- NULL
  
  for(i in 1:nrow(cal_probs))
    cal_scores[i] <- cal_probs[i,cal_labels[i]]
  
  # Validation data
  
  val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
  val_labels <- df_validation$label
  
  lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
  Q_hat_lac <- estimated_qhat(lac_conformity_scores, alphas[alpha])
  
  lac_prediction_set <- naive_lac_prediction_sets(val_probs,val_df,Q_hat_lac)
  
  for(class in 1:num_of_classes)
  {
    lac_conditional_cvg <- 0
    for(m in 1:nrow(lac_prediction_set))
    {
      if(val_labels[m] == class)
      {
        if(val_labels[m] %in% lac_prediction_set[m,])
          lac_conditional_cvg <- lac_conditional_cvg + 1
      }
    }
    cond_coverage_lac[class,alpha] <- lac_conditional_cvg/length(which(val_labels == class))
  }
  
  ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
  Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
  
  cdf_prediction_interval <- ordinal_CDF_prediction_sets(val_probs,val_labels,Q_hat_cdf)
  
  for(class in 1:num_of_classes)
  {
    cdf_conditional_cvg <- 0
    for(m in 1:nrow(cdf_prediction_interval))
    {
      if(val_labels[m] == class)
      {
        if(val_labels[m] %in% cdf_prediction_interval[m,])
          cdf_conditional_cvg <- cdf_conditional_cvg + 1
      }
    }
    cond_coverage_cdf[class,alpha] <- cdf_conditional_cvg/length(which(val_labels == class))
  }
  
  p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
  
  for(i in 1:nrow(val_probs))
  {
    for(j in 1:num_of_classes)
    {
      p.values[i,j] <- (sum(val_probs[i,j] >= cal_probs[which(cal_labels == j),j])+1)/(length(cal_probs[which(cal_labels == j),j]) + 1)
    }
  }
  
  ordinal_prediction_intervals <- opi_prediction_intervals(p.values,alphas[alpha])
  ordinal_prediction_sets <- opi_prediction_sets(p.values,alphas[alpha])
  
  for(class in 1:num_of_classes)
  {
    opi_cond_cvg <- 0
    ops_cond_cvg <- 0
    
    for(m in 1:nrow(ordinal_prediction_intervals))
    {
      if(val_labels[m] == class)
      {
        if(ordinal_prediction_intervals[m,1] <= val_labels[m] && val_labels[m] <= ordinal_prediction_intervals[m,2] )
          opi_cond_cvg <- opi_cond_cvg + 1
      }
    }
    for(m in 1:length(ordinal_prediction_sets))
    {
      if(val_labels[m] == class)
      {
        if(val_labels[m] %in% as.vector(ordinal_prediction_sets[[m]]))
          ops_cond_cvg <- ops_cond_cvg + 1
      }
    }
    
    cond_coverage_opi[class,alpha] <- opi_cond_cvg/length(which(val_labels == class))
    cond_coverage_ops[class,alpha] <- ops_cond_cvg/length(which(val_labels == class))
  }
}

aps_neg_cvg <- c(1.0,
                 0.9805505333493139,
                 0.9495251583574527,
                 0.9388899799869244,
                 0.9367870707609409,
                 0.9335108656023291,
                 0.9348922829944142,
                 0.9379334316650088,
                 0.9296380145808747,
                 0.9110370484717505,
                 0.9308883068965679,
                 0.9140365707346234,
                 0.9203211539854671,
                 0.9122764645680933,
                 0.901635438651503,
                 0.8957076141107511,
                 0.9165503246387704,
                 0.8765088530711644,
                 0.8560413255723212,
                 0.861330256414033)
aps_mil_cvg <- c(1.0,
                 0.9977777777777778,
                 0.9922834300675951,
                 0.9857812281002307,
                 0.9695734492565642,
                 0.9616129338820218,
                 0.9480210844355431,
                 0.9444387761118905,
                 0.9191389863690196,
                 0.9123146779640157,
                 0.9007283702315025,
                 0.8804195658089953,
                 0.8767139789234644,
                 0.8576704926015012,
                 0.8426680960601456,
                 0.837793709005747,
                 0.838173203015578,
                 0.79948266400602,
                 0.787175036565018,
                 0.7784801412240736)
aps_mod_cvg <- c(1.0,
                 1.0,
                 0.9913262977435151,
                 0.9843339131345544,
                 0.9645028754064752,
                 0.9566351883069922,
                 0.9510315570899515,
                 0.9562807051649559,
                 0.9272328417490465,
                 0.915869906096081,
                 0.892903549286223,
                 0.8792729905628262,
                 0.8665285474029186,
                 0.8552761966392299,
                 0.8512320516948393,
                 0.8122931485348823,
                 0.8229365457990392,
                 0.7814314482440757,
                 0.7497354660822287,
                 0.7508359077450943)
aps_sev_cvg <- c(1.0,
                 0.9728428009868619,
                 0.9693657895995452,
                 0.9598705502408139,
                 0.938323980735233,
                 0.9256528833530016,
                 0.9149927433890317,
                 0.9023686682970846,
                 0.890577697411613,
                 0.8923295745223369,
                 0.8899544673636924,
                 0.856063777760936,
                 0.8669428098064452,
                 0.8585057932903105,
                 0.8555595500403527,
                 0.8281662495073301,
                 0.827624000311779,
                 0.8004083236293997,
                 0.8091205024410899,
                 0.801044363688006)

cond_coverage_aps <- rbind(aps_neg_cvg,aps_mil_cvg,aps_mod_cvg,aps_sev_cvg)

lac_conditional_coverage <- NULL
cdf_conditional_coverage <- NULL
opi_conditional_coverage <- NULL
ops_conditional_coverage <- NULL
aps_conditional_coverage <- NULL

for(alpha in 1:length(alphas))
{
  lac_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_lac[,alpha],0)
  cdf_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_cdf[,alpha],0)
  opi_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_opi[,alpha],0)
  ops_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_ops[,alpha],0)
  aps_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_aps[,alpha],0)
}

conditional_coverage_gaussian <- cbind(cdf_conditional_coverage,lac_conditional_coverage,
                                       aps_conditional_coverage,opi_conditional_coverage,
                                       ops_conditional_coverage)

  
results_average_cvg_gaussian <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size_gaussian <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg_gaussian, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")



ggplt_2 <- ggplot(data = results_average_set_size_gaussian,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) + ylim(0,4) +
  labs(x = "alphas", y = "Set size", color = "Legend")

plt_multinom <- ggplt_1 / ggplt_2 

ggarrange(ggplt_1,ggplt_2,ncol = 2, nrow = 1,legend = "bottom")
