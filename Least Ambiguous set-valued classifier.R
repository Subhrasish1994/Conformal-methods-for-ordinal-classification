# Naive LAC method described in A.Angelopoulos, S.Bates (2020)

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("mvtnorm","ggplot2","caret","ipred","devtools","sm","hrbrthemes","tidyr","dplyr","nnet","viridis","plotly","scatterplot3d", "patchwork")
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
library(plotly)
library(viridis)
library(scatterplot3d)
library(patchwork)
library(e1071)

lac_score_function <- function(calibration.scores,calibration.labels)
{
  conformity_scores <- NULL
  for(i in 1:length(calibration.labels))
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

# set evaluation
evaluate_set_size_for_lac <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:nrow(prediction_set))
    set_size[i] <- sum(!is.na(prediction_set[i,]))
  return(mean(set_size))
}

estimated_qhat <- function(conformity_scores,a)
{
  q <- floor((cal_pct*(N-500) + 1)*(1 - a))/(cal_pct*(N-500) + 1)
  qhat <- quantile(conformity_scores, probs = c(0.25,q))
  return(qhat[2])
}

set.seed(12)

d <- 50
N <- 2000
num_of_classes <- 4
num_of_trials <- 100

alphas <- c(0.01,0.05,0.1,0.15,0.2)

covariance.matrix <- matrix(1/2,nrow = d,ncol = d)

for(i in 1:d)
{
  covariance.matrix[i,i] <- 1
}

X <- mvrnorm(N,rep(0,d),covariance.matrix)
X <- cbind(rep(1,nrow(X)),X)

beta_hat <- c(0,rep(1,3),-sqrt(2),1,rep(0,d-5))


product.function <- X %*% beta_hat
df <- data.frame(cbind(X[,-1],product.function))

df$odds <- exp(product.function)/(1+exp(product.function))
for(i in 1:length(df$odds))
{
  for(class in 0:num_of_classes)
  {
    if(df$odds[i] >= class/num_of_classes && df$odds[i] < (class+1)/num_of_classes)
      df$labels[i] <- class + 1
  }
}

df <- df[,-c(d+1,d+2)]

marginal_coverage_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

cond_coverage_lac <- matrix(0,nrow = num_of_classes,ncol = length(alphas))

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials)
  {
    training.index <- sample(1:N,500)
    df_training <- df[training.index,]
    df_test <- df[-training.index,]
    
    multinomial_model <- multinom(labels~.,data = df_training)
    
    cal_pct <- 0.35
    
    calibration.indexes <- createDataPartition(y=df_test$labels,p=cal_pct,list=FALSE)
    df_calibration <- df[calibration.indexes,]
    df_validation <- df[-calibration.indexes,]
    
    # Calibration data
    
    cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(d+1)],"probs")
    cal_labels <- df_calibration$labels
    
    cal_scores <- NULL
    
    for(i in 1:nrow(cal_probs))
      cal_scores[i] <- cal_probs[i,cal_labels[i]]
    
    # Validation data
    
    val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
    val_labels <- df_validation$labels
    
    lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
    Q_hat_lac <- estimated_qhat(lac_conformity_scores, alphas[alpha])
    
    lac_prediction_set <- naive_lac_prediction_sets(val_probs,val_df,Q_hat_lac)
    lac_coverage <- 0
    
    for(m in 1:nrow(lac_prediction_set))
    {
      if(val_labels[m] %in% lac_prediction_set[m,])
        lac_coverage <- lac_coverage + 1
    }
    
    marginal_coverage_lac[alpha,trial] <- lac_coverage/length(val_labels)
    set_size_lac[alpha,trial] <- evaluate_set_size_for_lac(lac_prediction_set)
  }
}


average_cvg_lac <- NULL
avg_set_size_lac <- NULL

for(alpha in 1:length(alphas))
{
  average_cvg_lac[alpha] <- mean(marginal_coverage_lac[alpha,])
  avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
}


# conditional coverage for Least ambiguous set-valued classifier

for(alpha in 1:length(alphas))
{
  training.index <- sample(1:N,500)
  df_training <- df[training.index,]
  df_test <- df[-training.index,]
  
  multinomial_model <- multinom(labels~.,data = df_training)
  
  cal_pct <- 0.35
  
  calibration.indexes <- createDataPartition(y=df_test$labels,p=cal_pct,list=FALSE)
  df_calibration <- df[calibration.indexes,]
  df_validation <- df[-calibration.indexes,]
  
  # Calibration data
  
  cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(d+1)],"probs")
  cal_labels <- df_calibration$labels
  
  cal_scores <- NULL
  
  for(i in 1:nrow(cal_probs))
    cal_scores[i] <- cal_probs[i,cal_labels[i]]
  
  # Validation data
  
  val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
  val_labels <- df_validation$labels
  
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
  
}


lac_conditional_coverage <- NULL

for(alpha in 1:length(alphas))
  lac_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_lac[,alpha],0)



