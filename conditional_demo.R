# conditional coverage with choice of base classifier to be logistic reg

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

set.seed(12)

dimensions <- 3
num_of_classes <- 3
num_of_trials <- 100

# alphas <- c(0.01,0.05,0.1,0.15,0.2)

alpha <- 0.1

N <- 1000 # 50 training observations

marginal_coverage_cdf <- matrix(0, nrow = num_of_classes, ncol = num_of_trials)
marginal_coverage_lac <- matrix(0, nrow = num_of_classes, ncol = num_of_trials)
marginal_coverage_opi <- matrix(0, nrow = num_of_classes, ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = num_of_classes, ncol = num_of_trials)

# set_size_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)


ordinal_CDF_score_function <- function(calibration.scores,calibration.labels)
{
  conformity_scores <- rep(0,nrow(calibration.scores))
  for(i in 1:nrow(calibration.scores))
  {
    pos.max <- which.is.max(calibration.scores[i,])
    conformity_scores[i] <- abs(sum(calibration.scores[i,1:pos.max]) - sum(calibration.scores[i,1:calibration.labels[i]]))
  }
  return(conformity_scores)
}

ordinal_CDF_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  prediction.interval <- matrix(0,nrow = nrow(validation.scores), ncol = num_of_classes)
  
  lower.edge <- NULL
  upper.edge <- NULL
  
  for(j in 1:nrow(validation.scores))
  {
    pos <- which.is.max(validation.scores[j,])
    for(k in 1:num_of_classes)
    {
      if(sum(validation.scores[j,1:k]) >= validation.scores[j,pos] - qhat && sum(validation.scores[j,1:k]) <= validation.scores[j,pos] + qhat)
        prediction.interval[j,k] <- k
      else
        prediction.interval[j,k] <- NA
    }
  }
  for(j in 1:nrow(prediction.interval))
  {
    lower.edge[j] <- min(prediction.interval[j,],na.rm = TRUE)
    upper.edge[j] <- max(prediction.interval[j,], na.rm = TRUE)
  }
  return(cbind(lower.edge,upper.edge))
}


lac_score_function <- function(calibration.scores,calibration.labels)
{
  conformity_scores <- rep(0,nrow(calibration.scores))
  for(i in 1:nrow(calibration.scores))
  {
    pos.max <- which.is.max(calibration.scores[i,])
    conformity_scores[i] <- 1 - sum(calibration.scores[i,1:calibration.labels[i]])
  }
  return(conformity_scores)
}

naive_lac_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  lac.prediction.set <- matrix(0,nrow = nrow(validation.scores), ncol = num_of_classes)
  for(j in 1:nrow(validation.scores))
  {
    for(k in 1:num_of_classes)
    {
      if(sum(validation.scores[j,1:k]) >= 1 - qhat)
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
  q <- floor((cal_pct*(N-50) + 1)*(1 - a))/(cal_pct*(N-50) + 1)
  qhat <- quantile(conformity_scores, probs = c(0.25,q))
  return(qhat[2])
}

for(class in seq(1,num_of_classes))
{
  for(trial in 1:num_of_trials)
  {
    
    labels.generated <- rep(0,N)
    for(i in 1:N)
      labels.generated[i] <- sample(1:num_of_classes,prob = rep(1/num_of_classes,num_of_classes))
    generated_data <- matrix(0,ncol = dimensions, nrow = N)
    for(i in 1:N)
    {
      if(labels.generated[i] == 1)
        generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(-1,dimensions),diag(dimensions))
      else if(labels.generated[i] == 2)
        generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(1,dimensions),diag(dimensions))
      else
        generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(sqrt(3),dimensions),diag(dimensions))
    }
    
    df <- data.frame(cbind(generated_data,labels.generated))
    names(df)[dimensions + 1] <- "labels"
    
    
    training.index <- sample(1:N,50)
    df_training <- df[training.index,]
    df_test <- df[-training.index,]
    
    # logistic regression
    
    multinomial_model <- multinom(labels~.,data = df_training)
    
    cal_pct <- 0.2
    
    calibration.indexes <- createDataPartition(y=df_test$labels,p=cal_pct,list=FALSE)
    df_calibration <- df[calibration.indexes,]
    df_validation <- df[-calibration.indexes,]
    
    # Calibration data
    
    cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(dimensions+1)],"probs")
    cal_labels <- df_calibration$labels
    
    # cal_scores <- NULL
    # 
    # for(i in 1:nrow(cal_probs))
    #   cal_scores[i] <- cal_probs[i,cal_labels[i]]
    # 
    # Validation data
    
    val_probs <- predict(multinomial_model,newdata = df_validation[,-(dimensions+1)],"probs")
    val_labels <- df_validation$labels
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] >= cal_probs[which(cal_labels == j),j])+1)/(length(cal_probs[which(cal_labels == j),j]) + 1)
      }
    }
    
    ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
    lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
    
    Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alpha)
    Q_hat_lac <- estimated_qhat(lac_conformity_scores, alpha)
    
    cdf_prediction_interval <- ordinal_CDF_prediction_sets(val_probs,val_labels,Q_hat_cdf)
    cdf_cond_cvg <- 0
    
    for(m in 1:nrow(cdf_prediction_interval))
    {
      if(val_labels[m] == class)
      {
        if(cdf_prediction_interval[m,1] <= val_labels[m] && val_labels[m] <= cdf_prediction_interval[m,2] )
          cdf_cond_cvg <- cdf_cond_cvg + 1
      }
    }
    marginal_coverage_cdf[class,trial] <- cdf_cond_cvg/length(which(val_labels == class))
    # set_size_cdf[alpha,trial] <- mean(cdf_prediction_interval[,2] - cdf_prediction_interval[,1] + 1)
    
    # naive lac
    
    lac_prediction_set <- naive_lac_prediction_sets(val_probs,val_df,Q_hat_lac)
    lac_cond_cvg <- 0
    
    for(m in 1:nrow(lac_prediction_set))
    {
      if(val_labels[m] == class)
      {
          if(val_labels[m] %in% lac_prediction_set[m,])
              lac_cond_cvg <- lac_cond_cvg + 1
      }
    }
    
    marginal_coverage_lac[class,trial] <- lac_cond_cvg/length(which(val_labels == class))
    # set_size_lac[alpha,trial] <- evaluate_set_size_for_lac(lac_prediction_set)
    
    # opi
    
    ordinal_prediction_intervals <- opi_prediction_intervals(p.values,alpha)
    opi_cond_cvg <- 0
    
    for(m in 1:nrow(ordinal_prediction_intervals))
    {
      if(val_labels[m] == class)
      {
        if(ordinal_prediction_intervals[m,1] <= val_labels[m] && val_labels[m] <= ordinal_prediction_intervals[m,2] )
             opi_cond_cvg <- opi_cond_cvg + 1
      }
    }
    marginal_coverage_opi[class,trial] <- opi_cond_cvg/length(which(val_labels == class))  
    # set_size_opi[alpha,trial] <- mean(ordinal_prediction_intervals[,2] - ordinal_prediction_intervals[,1] + 1)
    
    # ops
    
    ordinal_prediction_sets <- opi_prediction_sets(p.values,alpha)
    ops_cond_cvg <- 0
    
    for(m in 1:length(ordinal_prediction_sets))
    {
      if(val_labels[m] == class)
      {
      if(val_labels[m] %in% as.vector(ordinal_prediction_sets[[m]]))
        ops_cond_cvg <- ops_cond_cvg + 1
      }
    }
    marginal_coverage_ops[class,trial] <- ops_cond_cvg/length(which(val_labels == class))
    # set_size_ops[alpha,trial] <- evaluate_set_size_ops(ordinal_prediction_sets)
  }
}

plt <- plot_ly(alpha = 0.6)
plt <- plt %>% add_histogram(x=~marginal_coverage_cdf[1,],name = "Ordinal CDF")
plt <- plt %>% add_histogram(x=~marginal_coverage_lac[1,],name = "Naive LAC")
plt <- plt %>% add_histogram(x=~marginal_coverage_opi[1,],name = "Ordinal Prediction Interval")
plt <- plt %>% add_histogram(x=~marginal_coverage_ops[1,],name = "Ordinal Prediction Set")
plt <- plt %>% layout(barmode = "overlay",title = "conditional coverage on class 1",xaxis=list(title = "conditional coverage"))

plt <- plot_ly(alpha = 0.6)
plt <- plt %>% add_histogram(x=~marginal_coverage_cdf[2,],name = "Ordinal CDF")
plt <- plt %>% add_histogram(x=~marginal_coverage_lac[2,],name = "Naive LAC")
plt <- plt %>% add_histogram(x=~marginal_coverage_opi[2,],name = "Ordinal Prediction Interval")
plt <- plt %>% add_histogram(x=~marginal_coverage_ops[2,],name = "Ordinal Prediction Set")
plt <- plt %>% layout(barmode = "overlay",title = "conditional coverage on class 2",xaxis=list(title = "conditional coverage"))

plt <- plot_ly(alpha = 0.6)
plt <- plt %>% add_histogram(x=~marginal_coverage_cdf[3,],name = "Ordinal CDF")
plt <- plt %>% add_histogram(x=~marginal_coverage_lac[3,],name = "Naive LAC")
plt <- plt %>% add_histogram(x=~marginal_coverage_opi[3,],name = "Ordinal Prediction Interval")
plt <- plt %>% add_histogram(x=~marginal_coverage_ops[3,],name = "Ordinal Prediction Set")
plt <- plt %>% layout(barmode = "overlay",title = "conditional coverage on class 3",xaxis=list(title = "conditional coverage"))
