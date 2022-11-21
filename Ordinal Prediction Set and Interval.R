# Ordinal Prediction Set & Interval

# Based on multiple hypothesis testing formulation

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

opi_prediction_intervals <- function(p_values,level.of.significance)
{
  lower_edge <- NULL
  upper_edge <- NULL
  
  for(i in 1:nrow(p_values))
  {
    lower_edge[i] <- min(which(p_values[i,] > level.of.significance))
    upper_edge[i] <- max(which(p_values[i,] > level.of.significance))
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

evaluate_set_size_ops <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:length(prediction_set))
    set_size[i] <- length(as.vector(prediction_set[[i]]))
  return(mean(set_size))
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

marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

cond_coverage_opi <- matrix(0,nrow = num_of_classes,ncol = length(alphas))
cond_coverage_ops <- matrix(0,nrow = num_of_classes,ncol = length(alphas))

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials){
    
    training.index <- sample(1:N,500)
    df_training <- df[training.index,]
    df_test <- df[-training.index,]
    
    multinomial_model <- multinom(labels~.,data = df_training)
    
    cal_pct <- 0.35
    
    calibration.indexes <- createDataPartition(y=df_test$labels,p=cal_pct,list=FALSE)
    df_calibration <- df[calibration.indexes,]
    df_validation <- df[-calibration.indexes,]
    
    cal_probs <- predict(multinomial_model,newdata = df_calibration[,-(d+1)],"probs")
    cal_labels <- df_calibration$labels
    
    cal_scores <- NULL
    
    for(i in 1:nrow(cal_probs))
      cal_scores[i] <- cal_probs[i,cal_labels[i]]
    
    # Validation data
    
    val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
    val_labels <- df_validation$labels
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] >= cal_scores)+1)/(length(cal_scores) + 1)
      }
    }
    
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

average_cvg_opi <- NULL
average_cvg_ops <- NULL

avg_set_size_opi <- NULL
avg_set_size_ops <- NULL

for(alpha in 1:length(alphas))
{
  average_cvg_opi[alpha] <- mean(marginal_coverage_opi[alpha,])
  average_cvg_ops[alpha] <- mean(marginal_coverage_ops[alpha,])

  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,][is.finite(set_size_cdf[alpha,])])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

# conditional coverage for  Ordinal Prediction Interval and Sets

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
  
  # Validation data
  
  val_probs <- predict(multinomial_model,newdata = df_validation[,-(d+1)],"probs")
  val_labels <- df_validation$labels
  
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

opi_conditional_coverage <- NULL
ops_conditional_coverage <- NULL

for(alpha in 1:length(alphas))
{
  opi_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_opi[,alpha],0)
  ops_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_ops[,alpha],0)
}

library(ggpubr)

consolidated_data <- read.csv("/Users/subhrasishchakraborty/Desktop/consol.csv",header = TRUE)

plt_1 <- ggplot(consolidated_data,aes(x=alpha))+geom_line(aes(y = marginal.coverage,color=method))
plt_1 <- plt_1 + facet_grid(.~dimensions) + ylab("marginal coverage")

plt_2 <- ggplot(consolidated_data,aes(x=alpha))+geom_line(aes(y = set.size,color=method))
plt_2 <- plt_2 + facet_grid(.~dimensions) + ylab("set size")

plt_3 <- ggplot(consolidated_data,aes(x=alpha))+geom_line(aes(y = conditional.coverage,color=method))
plt_3 <- plt_3 + facet_grid(.~dimensions) + ylab("conditional coverage")

ggarrange(plt_1,plt_2,plt_3,ncol = 1, nrow = 3, common.legend = TRUE, legend = "right")
