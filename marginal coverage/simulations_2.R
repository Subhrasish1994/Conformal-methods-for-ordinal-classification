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

set.seed(1000)

dimensions <- 6
num_of_classes <- 4
num_of_trials <- 10

N <- 1000

alphas <- c(0.05, 0.1, 0.15, 0.2, 0.25)

#---------------------------------Ordinal CDF----------------------------------#

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


#----------------------------------Naive LAC-----------------------------------#

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

#---------------------------------OPI intervals--------------------------------#

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

#------------------------------------OPS---------------------------------------#

opi_prediction_sets <- function(p_values,level.of.significance)
{
  opi.prediction.set <- matrix(0,nrow = nrow(p_values), ncol = num_of_classes)
  for(j in 1:nrow(p_values))
  {
    for(k in 1:num_of_classes)
    {
      if(p_values[j,k] >= level.of.significance)
        opi.prediction.set[j,k] <- k
      else
        opi.prediction.set[j,k] <- NA
    }
  }
  return(opi.prediction.set)
}

#------------------------------------APS---------------------------------------#


#------------------------------------------------------------------------------#

marginal_coverage_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# marginal_coverage_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

set_size_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

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
  else if(labels.generated[i] == 3)
    generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(sqrt(3),dimensions),diag(dimensions))
  else
    generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(2,dimensions),diag(dimensions))
}

df <- data.frame(cbind(generated_data,labels.generated))
names(df)[dimensions + 1] <- "labels"

# logistic regression
multinomial_model <- multinom(labels~.,data = df)
posterior.probability <- predict(multinomial_model,newdata = df, "probs")
df_probs <- data.frame(posterior.probability,df$labels)
names(df_probs)[num_of_classes+1] <- "labels"

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials)
  {
    cal_pct <- 0.1
    cal_indexes <- sample(1:N,cal_pct*N)
    cal_labels <- df_probs$labels[cal_indexes]
    cal_probs <- df_probs[cal_indexes,-(num_of_classes+1)]
    
    val_df <- df_probs[-cal_indexes,]
    val_labels <- val_df$labels
    val_probs <- val_df[,-(num_of_classes+1)]
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] <= cal_probs[,j]) + 1)/(length(cal_probs[,j]) + 1)
      }
    }
    
    ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
    lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
    
    Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
    Q_hat_lac <- estimated_qhat(lac_conformity_scores, alphas[alpha])
    
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
    set_size_lac[alpha,trial] <- evaluate_set_size(lac_prediction_set)
    
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
    
    cbind(ordinal_prediction_intervals,val_labels)
    # ops
    
    ordinal_prediction_sets <- opi_prediction_sets(p.values,alphas[alpha])
    ops_coverage <- 0
    for(m in 1:nrow(ordinal_prediction_sets))
    {
      if(val_labels[m] %in% ordinal_prediction_sets[m,])
        ops_coverage <- ops_coverage + 1
    }
    marginal_coverage_ops[alpha,trial] <- ops_coverage/length(val_labels)
    set_size_ops[alpha,trial] <- evaluate_set_size(ordinal_prediction_sets)
  }
}

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
  
  avg_set_size_cdf[alpha] <- mean(set_size_cdf[alpha,])
  avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

results_average_cvg <- data.frame(cbind(alphas,average_cvg_opi,average_cvg_ops)) 
results_average_set_size <- data.frame(cbind(alphas,avg_set_size_opi,avg_set_size_ops))

fig <- plot_ly(results_average_cvg, x = ~alphas, y = ~average_cvg_opi, name = "<b>Ordinal Prediction Intervals", type = 'scatter', mode = 'lines') 
fig <- fig %>% add_trace(y = ~average_cvg_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
fig %>% layout(yaxis = list(title = "average marginal coverage"))

fig_1 <- plot_ly(results_average_set_size, x = ~alphas)
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
fig_1 %>% layout(yaxis = list(title = "average set size"))

#------------------------------------------------------------------------------#

set.seed(42)

dimensions <- 6
num_of_classes <- 4
N <- 50

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
  else if(labels.generated[i] == 3)
    generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(sqrt(3),dimensions),diag(dimensions))
  else
    generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(2,dimensions),diag(dimensions))
}

df <- data.frame(cbind(generated_data,labels.generated))
names(df)[dimensions + 1] <- "labels"

# logistic regression
multinomial_model <- multinom(labels~.,data = df)
posterior.probability <- predict(multinomial_model,newdata = df, "probs")
df_probs <- data.frame(posterior.probability,df$labels)
names(df_probs)[num_of_classes+1] <- "labels"

cal_pct <- 0.5
cal_indexes <- sample(1:N,cal_pct*N)
# sort(cal_indexes)
cal_labels <- df_probs$labels[cal_indexes]
cal_probs <- df_probs[cal_indexes,-(num_of_classes+1)]

val_df <- df_probs[-cal_indexes,]
val_labels <- val_df$labels
val_probs <- val_df[,-(num_of_classes+1)]

p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)

for(i in 1:nrow(val_probs))
{
  for(j in 1:num_of_classes)
  {
    p.values[i,j] <- (sum(val_probs[i,j] >= cal_probs[,j]) + 1)/(length(cal_probs[,j]) + 1)
  }
}

count <- 0
for(i in 1:nrow(p.values))
{
  idx <- which(p.values[i,] > 1)
  print(idx)
  if(val_labels[i] %in% idx)
    count <- count + 1
}

count/length(val_labels)
