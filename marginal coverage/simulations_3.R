# ipak <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg)) 
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, require, character.only = TRUE)
# }
# 
# packages <- c("mvtnorm","ggplot2","caret","ipred","devtools","sm","hrbrthemes","tidyr","dplyr","nnet","viridis","plotly","scatterplot3d", "patchwork")
# ipak(packages)
# 
# library(MASS)
# library(mvtnorm)
# library(ggplot2)
# library(caret)
# library(ipred)
# library(devtools)
# library(sm)
# library(hrbrthemes)
# library(tidyr)
# library(dplyr)
# library(plotly)
# library(viridis)
# library(scatterplot3d)
# library(patchwork)
# 
# set.seed(12)
# 
# dimensions <- 6
# num_of_classes <- 4
# num_of_trials <- 10
# 
# alphas <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
# 
# N <- 1050 # 50 training observations
# 
# labels.generated <- rep(0,N)
# for(i in 1:N)
#   labels.generated[i] <- sample(1:num_of_classes,prob = rep(1/num_of_classes,num_of_classes))
# generated_data <- matrix(0,ncol = dimensions, nrow = N)
# for(i in 1:N)
# {
#   if(labels.generated[i] == 1)
#     generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(-1,dimensions),diag(dimensions))
#   else if(labels.generated[i] == 2)
#     generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(1,dimensions),diag(dimensions))
#   else if(labels.generated[i] == 3)
#     generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(sqrt(3),dimensions),diag(dimensions))
#   else
#     generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(2,dimensions),diag(dimensions))
# }
# 
# df <- data.frame(cbind(generated_data,labels.generated))
# names(df)[dimensions + 1] <- "labels"
# 
# training.index <- sample(1:N,50)
# df_training <- df[training.index,]
# df_test <- df[-training.index,]
# 
# # logistic regression
# 
# multinomial_model <- multinom(labels~.,data = df_training)
# posterior.probability <- predict(multinomial_model,newdata = df_test, "probs")
# df_probs <- data.frame(posterior.probability,df_test$labels)
# names(df_probs)[num_of_classes+1] <- "labels"
# 
# #---------------------------------Ordinal CDF----------------------------------#
# 
# ordinal_CDF_score_function <- function(calibration.scores,calibration.labels)
# {
#   conformity_scores <- rep(0,nrow(calibration.scores))
#   for(i in 1:nrow(calibration.scores))
#   {
#     pos.max <- which.is.max(calibration.scores[i,])
#     conformity_scores[i] <- abs(sum(calibration.scores[i,1:pos.max]) - sum(calibration.scores[i,1:calibration.labels[i]]))
#   }
#   return(conformity_scores)
# }
# 
# ordinal_CDF_prediction_sets <- function(validation.scores,validation.labels,qhat)
# {
#   prediction.interval <- matrix(0,nrow = nrow(validation.scores), ncol = num_of_classes)
#   
#   lower.edge <- NULL
#   upper.edge <- NULL
#   
#   for(j in 1:nrow(validation.scores))
#   {
#     pos <- which.is.max(validation.scores[j,])
#     for(k in 1:num_of_classes)
#     {
#       if(sum(validation.scores[j,1:k]) >= validation.scores[j,pos] - qhat && sum(validation.scores[j,1:k]) <= validation.scores[j,pos] + qhat)
#         prediction.interval[j,k] <- k
#       else
#         prediction.interval[j,k] <- NA
#     }
#   }
#   for(j in 1:nrow(prediction.interval))
#   {
#     lower.edge[j] <- min(prediction.interval[j,],na.rm = TRUE)
#     upper.edge[j] <- max(prediction.interval[j,], na.rm = TRUE)
#   }
#   return(cbind(lower.edge,upper.edge))
# }
# 
# 
# #----------------------------------Naive LAC-----------------------------------#
# 
# lac_score_function <- function(calibration.scores,calibration.labels)
# {
#   conformity_scores <- rep(0,nrow(calibration.scores))
#   for(i in 1:nrow(calibration.scores))
#   {
#     pos.max <- which.is.max(calibration.scores[i,])
#     conformity_scores[i] <- 1 - sum(calibration.scores[i,1:calibration.labels[i]])
#   }
#   return(conformity_scores)
# }
# 
# naive_lac_prediction_sets <- function(validation.scores,validation.labels,qhat)
# {
#   lac.prediction.set <- matrix(0,nrow = nrow(validation.scores), ncol = num_of_classes)
#   for(j in 1:nrow(validation.scores))
#   {
#     for(k in 1:num_of_classes)
#     {
#       if(sum(validation.scores[j,1:k]) >= 1 - qhat)
#         lac.prediction.set[j,k] <- k
#       else
#         lac.prediction.set[j,k] <- NA
#     }
#   }
#   return(lac.prediction.set)
# }
# 
#---------------------------------OPI intervals--------------------------------#

# opi_prediction_intervals <- function(p_values,level.of.significance)
# {
#   lower_edge <- NULL
#   upper_edge <- NULL
# 
#   for(i in 1:nrow(p_values))
#   {
#     lower_edge[i] <- min(which(p_values[i,] > level.of.significance))
#     upper_edge[i] <- max(which(p_values[i,] > level.of.significance))
#   }
#   return(cbind(lower_edge,upper_edge))
# }
# 
# #------------------------------------OPS---------------------------------------#
# 
# opi_prediction_sets <- function(p_values,level.of.significance)
# {
#   opi.prediction.set <- matrix(0,nrow = nrow(p_values), ncol = num_of_classes)
#   for(j in 1:nrow(p_values))
#   {
#     for(k in 1:num_of_classes)
#     {
#       if(p_values[j,k] > level.of.significance)
#         opi.prediction.set[j,k] <- k
#       else
#         opi.prediction.set[j,k] <- NA
#     }
#   }
#   return(opi.prediction.set)
# }

# #------------------------------------------------------------------------------#
# 
# evaluate_set_size <- function(prediction_set)
# {
#   set_size <- NULL
#   for(i in 1:nrow(prediction_set))
#     set_size[i] <- sum(!is.na(prediction_set[i,]))
#   return(mean(set_size))
# }
# 
# estimated_qhat <- function(conformity_scores,a)
# {
#   q <- floor((cal_pct*N + 1)*(1 - a))/(cal_pct*N + 1)
#   qhat <- quantile(conformity_scores, probs = c(0.25,q))
#   return(qhat[2])
# }
# 
# marginal_coverage_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# marginal_coverage_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# # marginal_coverage_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# 
# set_size_cdf <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_lac <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# # set_size_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
# 
# for(alpha in 1:length(alphas))
# {
#   for(trial in 1:num_of_trials)
#   {
#     cal_pct <- 0.2
#     
#     calibration.indexes <- createDataPartition(y=df_probs$labels,p=cal_pct,list=FALSE)
#     df_cal <- df_probs[calibration.indexes,]
#     df_val <- df_probs[-calibration.indexes,]
#     
#     cal_labels <- df_cal$labels
#     cal_probs <- df_cal[,-(num_of_classes+1)]
#     
#     val_labels <- df_val$labels
#     val_probs <- df_val[,-(num_of_classes+1)]
#     
#     p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
#     
#     for(i in 1:nrow(val_probs))
#     {
#       for(j in 1:num_of_classes)
#       {
#         p.values[i,j] <- (sum(val_probs[i,j] <= cal_probs[,j])+1)/(length(cal_probs[,j]) + 1)
#       }
#     }
#     ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
#     lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
#     
#     Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
#     Q_hat_lac <- estimated_qhat(lac_conformity_scores, alphas[alpha])
#     
#     cdf_prediction_interval <- ordinal_CDF_prediction_sets(val_probs,val_labels,Q_hat_cdf)
#     cdf_coverage <- 0
#     
#     for(m in 1:nrow(cdf_prediction_interval))
#     {
#       if(cdf_prediction_interval[m,1] <= val_labels[m] && val_labels[m] <= cdf_prediction_interval[m,2] )
#         cdf_coverage <- cdf_coverage + 1
#     }
#     marginal_coverage_cdf[alpha,trial] <- cdf_coverage/length(val_labels)
#     set_size_cdf[alpha,trial] <- mean(cdf_prediction_interval[,2] - cdf_prediction_interval[,1] + 1)
#     
#     # naive lac
#     
#     lac_prediction_set <- naive_lac_prediction_sets(val_probs,val_df,Q_hat_lac)
#     lac_coverage <- 0
#     for(m in 1:nrow(lac_prediction_set))
#     {
#       if(val_labels[m] %in% lac_prediction_set[m,])
#         lac_coverage <- lac_coverage + 1
#     }
#     
#     marginal_coverage_lac[alpha,trial] <- lac_coverage/length(val_labels)
#     set_size_lac[alpha,trial] <- evaluate_set_size(lac_prediction_set)
#     
#     # opi
#     
#     ordinal_prediction_intervals <- opi_prediction_intervals(p.values,alphas[alpha])
#     opi_coverage <- 0
#     for(m in 1:nrow(ordinal_prediction_intervals))
#     {
#       if(ordinal_prediction_intervals[m,1] <= val_labels[m] && val_labels[m] <= ordinal_prediction_intervals[m,2] )
#         opi_coverage <- opi_coverage + 1
#     }
#     marginal_coverage_opi[alpha,trial] <- opi_coverage/length(val_labels)  
#     set_size_opi[alpha,trial] <- mean(ordinal_prediction_intervals[,2] - ordinal_prediction_intervals[,1] + 1)
#     
#     # ops
#     
#     ordinal_prediction_sets <- opi_prediction_sets(p.values,alphas[alpha])
#     ops_coverage <- 0
#     for(m in 1:nrow(ordinal_prediction_sets))
#     {
#       if(val_labels[m] %in% ordinal_prediction_sets[m,])
#         ops_coverage <- ops_coverage + 1
#     }
#     marginal_coverage_ops[alpha,trial] <- ops_coverage/length(val_labels)
#     set_size_ops[alpha,trial] <- evaluate_set_size(ordinal_prediction_sets)
#     
#   }
# }
# 
# average_cvg_cdf <- NULL
# average_cvg_lac <- NULL
# average_cvg_opi <- NULL
# average_cvg_ops <- NULL
# 
# avg_set_size_cdf <- NULL
# avg_set_size_lac <- NULL
# avg_set_size_opi <- NULL
# avg_set_size_ops <- NULL
# 
# for(alpha in 1:length(alphas))
# {
#   average_cvg_cdf[alpha] <- mean(marginal_coverage_cdf[alpha,])
#   average_cvg_lac[alpha] <- mean(marginal_coverage_lac[alpha,])
#   average_cvg_opi[alpha] <- mean(marginal_coverage_opi[alpha,])
#   average_cvg_ops[alpha] <- mean(marginal_coverage_ops[alpha,])
#   
#   avg_set_size_cdf[alpha] <- mean(set_size_cdf[alpha,])
#   avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
#   avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
#   avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
# }
# 
# results_average_cvg <- data.frame(cbind(alphas,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops)) 
# results_average_set_size <- data.frame(cbind(alphas, avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))
# 
# fig <- plot_ly(results_average_cvg, x = ~alphas, y = ~average_cvg_cdf, name = "<b>Ordinal CDF", type = 'scatter', mode = 'lines') 
# fig <- fig %>% add_trace(y = ~average_cvg_lac, name = "<b>Naive LAC", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig %>% layout(yaxis = list(title = "average marginal coverage"))
# 
# fig_1 <- plot_ly(results_average_set_size, x = ~alphas) 
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_cdf, name = "<b>Ordinal CDF", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_lac, name = "<b>Naive LAC", connectgaps = TRUE)
# # fig_1 <- fig_1 %>% add_lines(y = ~V4, name = "<b>APS", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig_1 %>% layout(yaxis = list(title = "average set size"))
# 
# par(mfrow = c(2,2))
# 
# hist(p.values[,1],main = "p1")
# hist(p.values[,2],main = "p2")
# hist(p.values[,3],main = "p3")
# hist(p.values[,4],main = "p4")

# count <- 0
# prediction.sets <- list()
# prediction.sets <- vector("list",length = length(val_labels))
# for(i in 1:nrow(p.values))
# {
#    idx <- which(p.values[i,] >= 0.2)
#    prediction.sets[[i]] <- idx
# }
# 
# length(prediction.sets)
# count <- 0
# for(i in 1:length(val_labels))
# {
#   if(val_labels[i] %in% as.vector(prediction.sets[[i]]))
#     count <- count + 1
# }
# count/length(val_labels)

#==============================================================================#
#==============================================================================#
#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#
#==============================================================================#
#==============================================================================#

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

set.seed(12)

dimensions <- 6
num_of_classes <- 4
num_of_trials <- 10

alphas <- c(0.05,0.1,0.15,0.2,0.25)

N <- 1050 # 50 training observations

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
  # opi.prediction.set <- matrix(0,nrow = nrow(p_values), ncol = num_of_classes)
  opi.prediction.set = list()
  opi.prediction.set = vector("list",length = nrow(p_values))
  for(i in 1:nrow(p.values))
  {
    idx <- which(p.values[i,] >= level.of.significance)
    opi.prediction.set[[i]] <- idx
  }
  return(opi.prediction.set)
}

#------------------------------------------------------------------------------#

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

for(alpha in 1:length(alphas))
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
      else if(labels.generated[i] == 3)
        generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(sqrt(3),dimensions),diag(dimensions))
      else
        generated_data[i,] <- mvrnorm(1,(1/sqrt(dimensions))*rep(2,dimensions),diag(dimensions))
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
    
    # Validation data
    
    val_probs <- predict(multinomial_model,newdata = df_validation[,-(dimensions+1)],"probs")
    val_labels <- df_validation$labels
    
    write.csv(data.frame(cbind(cal_probs,cal_labels)),"/Users/subhrasishchakraborty/Desktop/calibration_data.csv")
    write.csv(data.frame(cbind(val_probs,val_labels)), "/Users/subhrasishchakraborty/Desktop/validation_data.csv")
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] <= cal_probs[,j])+1)/(length(cal_probs[,j]) + 1)
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
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

average_cvg_aps <- c(0.935, 0.905,0.867,0.786,0.717)
avg_set_size_aps <- c(3, 2.53, 2.30,2.02, 1.81)
results_average_cvg <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

# fig <- plot_ly(results_average_cvg, x = ~alphas, y = ~average_cvg_cdf, name = "<b>Ordinal CDF", type = 'scatter', mode = 'lines')
# fig <- fig %>% add_trace(y = ~average_cvg_lac, name = "<b>Naive LAC", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig %>% layout(yaxis = list(title = "average marginal coverage"))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

ggplt_1

# fig_1 <- plot_ly(results_average_set_size, x = ~alphas)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_cdf, name = "<b>Ordinal CDF", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_lac, name = "<b>Naive LAC", connectgaps = TRUE)
# # fig_1 <- fig_1 %>% add_lines(y = ~V4, name = "<b>APS", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig_1 %>% layout(yaxis = list(title = "average set size"))

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_2 

ggplt_1 / ggplt_2

#==============================================================================#
#==============================================================================#
#==============================================================================#

ordinal_aps_function <- function(validation.scores,validation.labels,lambda)
{
  starting.point <- NULL
  ending.point <- NULL
  
  for(i in 1:nrow(validation.scores))
  {
    maximum.position <- which.is.max(validation.scores[i,])
    begin <- max(maximum.position-1,1)
    end <- min(maximum.position+1,num_of_classes)
    if(begin<end)
    {
      q <- 0
      while(q <= lambda)
      {
        q = sum(validation.scores[i,begin:end])
        end <- min(end + 1, num_of_classes)
      }
      starting.point[i] <- begin
      ending.point[i] <- min(end+1,num_of_classes)
    }
    else
    {
      starting.point[i] <- maximum.position
      ending.point[i] <- maximum.position
    }
  }
  return(cbind(starting.point,ending.point))
}

pred_int <- ordinal_aps_function(val_probs,val_labels,0.8)

count <- 0
for(m in 1:length(val_labels))
{
  if(val_labels[m] >= pred_int[m,1] && val_labels[m] <= pred_int[m,2])
    count <- count + 1
}
count/length(val_labels)