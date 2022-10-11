# ipak <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg)) 
#     install.packages(new.pkg, dependencies = TRUE)
#   sapply(pkg, require, character.only = TRUE)
# }
# 
# packages <- c("mvtnorm","ggplot2","caret","ipred","devtools","sm","hrbrthemes","tidyr","dplyr","nnet","viridis","plotly","scatterplot3d", "patchwork")
# ipak(packages)

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
library(PerformanceAnalytics)

getwd()
setwd("/Users/subhrasishchakraborty/Desktop/dataset")

traffic_accident_leeds <- read.csv("Traffic accidents_2019_Leeds.csv",header = TRUE)

# View(traffic_accident_leeds)
earthquake_nz_intensity <- read.csv("earthquakes_NZ.csv", header = TRUE)

names(earthquake_nz_intensity)
# install.packages("PerformanceAnalytics")

# chart.Correlation(traffic_accident_leeds, histogram = TRUE, method = "pearson")
# chart.Correlation(earthquake_nz_intensity, histogram = TRUE, method = "pearson")

# View(earthquake_nz_intensity)
earthquake_nz_intensity <- earthquake_nz_intensity[,-1]

# pairs(earthquake_nz_intensity)
for(label in 1:nrow(earthquake_nz_intensity))
{
  if(earthquake_nz_intensity$magnitude[label] < 1.5)
    earthquake_nz_intensity$labels[label] <- 3
  else if(1.5 <= earthquake_nz_intensity$magnitude[label] && earthquake_nz_intensity$magnitude[label] <= 3)
    earthquake_nz_intensity$labels[label] <- 2
  else
    earthquake_nz_intensity$labels[label] <- 1
}
model_accident <- multinom(Casualty.Severity~.,data = traffic_accident_leeds)
# model_earthquake <- multinom(labels~latitude+longitude+depth,data = earthquake_nz_intensity)
posterior.probability <- predict(model_accident,newdata = traffic_accident_leeds, "probs")
# posterior.probability <- predict(model_earthquake,newdata = earthquake_nz_intensity, "probs")
df_probs_1 <- data.frame(posterior.probability,traffic_accident_leeds$Casualty.Severity)
# df_probs_1 <- data.frame(posterior.probability,earthquake_nz_intensity$labels)
names(df_probs_1)[length(unique(traffic_accident_leeds$Casualty.Severity))+1] <- "labels"
# names(df_probs_1)[length(unique(earthquake_nz_intensity$labels))+1] <- "labels"

cal_pct <- 0.4

N <- nrow(traffic_accident_leeds)

alphas <- seq(0.01,0.15,0.01)

# N <- nrow(earthquake_nz_intensity)

num_of_classes <- length(unique(traffic_accident_leeds$Casualty.Severity))

# num_of_classes <- length(unique(earthquake_nz_intensity$labels))
num_of_trials <- 10

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
for(trial in 1:num_of_trials)
{

cal_indexes <- sample(1:N,cal_pct*N)
cal_labels <- df_probs_1$labels[cal_indexes]
cal_probs <- df_probs_1[cal_indexes,-(num_of_classes+1)]

val_df <- df_probs_1[-cal_indexes,]
val_labels <- val_df$labels
val_probs <- val_df[,-(num_of_classes+1)]

p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)

ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)

for(i in 1:nrow(val_probs))
{
  for(j in 1:num_of_classes)
  {
    p.values[i,j] <- (sum(val_probs[i,j] <= cal_probs[,j]) + 1)/(length(cal_probs[,j]) + 1)
  }
}


for(alpha in 1:length(alphas))
{
  # ordinal cdf
  
  Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
  
  Q_hat_lac <- estimated_qhat(lac_conformity_scores,alphas[alpha])
  
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
  
  # aps
  
  # aps_prediction_interval <- ordinal_aps_function(val_probs,val_labels,1-alphas[alpha])
  # aps_coverage <- 0
  # for(m in 1:nrow(aps_prediction_interval))
  # {
  #   if(aps_prediction_interval[m,1] <= val_labels[m] && val_labels[m] <= aps_prediction_interval[m,2] )
  #     aps_coverage <- aps_coverage + 1
  # }
  # marginal_coverage_aps[alpha,trial] <- aps_coverage/length(val_labels)
  # set_size_aps[alpha,trial] <- mean(aps_prediction_interval[,2] - aps_prediction_interval[,1] + 1)
  
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

results_average_cvg <- data.frame(cbind(alphas,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops)) 
results_average_set_size <- data.frame(cbind(alphas, avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

fig <- plot_ly(results_average_cvg, x = ~alphas, y = ~average_cvg_cdf, name = "<b>Ordinal CDF", type = 'scatter', mode = 'lines') 
fig <- fig %>% add_trace(y = ~average_cvg_lac, name = "<b>Naive LAC", connectgaps = TRUE)
fig <- fig %>% add_trace(y = ~average_cvg_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
fig <- fig %>% add_trace(y = ~average_cvg_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
fig %>% layout(yaxis = list(title = "average marginal coverage"))

fig_1 <- plot_ly(results_average_set_size, x = ~alphas) 
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_cdf, name = "<b>Ordinal CDF", connectgaps = TRUE)
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_lac, name = "<b>Naive LAC", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~V4, name = "<b>APS", connectgaps = TRUE)
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
fig_1 %>% layout(yaxis = list(title = "average set size"))

# marginal_coverage_rs <- matrix(0,nrow = length(alphas), ncol = 5)
# 
# for(alpha in 1:length(alphas))
# {
#   marginal_coverage_rs[alpha,1] <- mean(marginal_coverage_cdf[alpha,])
#   marginal_coverage_rs[alpha,2] <- mean(marginal_coverage_lac[alpha,])
#   marginal_coverage_rs[alpha,3] <- mean(marginal_coverage_aps[alpha,])
#   marginal_coverage_rs[alpha,4] <- mean(marginal_coverage_opi[alpha,])
#   marginal_coverage_rs[alpha,5] <- mean(marginal_coverage_ops[alpha,])
# }
# 
# marginal_coverage_rs <- data.frame(cbind(alphas,marginal_coverage_rs))
# View(marginal_coverage_rs)

# fig <- plot_ly(marginal_coverage_rs, x = ~alphas, y = ~V2, name = "<b>Ordinal CDF", type = 'scatter', mode = 'lines') 
# fig <- fig %>% add_trace(y = ~V3, name = "<b>Naive LAC", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~V4, name = "<b>APS", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~V5, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~V6, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig %>% layout(yaxis = list(title = "marginal coverage"))
# 
# set_size_rs <- matrix(0,nrow = length(alphas), ncol = 5)
# for(alpha in 1:length(alphas))
# {
#   set_size_rs[alpha,1] <- mean(set_size_cdf[alpha,])
#   set_size_rs[alpha,2] <- mean(set_size_lac[alpha,])
#   set_size_rs[alpha,3] <- mean(set_size_aps[alpha,])
#   set_size_rs[alpha,4] <- mean(set_size_opi[alpha,])
#   set_size_rs[alpha,5] <- mean(set_size_ops[alpha,])
# }
# set_size_rs <- data.frame(cbind(alphas,set_size_rs))
# # View(set_size_rs)
# 
# fig_1 <- plot_ly(set_size_rs, x = ~alphas) 
# fig_1 <- fig_1 %>% add_lines(y = ~V2, name = "<b>Ordinal CDF", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~V3, name = "<b>Naive LAC", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~V4, name = "<b>APS", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~V5, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~V6, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig_1 %>% layout(yaxis = list(title = "average set size"))

getwd()
setwd("/Users/subhrasishchakraborty/Desktop")

calibration_data <- read.csv("cal_df.csv", header = TRUE)
calibration_data[,1] <- calibration_data[,1] + 1
names(calibration_data)[1] <- "labels"

validation_data <- read.csv("val_df.csv", header = TRUE)
validation_data[,1] <- validation_data[,1] + 1
names(validation_data)[1] <- "labels"

alphas <- c(0.01, 0.05, 0.1, 0.15)

num_of_trials <- 10
num_of_classes <- 4

marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

cal_labels <- calibration_data$labels
cal_probs <- calibration_data[,-1]

val_labels <- validation_data$labels
val_probs <- validation_data[,-1]

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials)
  {
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] >= cal_probs[,j]) + 1)/(length(cal_probs[,j]) + 1)
      }
    }
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

average_cvg_opi <- NULL
average_cvg_ops <- NULL

avg_set_size_opi <- NULL
avg_set_size_ops <- NULL

for(alpha in 1:length(alphas))
{
#  average_cvg_cdf[alpha] <- mean(marginal_coverage_cdf[alpha,])
#  average_cvg_lac[alpha] <- mean(marginal_coverage_lac[alpha,])
  average_cvg_opi[alpha] <- mean(marginal_coverage_opi[alpha,])
  average_cvg_ops[alpha] <- mean(marginal_coverage_ops[alpha,])
  
#  avg_set_size_cdf[alpha] <- mean(set_size_cdf[alpha,])
#  avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

agg_lac_cvg <- c(0.9936963407293489,
                 0.9520577335107429,
                 0.898688089436437,
                 0.8505620678624405)
agg_aps_cvg <- c(0.9924839057548822,
                 0.9522103858879243,
                 0.8991735767097603,
                 0.8519249218571119)
agg_cdf_cvg <- c(0.9943379531149088,
                 0.9527980002038952,
                 0.9008392450819224,
                 0.8542435116016236)

agg_lac_set <- c(2.386364065803044,
                 1.5278388058476216,
                 1.2699259168792847,
                 1.1208557153971348)

agg_cdf_set <- c(3.8860329796245088,
                 3.662968752115455,
                 3.4412619567857057,
                 2.996979610007956)

agg_aps_set <- c(2.1759808121890596,
                 1.5296141501919995,
                 1.2750213481039598,
                 1.128638955645595)

hist(p.values[,1],xlab = "p-value",main = expression("Histogram of P"[1]))
hist(p.values[,2],xlab = "p-value",main = expression("Histogram of P"[2]))
hist(p.values[,3],xlab = "p-value",main = expression("Histogram of P"[3]))
hist(p.values[,4],xlab = "p-value",main = expression("Histogram of P"[4]))

results_average_cvg <- data.frame(cbind(alphas,agg_lac_cvg,agg_aps_cvg,agg_cdf_cvg,average_cvg_opi,average_cvg_ops)) 
results_average_set_size <- data.frame(cbind(alphas, agg_lac_set, agg_cdf_set, agg_aps_set, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = agg_cdf_cvg,color = "Ordinal CDF")) +
  geom_line(aes(y = agg_lac_cvg,color = "Naive LAC")) +
  geom_line(aes(y = agg_aps_cvg, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

ggplt_1

# fig <- plot_ly(results_average_cvg, x = ~alphas, y = ~agg_cdf_cvg, name = "<b>Ordinal CDF", type = 'scatter', mode = 'lines') 
# fig <- fig %>% add_trace(y = ~agg_lac_cvg, name = "<b>Naive LAC", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~agg_aps_cvg, name = "<b>APS", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig <- fig %>% add_trace(y = ~average_cvg_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig %>% layout(yaxis = list(title = "marginal coverage"))

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = agg_cdf_set,color = "Ordinal CDF")) +
  geom_line(aes(y = agg_lac_set,color = "Naive LAC")) +
  geom_line(aes(y = agg_aps_set, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_1/ggplt_2 

# fig_1 <- plot_ly(results_average_set_size, x = ~alphas) 
# fig_1 <- fig_1 %>% add_lines(y = ~agg_cdf_set, name = "<b>Ordinal CDF", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~agg_lac_set, name = "<b>Naive LAC", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~agg_aps_set, name = "<b>APS", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_opi, name = "<b>Ordinal Prediction Intervals", connectgaps = TRUE)
# fig_1 <- fig_1 %>% add_lines(y = ~avg_set_size_ops, name = "<b>Ordinal Prediction Sets", connectgaps = TRUE)
# fig_1 %>% layout(yaxis = list(title = "average set size"))


#////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#=============================REAL DATA ANALYSIS===============================#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////////#

setwd("/Users/subhrasishchakraborty/Desktop/dataset")

traffic_accident_leeds <- read.csv("Traffic accidents_2019_Leeds.csv",header = TRUE)
training_pct <- 0.25
training_index <- sample(1:nrow(traffic_accident_leeds),training_pct*nrow(traffic_accident_leeds))
training_data <- traffic_accident_leeds[training_index,]
test_data <- traffic_accident_leeds[-training.index,]

model_accident <- multinom(Casualty.Severity~.,data = training_data)
posterior.probability <- predict(model_accident,newdata = test_data, "probs")
df_probs_1 <- data.frame(posterior.probability,test_data$Casualty.Severity)
names(df_probs_1)[4] <- "labels"

cal_pct <- 0.4

cal_indexes <- sample(1:nrow(df_probs_1),cal_pct*nrow(df_probs_1))
cal_labels <- df_probs_1$labels[cal_indexes]
cal_probs <- df_probs_1[cal_indexes,-4]

val_df <- df_probs_1[-cal_indexes,]
val_labels <- val_df$labels
val_probs <- val_df[,-4]

p.values <- matrix(0,nrow = nrow(val_probs),ncol = 3)

ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)

for(i in 1:nrow(val_probs))
{
  for(j in 1:3)
  {
    p.values[i,j] <- (sum(val_probs[i,j] <= cal_probs[,j]) + 1)/(length(cal_probs[,j]) + 1)
  }
}
par(mfrow=c(3,1))
hist(p.values[,1],xlab = "p-value",main = expression("Histogram of P"[1]))
hist(p.values[,2],xlab = "p-value",main = expression("Histogram of P"[2]))
hist(p.values[,3],xlab = "p-value",main = expression("Histogram of P"[3]))
#hist(p.values[,4],xlab = "p-value",main = expression("Histogram of P"[4]))

write.csv(data.frame(cal_probs,cal_labels),"/Users/subhrasishchakraborty/Desktop/calibration_data_acc_leeds.csv")
write.csv(data.frame(val_probs,val_labels),"/Users/subhrasishchakraborty/Desktop/validation_data_acc_leeds.csv")

alphas <- seq(0.05,0.15,0.01)

num_of_classes <- 3

marginal_coverage_cdf <- matrix(0, nrow = length(alphas), ncol = 5)
marginal_coverage_lac <- matrix(0, nrow = length(alphas), ncol = 5)
# marginal_coverage_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = 5)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = 5)

set_size_cdf <- matrix(0, nrow = length(alphas), ncol = 5)
set_size_lac <- matrix(0, nrow = length(alphas), ncol = 5)
# set_size_aps <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_opi <- matrix(0, nrow = length(alphas), ncol = 5)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = 5)

for(alpha in 1:length(alphas))
{
  for(trial in 1:5)
  {
    # ordinal cdf
    
    Q_hat_cdf <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[alpha])
    
    Q_hat_lac <- estimated_qhat(lac_conformity_scores,alphas[alpha])
    
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
    
    # aps
    
    # aps_prediction_interval <- ordinal_aps_function(val_probs,val_labels,1-alphas[alpha])
    # aps_coverage <- 0
    # for(m in 1:nrow(aps_prediction_interval))
    # {
    #   if(aps_prediction_interval[m,1] <= val_labels[m] && val_labels[m] <= aps_prediction_interval[m,2] )
    #     aps_coverage <- aps_coverage + 1
    # }
    # marginal_coverage_aps[alpha,trial] <- aps_coverage/length(val_labels)
    # set_size_aps[alpha,trial] <- mean(aps_prediction_interval[,2] - aps_prediction_interval[,1] + 1)
    
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
      if(val_labels[m] %in% ordinal_prediction_sets[[m]])
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
  
  avg_set_size_cdf[alpha] <- mean(set_size_cdf[alpha,])
  avg_set_size_lac[alpha] <- mean(set_size_lac[alpha,])
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}

average_cvg_aps <- c(0.9479820627802692,
                     0.9426008968609866,
                     0.9336322869955158,
                     0.9246636771300448,
                     0.9147982062780269,
                     0.8941704035874439,
                     0.8932735426008968,
                     0.8681614349775785,
                     0.863677130044843,
                     0.8573991031390135,
                     0.8484304932735427)
results_average_cvg <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

avg_set_size_aps <- c(1.6188340807174888,
                      1.5865470852017938,
                      1.5461883408071748,
                      1.495964125560538,
                      1.4385650224215247,
                      1.3551569506726457,
                      1.3300448430493275,
                      1.240358744394619,
                      1.209865470852018,
                      1.1829596412556054,
                      1.1497757847533632)

results_average_set_size <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_1 / ggplt_2
