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

# 1. Implement SVM, QDA instead of logistic regression

# 2. Compare the marginal coverage and set sizes

#---------------------------Traffic accident data UK---------------------------#

setwd("/Users/subhrasishchakraborty/Desktop/dataset")

traffic_accident_leeds <- read.csv("Traffic accidents_2019_Leeds.csv",header = TRUE)

model_accident <- multinom(Casualty.Severity~.,data = traffic_accident_leeds)
posterior.probability <- predict(model_accident,newdata = traffic_accident_leeds, "probs")
df_probs_1 <- data.frame(posterior.probability,traffic_accident_leeds$Casualty.Severity)
names(df_probs_1)[length(unique(traffic_accident_leeds$Casualty.Severity))+1] <- "labels"

cal_pct <- 0.4

N <- nrow(traffic_accident_leeds)

alphas <- seq(0.01,0.15,0.01)

num_of_classes <- length(unique(traffic_accident_leeds$Casualty.Severity))

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
  
  cal_scores <- NULL
  
  for(i in 1:nrow(cal_probs))
    cal_scores[i] <- cal_probs[i,cal_labels[i]]
  
  for(i in 1:nrow(val_probs))
  {
    for(j in 1:num_of_classes)
    {
      p.values[i,j] <- (sum(val_probs[i,j] >= cal_scores) + 1)/(length(cal_scores) + 1)
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
      if(val_labels[m] %in% ordinal_prediction_sets[[m]])
        ops_coverage <- ops_coverage + 1
    }
    marginal_coverage_ops[alpha,trial] <- ops_coverage/length(val_labels)
    set_size_ops[alpha,trial] <- evaluate_set_size_ops(ordinal_prediction_sets)
  }
}

write.csv(data.frame(cbind(cal_probs,cal_labels)),"/Users/subhrasishchakraborty/Desktop/calibration_data.csv")
write.csv(data.frame(cbind(val_probs,val_labels)), "/Users/subhrasishchakraborty/Desktop/validation_data.csv")

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

average_cvg_aps <- c(0.9930131004366812,
                     0.9851528384279475,
                     0.9694323144104805,
                     0.9537117903930131,
                     0.9432314410480348,
                     0.9310043668122271,
                     0.9161572052401746,
                     0.9109170305676855,
                     0.9091703056768559,
                     0.8890829694323145,
                     0.8786026200873363,
                     0.8759825327510917,
                     0.8655021834061134,
                     0.856768558951965,
                     0.8541484716157205)
avg_set_size_aps <- c(2.0890829694323143,
                      1.9004366812227071,
                      1.7318777292576422,
                      1.6087336244541486,
                      1.5275109170305678,
                      1.4768558951965063,
                      1.36943231441048,
                      1.3301310043668122,
                      1.3187772925764192,
                      1.2558951965065501,
                      1.227947598253275,
                      1.2034934497816594,
                      1.1685589519650654,
                      1.1318777292576419,
                      1.1056768558951968)
results_average_cvg <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

# results_average_cvg <- data.frame(cbind(alphas,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
# results_average_set_size <- data.frame(cbind(alphas,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_1 / ggplt_2


#--------------------------New Zealand earthquake data-------------------------#

earthquake_nz_intensity <- read.csv("earthquakes_NZ.csv", header = TRUE)

names(earthquake_nz_intensity)
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
model_earthquake <- multinom(labels~latitude+longitude+depth,data = earthquake_nz_intensity)
posterior.probability <- predict(model_earthquake,newdata = earthquake_nz_intensity, "probs")
df_probs_1 <- data.frame(posterior.probability,earthquake_nz_intensity$labels)

names(df_probs_1)[length(unique(earthquake_nz_intensity$labels))+1] <- "labels"

N <- nrow(earthquake_nz_intensity)
num_of_classes <- length(unique(earthquake_nz_intensity$labels))

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
  
  cal_scores <- NULL
  
  for(i in 1:nrow(cal_probs))
    cal_scores[i] <- cal_probs[i,cal_labels[i]]
  
  for(i in 1:nrow(val_probs))
  {
    for(j in 1:num_of_classes)
    {
      p.values[i,j] <- (sum(val_probs[i,j] >= cal_scores) + 1)/(length(cal_scores) + 1)
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
      if(val_labels[m] %in% ordinal_prediction_sets[[m]])
        ops_coverage <- ops_coverage + 1
    }
    marginal_coverage_ops[alpha,trial] <- ops_coverage/length(val_labels)
    set_size_ops[alpha,trial] <- evaluate_set_size_ops(ordinal_prediction_sets)
  }
}

write.csv(data.frame(cbind(cal_probs,cal_labels)),"/Users/subhrasishchakraborty/Desktop/calibration_data.csv")
write.csv(data.frame(cbind(val_probs,val_labels)), "/Users/subhrasishchakraborty/Desktop/validation_data.csv")

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

average_cvg_aps <- c(0.9901525546856081,
                     0.9828073290822503,
                     0.9729598837678586,
                     0.9631931552183387,
                     0.9543950278472838,
                     0.945596900476229,
                     0.9350230042779886,
                     0.9217854548389701,
                     0.9125837436435548,
                     0.9032205989183953,
                     0.8933731536040035,
                     0.8827992574057631,
                     0.872467511502139,
                     0.8621357655985149,
                     0.8516425861651464)
avg_set_size_aps <- c(2.0890829694323143,
                      1.9004366812227071,
                      1.7318777292576422,
                      1.6087336244541486,
                      1.5275109170305678,
                      1.4768558951965063,
                      1.36943231441048,
                      1.3301310043668122,
                      1.3187772925764192,
                      1.2558951965065501,
                      1.227947598253275,
                      1.2034934497816594,
                      1.1685589519650654,
                      1.1318777292576419,
                      1.1056768558951968)
results_average_cvg <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

# results_average_cvg <- data.frame(cbind(alphas,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
# results_average_set_size <- data.frame(cbind(alphas,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_1 / ggplt_2

#-----------------------------Lumbar conformal data----------------------------#

alphas <- c(0.01, 0.05, 0.1 , 0.15, 0.2)
num_of_trials <- 100


marginal_coverage_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
marginal_coverage_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

set_size_opi <- matrix(0, nrow = length(alphas), ncol = num_of_trials)
set_size_ops <- matrix(0, nrow = length(alphas), ncol = num_of_trials)

for(alpha in 1:length(alphas))
{
  for(trial in 1:num_of_trials)
  {
    cal_indexes <- sample(1:N,cal_pct*N)
    cal_labels <- df_probs_1$labels[cal_indexes]
    cal_probs <- df_probs_1[cal_indexes,-(num_of_classes+1)]
    
    val_df <- df_probs_1[-cal_indexes,]
    val_labels <- val_df$labels
    val_probs <- val_df[,-(num_of_classes+1)]
    
    p.values <- matrix(0,nrow = nrow(val_probs),ncol = num_of_classes)
    
    # ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
    # lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)
    
    cal_scores <- NULL
    
    for(i in 1:nrow(cal_probs))
      cal_scores[i] <- cal_probs[i,cal_labels[i]]
    
    for(i in 1:nrow(val_probs))
    {
      for(j in 1:num_of_classes)
      {
        p.values[i,j] <- (sum(val_probs[i,j] >= cal_scores) + 1)/(length(cal_scores) + 1)
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
      if(val_labels[m] %in% ordinal_prediction_sets[[m]])
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
  
  avg_set_size_opi[alpha] <- mean(set_size_opi[alpha,])
  avg_set_size_ops[alpha] <- mean(set_size_ops[alpha,])
}
average_cvg_aps <- c(0.9936963407293489,
                     0.9520577335107429,
                     0.898688089436437,
                     0.8505620678624405,
                     0.8070639618855242)
average_cvg_cdf <- c(0.9924839057548822,
                     0.9522103858879243,
                     0.8991735767097603,
                     0.8519249218571119,
                     0.82071381883283)
average_cvg_lac <- c(0.9943379531149088,
                     0.9527980002038952,
                     0.9008392450819224,
                     0.8542435116016236,
                     0.8195292178479257)

avg_set_size_aps <- c(2.1759808121890596,
                      1.5296141501919995,
                      1.2750213481039598,
                      1.128638955645595,
                      1.0446423022825382)
avg_set_size_cdf <- c(3.8860329796245088,
                      3.662968752115455,
                      3.4412619567857057,
                      2.996979610007956,
                      2.0600152380719168)
avg_set_size_lac <- c(2.386364065803044,
                      1.5278388058476216,
                      1.2699259168792847,
                      1.1208557153971348,
                      1.0148694152140951)
results_average_cvg <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")

ggplt_2 <- ggplot(data = results_average_set_size,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Set size", color = "Legend")

ggplt_1 / ggplt_2

