install.packages("mvnfast")
install.packages("ggplot2")
install.packages("mvtnorm")
install.packages("patchwork")
install.packages("plotly")

library(mvnfast)
library(ggplot2)
library(mvtnorm)
library(patchwork)
library(plotly)

library(MASS)

set.seed(1234)
install.packages("HTLR")
library(HTLR)

N <- 1000 # number of observations
d <- 10 # dimensions 5, 10
K <- 8 # number of classes 4, 8
df <- gendata_MLR(N,d,K)
#qplot(df$y, bins = 6, xlab = "labels")
#corrplot::corrplot(cor(df$X))

simulated.data <- data.frame(df$X,df$y)
#View(simulated.data)

names(simulated.data)[d+1] <- "y"

#require(nnet)

nn.model <- multinom(y~.,data = simulated.data)
posterior.probs <- predict(nn.model,newdata = simulated.data, "probs")
df_probs <- data.frame(posterior.probs,simulated.data$y)

cal_pct <- 0.05

cal_index <- sample(1:N,cal_pct*N)
cal_df <- df_probs[cal_index,]
cal_labels <- cal_df$simulated.data.y
cal_probs <- cal_df[,-(K+1)]

val_df <- df_probs[-cal_index,]
val_labels <- val_df$simulated.data.y
val_probs <- val_df[,-(K+1)]

#----------------------------ORDINAL CDF----------------------------#

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

#----------------------------Naive LAC-------------------------------#

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

#----------------------------Ordinal APS-----------------------------#

ordinal_aps_function <- function(validation.scores,validation.labels,lambda)
{
  # returns a prediction interval
  starting.point <- NULL
  ending.point <- NULL
  
  for(i in 1:nrow(validation.scores))
  {
    maximum.position <- which.is.max(validation.scores[i,])
    begin <- max(maximum.position - 1, 1)
    end <- min(maximum.position + 1,K)
    if(begin < end)
    {
      q <- 0
      # starting.point[i] <- begin
      # ending.point[i] <- begin
      while(q > lambda)
      {
        q = sum(validation.scores[i,begin:end])
        end <- end - 1
      }
      starting.point[i] <- begin
      ending.point[i] <- min(end + 1,K)
    }
    else{
      starting.point[i] <- maximum.position
      ending.point[i] <- maximum.position
    }
  }
  return(cbind(starting.point,ending.point))
  
}

alphas <- seq(0.01,0.3,0.01)

ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)

p.values <- matrix(0,nrow = nrow(val_probs),ncol = K)

for(i in 1:nrow(val_probs))
{
  for(j in 1:K)
  {
    p.values[i,j] <- (sum(val_probs[i,j] >= cal_probs[,j]) + 1)/(nrow(cal_df) + 1)
  }
}

marginal.results <- data.frame(cbind(alphas,rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas))))
names(marginal.results)[2] <- "OAPS"
names(marginal.results)[3] <- "Ordinal CDF"
names(marginal.results)[4] <- "Naive LAC"
names(marginal.results)[5] <- "OPI" # Ordinal prediction intervals
names(marginal.results)[6] <- "OPS" # Prediction sets

set.size_results <- marginal.results <- data.frame(cbind(alphas,rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas)),rep(0,length(alphas))))
names(set.size_results)[2] <- "OAPS"
names(set.size_results)[3] <- "Ordinal CDF"
names(set.size_results)[4] <- "Naive LAC"
names(set.size_results)[5] <- "OPI" # Ordinal prediction intervals
names(set.size_results)[6] <- "OPS" # Prediction sets

estimated_qhat <- function(conformity_scores,a)
  {
   q <- floor((cal_pct*N + 1)*(1 - a))/(cal_pct*N + 1)
   qhat <- quantile(conformity_scores, probs = c(0.25,q))
   return(qhat[2])
  }

#---------------------------ORDINAL CDF-----------------------------#

ordinal_CDF_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  prediction.interval <- matrix(0,nrow = nrow(validation.scores), ncol = K)
  
  lower.edge <- NULL
  upper.edge <- NULL
  
  for(j in 1:nrow(validation.scores))
  {
    pos <- which.is.max(validation.scores[j,])
    for(k in 1:K)
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

#------------------------------Naive LAC----------------------------#

naive_lac_prediction_sets <- function(validation.scores,validation.labels,qhat)
{
  lac.prediction.set <- matrix(0,nrow = nrow(validation.scores), ncol = K)
  for(j in 1:nrow(validation.scores))
  {
    for(k in 1:K)
    {
      if(sum(validation.scores[j,1:k]) >= 1 - qhat)
        lac.prediction.set[j,k] <- k
      else
        lac.prediction.set[j,k] <- NA
    }
  }
  return(lac.prediction.set)
}

#----------------------------OPI intervals--------------------------#

opi_prediction_intervals <- function(p_values,level.of.significance)
{
  lower_edge <- NULL
  upper_edge <- NULL
  
  for(i in 1:nrow(p_values))
  {
    for(k in 1:K)
    {
      if(p_values[i,k] > level.of.significance)
      {
        lower_edge[i] <- k
        break
      }
    }
    for(k in K:1)
    {
      if(p_values[i,k] > level.of.significance)
      {
        upper_edge[i] <- k
        break
      }
    }
  }
  return(cbind(lower_edge,upper_edge))
}

#------------------------Prediction sets----------------------------#

opi_prediction_sets <- function(p_values,level.of.significance)
{
  opi.prediction.set <- matrix(0,nrow = nrow(p_values), ncol = K)
  for(j in 1:nrow(p_values))
  {
    for(k in 1:K)
    {
      if(p_values[j,k] > level.of.significance)
        opi.prediction.set[j,k] <- k
      else
        opi.prediction.set[j,k] <- NA
    }
  }
  return(opi.prediction.set)
}

#------------------------------------------------------------------#

evaluate_set_size <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:nrow(prediction_set))
      set_size[i] <- sum(!is.na(prediction_set[i,]))
  return(mean(set_size))
}


for(l in 1:length(alphas))
{
  
  #OAPS
  
  pred_tmp_1 <- ordinal_aps_function(val_probs,val_labels,1-alphas[l])
  count_1 <- 0
  for(m in 1:nrow(pred_tmp_1))
  {
    if(pred_tmp_1[m,1] <= val_labels[m] && val_labels[m] <= pred_tmp_1[m,2] )
      count_1 <- count_1 + 1
  }
  marginal.results[l,2] <- count_1/length(val_labels)
  set.size_results[l,2] <- mean(pred_tmp_1[,2] - pred_tmp_1[,1] + 1)
  
  # Ordinal CDF
  
  Q_hat <- estimated_qhat(ordinal_cdf_conformity_scores,alphas[l])
  pred_tmp_2 <- ordinal_CDF_prediction_sets(val_probs,val_labels,Q_hat)
  count_2 <- 0
  for(m in 1:nrow(pred_tmp_1))
  {
    if(pred_tmp_2[m,1] <= val_labels[m] && val_labels[m] <= pred_tmp_2[m,2] )
      count_2 <- count_2 + 1
  }
  marginal.results[l,3] <- count_2/length(val_labels)
  set.size_results[l,3] <- mean(pred_tmp_2[,2] - pred_tmp_2[,1] + 1)
  
  # Naive LAC
  
  pred_tmp_3 <- naive_lac_prediction_sets(val_probs,val_df,Q_hat)
  count_3 <- 0
  for(m in 1:nrow(pred_tmp_3))
  {
    if(val_labels[m] %in% pred_tmp_3[m,])
      count_3 <- count_3 + 1
  }
  marginal.results[l,4] <- count_3/length(val_labels)
  set.size_results[l,4] <- evaluate_set_size(pred_tmp_3)
  # nrow(pred_tmp_3)
  
  # OPI intervals
  
  pred_tmp_4 <- opi_prediction_intervals(p.values,alphas[l])
  count_4 <- 0
  for(m in 1:nrow(pred_tmp_4))
  {
    if(pred_tmp_4[m,1] <= val_labels[m] && val_labels[m] <= pred_tmp_4[m,2] )
      count_4 <- count_4 + 1
  }
  marginal.results[l,5] <- count_4/length(val_labels)  
  set.size_results[l,5] <- mean(pred_tmp_4[,2] - pred_tmp_4[,1] + 1)
  
  # OPS
  pred_tmp_5 <- opi_prediction_sets(p.values,alphas[l])
  count_5 <- 0
  for(m in 1:nrow(pred_tmp_5))
  {
    if(val_labels[m] %in% pred_tmp_5[m,])
      count_5 <- count_5 + 1
  }
  marginal.results[l,6] <- count_5/length(val_labels)
  set.size_results[l,6] <- evaluate_set_size(pred_tmp_5)
}

for(i in 1:nrow(set.size_results))
{
  for(j in 1:ncol(set.size_results))
  {
    if(set.size_results[i,j] == -Inf)
      set.size_results[i,j] <- 1
  }
}

res1 <- data.frame(cbind(rep(marginal.results$alphas,5),c(marginal.results[,2],marginal.results[,3],marginal.results[,4],marginal.results[,5],marginal.results[,6])))
names(res1)[1] <- "alpha"
names(res1)[2] <- "marginal.coverage"
res1 <- cbind(res1,c(rep("OAPS",length(alphas)),rep("Ordinal CDF",length(alphas)),rep("Naive LAC",length(alphas)),rep("OPI",length(alphas)),rep("OPS",length(alphas))))
names(res1)[3] <- "method"

res1

ggplt_1 <- ggplot(data = res1, aes(x = alpha)) + geom_line(aes(y = marginal.coverage, colour = method))
ggplt_1 + ylab("marginal coverage")

res2 <- data.frame(cbind(rep(set.size_results$alphas,5),c(set.size_results[,2],set.size_results[,3],set.size_results[,4],set.size_results[,5],set.size_results[,6])))
names(res2)[1] <- "alpha"
names(res2)[2] <- "set.size"
res2 <- cbind(res2,c(rep("OAPS",length(alphas)),rep("Ordinal CDF",length(alphas)),rep("Naive LAC",length(alphas)),rep("OPI",length(alphas)),rep("OPS",length(alphas))))
names(res2)[3] <- "method"

res2

ggplt_2 <- ggplot(data = res2, aes(x = alpha)) + geom_line(aes(y = set.size, colour = method))
ggplt_2 + ylab("set size")

ggplt_1/ggplt_2
