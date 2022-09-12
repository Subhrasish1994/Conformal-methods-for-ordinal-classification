# Implementation of OAPS, naive LAC, ordinal CDF comparing the average coverage and set size with proposed formulation

install.packages("mvtnorm")
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

set.seed(12)

alpha <- 0.01
N <- 1000
labels.generated <- rep(0,N)
for(i in 1:N)
{
  labels.generated[i] <- sample(1:4,prob = rep(1/4,4))
}
df <- rep(0,N)
for(i in 1:N)
{
  if(labels.generated[i] == 1)
    df[i] <- rnorm(1,1,2)
  else if(labels.generated[i] == 2)
    df[i] <- rnorm(1,3,2)
  else if(labels.generated[i] == 3)
    df[i] <- rnorm(1,5,1)
  else
    df[i] <- rnorm(1,0,1)
}
df <- data.frame(df,labels.generated)
names(df)[2] <- "class"

require(nnet)

nn.model <- multinom(class~.,data = df)
posterior.probs <- predict(nn.model,newdata = df, "probs")
df_probs <- data.frame(posterior.probs,labels.generated)

# View(df_probs)

cal_pct <- 0.05
cal_index <- sample(1:N,cal_pct*N)
cal_df <- df_probs[cal_index,]
cal_labels <- cal_df$labels.generated
cal_probs <- cal_df[,-5]

val_df <- df_probs[-cal_index,]
val_labels <- val_df$labels.generated
val_probs <- val_df[,-5]

ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)

q <- floor((cal_pct*N + 1)*(1 - alpha))/(cal_pct*N + 1)

qhat_ordinal_cdf <- quantile(ordinal_cdf_conformity_scores, probs = c(0.25,q))
qhat_ordinal_cdf

qhat_lac <- quantile(lac_conformity_scores,probs = c(0.25,q))
qhat_lac

p1 <- rep(0,nrow(val_df))
p2 <- rep(0,nrow(val_df))
p3 <- rep(0,nrow(val_df))
p4 <- rep(0,nrow(val_df))

for(j in 1:nrow(val_df))
{
  p1[j] <- (sum(val_probs[j,1] >= cal_df[,1]) + 1)/(nrow(cal_df) + 1)
  p2[j] <- (sum(val_probs[j,2] >= cal_df[,2]) + 1)/(nrow(cal_df) + 1)
  p3[j] <- (sum(val_probs[j,3] >= cal_df[,3]) + 1)/(nrow(cal_df) + 1)
  p4[j] <- (sum(val_probs[j,4] >= cal_df[,4]) + 1)/(nrow(cal_df) + 1)
}

prediction.set.lac <- matrix(0, ncol = 4, nrow = nrow(val_df))
prediction.set.ordinal.cdf <- matrix(0, ncol = 4, nrow = nrow(val_df))
prediction.set <- matrix(0, ncol = 4, nrow = nrow(val_df))

for(j in 1:nrow(val_df))
{
  pos_max <- which.is.max(val_probs[j,])
  for(k in 1:4)
  {
    if(sum(val_probs[j,1:k]) >= val_probs[j,pos_max] - qhat_ordinal_cdf[2] && sum(val_probs[j,1:k]) <= val_probs[j,pos_max] + qhat_ordinal_cdf[2])
      prediction.set.ordinal.cdf[j,k] <- k
    else
      prediction.set.ordinal.cdf[j,k] <- NA
  }
}

# View(prediction.set.ordinal.cdf)

for(j in 1:nrow(val_df))
{
  for(k in 1:4)
  {
    if(sum(val_probs[j,1:k]) >= 1 - qhat_lac[2])
      prediction.set.lac[j,k] <- k
    else
      prediction.set.lac[j,k] <- NA
  }
}

# View(prediction.set.lac)

for(j in 1:nrow(val_probs))
{
  if(p1[j] > alpha)
    prediction.set[j,1] <- 1
  else
    prediction.set[j,1] <- NA
  
  if(p2[j] > alpha)
    prediction.set[j,2] <- 2
  else
    prediction.set[j,2] <- NA

  if(p3[j] > alpha)
    prediction.set[j,3] <- 3
  else
    prediction.set[j,3] <- NA
  
  if(p4[j] > alpha)
    prediction.set[j,4] <- 4
  else
    prediction.set[j,4] <- NA
}

# View(prediction.set)

prediction.interval = matrix(0, nrow = nrow(val_probs), ncol = 2)
probs <- data.frame(cbind(p1,p2,p3,p4))

for(j in 1:nrow(val_probs))
{
  for(k in 1:4)
  {
    if(probs[j,k] > alpha)
    {
      prediction.interval[j,1] <- k
      break
    }
  }
  
  for(k1 in 4:1)
  {
    if(probs[j,k1] > alpha)
    {
      prediction.interval[j,2] <- k1
      break
    }
  }
}

# View(prediction.interval)

count1 <- 0
count2 <- 0
count3 <- 0 
count4 <- 0

for(i in 1:length(val_labels))
{
  if(val_labels[i] %in% prediction.set.ordinal.cdf[i,])
    count1 <- count1 + 1
  if(val_labels[i] %in% prediction.set.lac[i,])
    count2 <- count2 + 1
  if(val_labels[i] %in% prediction.set[i,])
    count3 <- count3 + 1
  if(val_labels[i] >= prediction.interval[i,1] && val_labels[i] <= prediction.interval[i,2])
    count4 <- count4 + 1
}
 count1/length(val_labels)
 count2/length(val_labels)
 count3/length(val_labels)
 count4/length(val_labels)
 
s1 <- NULL
s2 <- NULL
s3 <- NULL
s4 <- NULL

for(i in 1:length(val_labels))
{
  s1[i] <- sum(!is.na(prediction.set.ordinal.cdf[i,]))
  s2[i] <- sum(!is.na(prediction.set.lac[i,]))
  s3[i] <- sum(!is.na(prediction.set[i,]))
  s4[i] <- prediction.interval[i,2] - prediction.interval[i,1]
 }
mean(s1)
mean(s2)
mean(s3)
mean(s4)

#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#


set.seed(12)

alpha <- 0.01
d <- 100 # dimension should be varied from 2, 3, 6, 9, 15, 20, 30, 40, 50
N <- 1000
labels.generated <- rep(0,N)
for(i in 1:N)
{
  labels.generated[i] <- sample(1:4,prob = rep(1/4,4))
}
df <- matrix(0, nrow = N, ncol = d)
for(i in 1:N)
{ 
  if(labels.generated[i] == 1)
    df[i,] <- rmvnorm(1,c(-1,rep(1,d-1)),diag(d))
  else if(labels.generated[i] == 2)
    df[i,] <- rmvnorm(1,0.6*c(-1,rep(1,d-1)),diag(d))
  else if(labels.generated[i] == 3)
    df[i,] <- rnorm(1,-1*c(-1,rep(1,d-1)),diag(d))
  else
    df[i,] <- rnorm(1,2*c(-1,rep(1,d-1)),diag(d))
}
#View(df)
df <- data.frame(df,labels.generated)
names(df)[d+1] <- "class"

nn.model <- multinom(class~.,data = df)
posterior.probs <- predict(nn.model,newdata = df, "probs")
df_probs <- data.frame(posterior.probs,labels.generated)

# View(df_probs)

cal_pct <- 0.05
cal_index <- sample(1:N,cal_pct*N)
cal_df <- df_probs[cal_index,]
cal_labels <- cal_df$labels.generated
cal_probs <- cal_df[,-5]

val_df <- df_probs[-cal_index,]
val_labels <- val_df$labels.generated
val_probs <- val_df[,-5]

ordinal_cdf_conformity_scores <- ordinal_CDF_score_function(cal_probs,cal_labels)
lac_conformity_scores <- lac_score_function(cal_probs,cal_labels)

q <- floor((cal_pct*N + 1)*(1 - alpha))/(cal_pct*N + 1)

qhat_ordinal_cdf <- quantile(ordinal_cdf_conformity_scores, probs = c(0.25,q))
qhat_ordinal_cdf

qhat_lac <- quantile(lac_conformity_scores,probs = c(0.25,q))
qhat_lac

p1 <- rep(0,nrow(val_df))
p2 <- rep(0,nrow(val_df))
p3 <- rep(0,nrow(val_df))
p4 <- rep(0,nrow(val_df))

for(j in 1:nrow(val_df))
{
  p1[j] <- (sum(val_probs[j,1] >= cal_df[,1]) + 1)/(nrow(cal_df) + 1)
  p2[j] <- (sum(val_probs[j,2] >= cal_df[,2]) + 1)/(nrow(cal_df) + 1)
  p3[j] <- (sum(val_probs[j,3] >= cal_df[,3]) + 1)/(nrow(cal_df) + 1)
  p4[j] <- (sum(val_probs[j,4] >= cal_df[,4]) + 1)/(nrow(cal_df) + 1)
}

prediction.set.lac <- matrix(0, ncol = 4, nrow = nrow(val_df))
prediction.set.ordinal.cdf <- matrix(0, ncol = 4, nrow = nrow(val_df))
prediction.set <- matrix(0, ncol = 4, nrow = nrow(val_df))

for(j in 1:nrow(val_df))
{
  pos_max <- which.is.max(val_probs[j,])
  for(k in 1:4)
  {
    if(sum(val_probs[j,1:k]) >= val_probs[j,pos_max] - qhat_ordinal_cdf[2] && sum(val_probs[j,1:k]) <= val_probs[j,pos_max] + qhat_ordinal_cdf[2])
      prediction.set.ordinal.cdf[j,k] <- k
    else
      prediction.set.ordinal.cdf[j,k] <- NA
  }
}

# View(prediction.set.ordinal.cdf)

for(j in 1:nrow(val_df))
{
  for(k in 1:4)
  {
    if(sum(val_probs[j,1:k]) >= 1 - qhat_lac[2])
      prediction.set.lac[j,k] <- k
    else
      prediction.set.lac[j,k] <- NA
  }
}

# View(prediction.set.lac)

for(j in 1:nrow(val_probs))
{
  if(p1[j] > alpha)
    prediction.set[j,1] <- 1
  else
    prediction.set[j,1] <- NA
  
  if(p2[j] > alpha)
    prediction.set[j,2] <- 2
  else
    prediction.set[j,2] <- NA
  
  if(p3[j] > alpha)
    prediction.set[j,3] <- 3
  else
    prediction.set[j,3] <- NA
  
  if(p4[j] > alpha)
    prediction.set[j,4] <- 4
  else
    prediction.set[j,4] <- NA
}

# View(prediction.set)

prediction.interval = matrix(0, nrow = nrow(val_probs), ncol = 2)
probs <- data.frame(cbind(p1,p2,p3,p4))

for(j in 1:nrow(val_probs))
{
  for(k in 1:4)
  {
    if(probs[j,k] > alpha)
    {
      prediction.interval[j,1] <- k
      break
    }
  }
  
  for(k1 in 4:1)
  {
    if(probs[j,k1] > alpha)
    {
      prediction.interval[j,2] <- k1
      break
    }
  }
}

# View(prediction.interval)

count1 <- 0
count2 <- 0
count3 <- 0 
count4 <- 0

for(i in 1:length(val_labels))
{
  if(val_labels[i] %in% prediction.set.ordinal.cdf[i,])
    count1 <- count1 + 1
  if(val_labels[i] %in% prediction.set.lac[i,])
    count2 <- count2 + 1
  if(val_labels[i] %in% prediction.set[i,])
    count3 <- count3 + 1
  if(val_labels[i] >= prediction.interval[i,1] && val_labels[i] <= prediction.interval[i,2])
    count4 <- count4 + 1
}
count1/length(val_labels)
count2/length(val_labels)
count3/length(val_labels)
count4/length(val_labels)

s1 <- NULL
s2 <- NULL
s3 <- NULL
s4 <- NULL

for(i in 1:length(val_labels))
{
  s1[i] <- sum(!is.na(prediction.set.ordinal.cdf[i,]))
  s2[i] <- sum(!is.na(prediction.set.lac[i,]))
  s3[i] <- sum(!is.na(prediction.set[i,]))
  s4[i] <- prediction.interval[i,2] - prediction.interval[i,1]
}
mean(s1)
mean(s2)
mean(s3)
mean(s4)

