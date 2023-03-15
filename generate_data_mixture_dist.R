ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("mvtnorm","ggplot2","caret","ipred","devtools","sm","hrbrthemes","tidyr","dplyr","nnet","viridis","scatterplot3d", "patchwork")
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
    df[index,] <- 0.2*mvrnorm(1,c(-1,0),sigma)+0.8*mvrnorm(1,c(-1,-1),sigma)
  else if(labels.generated[index] == 2)
    df[index,] <- 0.2*mvrnorm(1,c(-1,-1),sigma)+0.8*mvrnorm(1,c(0,-1),sigma)
  else if(labels.generated[index] == 3)
    df[index,] <- 0.2*mvrnorm(1,c(0,-1),sigma)+0.8*mvrnorm(1,c(1,-1),sigma)
  else
    df[index,] <- 0.2*mvrnorm(1,c(1,-1),sigma)+0.8*mvrnorm(1,c(-1,0),sigma)
}

df <- data.frame(cbind(df,labels.generated))
names(df)[d+1] <- "label"
# View(df)

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
                     0.9897660818713451,
                     0.9853801169590642,
                     0.9678362573099415,
                     0.955263157894737,
                     0.9552631578947368,
                     0.9394736842105263,
                     0.9251461988304092,
                     0.9225146198830408,
                     0.9058479532163742,
                     0.9032163742690058,
                     0.8964912280701753,
                     0.8801169590643274,
                     0.8807017543859649,
                     0.8640350877192983,
                     0.8538011695906433,
                     0.854970760233918,
                     0.8210526315789475,
                     0.8160818713450292,
                     0.7929824561403508)
avg_set_size_aps <- c(3.9859649122807013,
                      3.6131578947368426,
                      3.4985380116959064,
                      3.3058479532163743,
                      3.2219298245614034,
                      3.2026315789473685,
                      3.1108187134502923,
                      3.0514619883040934,
                      3.022514619883041,
                      2.962573099415205,
                      2.9257309941520466,
                      2.910818713450292,
                      2.817836257309941,
                      2.826315789473685,
                      2.739181286549708,
                      2.6631578947368424,
                      2.651169590643275,
                      2.4997076023391815,
                      2.4713450292397656,
                      2.3929824561403508)

results_average_cvg_mixture <- data.frame(cbind(alphas,average_cvg_aps,average_cvg_cdf,average_cvg_lac,average_cvg_opi,average_cvg_ops))
results_average_set_size_mixture <- data.frame(cbind(alphas, avg_set_size_aps,avg_set_size_cdf, avg_set_size_lac, avg_set_size_opi, avg_set_size_ops))

ggplt_1 <- ggplot(data = results_average_cvg_mixture, aes(x = alphas))+
  geom_line(aes(y = average_cvg_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = average_cvg_lac,color = "Naive LAC")) +
  geom_line(aes(y = average_cvg_aps, color = "Ordinal APS")) +
  geom_line(aes(y = average_cvg_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = average_cvg_ops,color = "Ordinal Prediction Sets")) +
  labs(x = "alphas", y = "Marginal coverage", color = "Legend")



ggplt_2 <- ggplot(data = results_average_set_size_mixture,aes(x = alphas)) +
  geom_line(aes(y = avg_set_size_cdf,color = "Ordinal CDF")) +
  geom_line(aes(y = avg_set_size_lac,color = "Naive LAC")) +
  geom_line(aes(y = avg_set_size_aps, color = "Ordinal APS")) +
  geom_line(aes(y = avg_set_size_opi,color = "Ordinal Prediction Interval")) +
  geom_line(aes(y = avg_set_size_ops,color = "Ordinal Prediction Sets")) + ylim(0,4) +
  labs(x = "alphas", y = "Set size", color = "Legend")

plt_multinom <- ggplt_1 / ggplt_2 

ggarrange(ggplt_1,ggplt_2,ncol = 2, nrow = 1,legend = "bottom")

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
                 0.9977883675464321,
                 0.9936382627854938,
                 0.9654044686349561,
                 0.9508122918568036,
                 0.9531242347327493,
                 0.9258135332682766,
                 0.8996399947094315,
                 0.8945678829398229,
                 0.875515140210533,
                 0.8710779248171108,
                 0.8632774916833575,
                 0.8498310044704713,
                 0.8343308785128022,
                 0.8081876192143669,
                 0.8006293671824135,
                 0.8034072356722479,
                 0.7654318475410117,
                 0.7690216541633867,
                 0.7471316544668876)
aps_mil_cvg <- c(1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 1.0,
                 0.9989473684210527,
                 0.9943992689419938)
aps_mod_cvg <- c(1.0,
                 0.9818240333904349,
                 0.968607985944692,
                 0.9539413616291551,
                 0.9450731614269895,
                 0.9441124241413583,
                 0.9342393948758249,
                 0.9187881242828324,
                 0.9126106498381661,
                 0.9031929243534098,
                 0.9032117727991895,
                 0.9017656140763108,
                 0.8723237603973999,
                 0.8905557833944464,
                 0.8696474432067601,
                 0.8607197642714108,
                 0.867383026318097,
                 0.8419375735513531,
                 0.8302017951669587,
                 0.8122116579785658)
aps_sev_cvg <- c(1.0,
                 0.9767624687025631,
                 0.9743363728240189,
                 0.9475490831200112,
                 0.9201321502773944,
                 0.9170628923993481,
                 0.8912541473719351,
                 0.8731465995382959,
                 0.8735522361454408,
                 0.8367742284772595,
                 0.8297291709040178,
                 0.815481700072367,
                 0.7840131236637672,
                 0.7868124041286821,
                 0.7684024538756248,
                 0.7319658795271187,
                 0.7330397797096886,
                 0.6587511673884522,
                 0.6542688807961775,
                 0.5999509110266862)

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

conditional_coverage_mixture <- cbind(cdf_conditional_coverage,lac_conditional_coverage,
                                      aps_conditional_coverage,opi_conditional_coverage,
                                      ops_conditional_coverage)

res_2 <- read.csv("/Users/Owner/Desktop/result_2.csv")

n <- length(alphas)*5

distribution <- c(rep('gaussian',n),rep('mixture gaussian',n))
method <- rep(c(rep('Ordinal CDF',length(alphas)),rep('Ordinal LAC',length(alphas)),
                rep('Ordinal APS',length(alphas)),rep('Ordinal Prediction Interval',length(alphas)),
                rep('Ordinal Prediction Set',length(alphas))),2)
miscoverage <- rep(rep(alphas,5),2)

marginal.coverage <- c(results_average_cvg_gaussian[,3],results_average_cvg_gaussian[,4],
                       results_average_cvg_gaussian[,2],results_average_cvg_gaussian[,5],
                       results_average_cvg_gaussian[,6],results_average_cvg_mixture[,3],results_average_cvg_mixture[,4],
                       results_average_cvg_mixture[,2],results_average_cvg_mixture[,5],
                       results_average_cvg_mixture[,6])
set.size <- c(results_average_set_size_gaussian[,3],results_average_set_size_gaussian[,4],
              results_average_set_size_gaussian[,2],results_average_set_size_gaussian[,5],
              results_average_set_size_gaussian[,6],results_average_set_size_mixture[,3],
              results_average_set_size_mixture[,4],results_average_set_size_mixture[,2],
              results_average_set_size_mixture[,5],results_average_set_size_mixture[,6])
conditional.coverage <- c(conditional_coverage_gaussian[,1],conditional_coverage_gaussian[,2],
                          conditional_coverage_gaussian[,3],conditional_coverage_gaussian[,4],
                          conditional_coverage_gaussian[,5],conditional_coverage_mixture[,1],
                          conditional_coverage_mixture[,2],conditional_coverage_mixture[,3],
                          conditional_coverage_mixture[,4],conditional_coverage_mixture[,5])
length(set.size)
res_3 <- data.frame(cbind(distribution,method,as.numeric(miscoverage),as.numeric(marginal.coverage),as.numeric(set.size),as.numeric(conditional.coverage)))
write.csv(res_3,"/Users/Owner/Desktop/res.csv")
names(res_3)[3] <- "miscoverage"
names(res_3)[4] <- "marginal.coverage"
names(res_3)[5] <- "set.size"
names(res_3)[6] <- "conditional.coverage"
View(res_3)

res_3 <- read.csv("/Users/Owner/Desktop/res.csv",header = TRUE)

# cols <- c("Ordinal APS" = "red",
#           "Ordinal CDF" = "green",
#           "Ordinal LAC" = "blue",
#           "Ordinal Prediction Interval" = "darkgreen",
#           "Ordinal Prediction Set" = "orange")

plt_1 <- ggplot(res_3,aes(x=miscoverage))+geom_line(aes(y = marginal.coverage,color=method))
plt_1 <- plt_1 + facet_grid(.~distribution) + ylab("marginal coverage") 

plt_2 <- ggplot(res_3,aes(x=miscoverage))+geom_line(aes(y = set.size,color=method))
plt_2 <- plt_2 + facet_grid(.~distribution) + ylab("set size") 

plt_3 <- ggplot(res_3,aes(x=miscoverage))+geom_line(aes(y = conditional.coverage,color=method))
plt_3 <- plt_3 + facet_grid(.~distribution) + ylab("conditional coverage") 

ggarrange(plt_1,plt_2,plt_3,ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
