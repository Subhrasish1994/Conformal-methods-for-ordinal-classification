# conditional coverage of lumbar conformal data

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

cal_df <- read.csv("/Users/subhrasishchakraborty/Desktop/cal_df.csv")
val_df <- read.csv("/Users/subhrasishchakraborty/Desktop/val_df.csv")

cal_scores <- cal_df[,2:5]
cal_labels <- cal_df[,1]+1

val_scores <- val_df[,2:5]
val_labels <- val_df[,1]+1

hist(c(cal_probs[,1],cal_probs[,2],cal_probs[,3],cal_probs[,4],val_probs[,1],val_probs[,2],val_probs[,3],val_probs[,4]),xlab = "Posterior probabilities",main = "Histogram of Posterior probabilities of lumbar conformal data")

alphas <- c(0.01,0.05,0.1,0.15,0.2)

cond_coverage_opi <- matrix(0,nrow = 4,ncol = length(alphas))
cond_coverage_ops <- matrix(0,nrow = 4,ncol = length(alphas))

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

evaluate_set_size_ops <- function(prediction_set)
{
  set_size <- NULL
  for(i in 1:length(prediction_set))
    set_size[i] <- length(as.vector(prediction_set[[i]]))
  return(mean(set_size))
}

for(alpha in 1:length(alphas))
{
  p.values <- matrix(0,nrow = nrow(val_scores),ncol = 4)
  
  for(i in 1:nrow(val_scores))
  {
    for(j in 1:4)
    {
      p.values[i,j] <- (sum(val_scores[i,j] >= cal_scores[which(cal_labels == j),j])+1)/(length(cal_scores[which(cal_labels == j),j]) + 1)
    }
  }
  
  ordinal_prediction_intervals <- opi_prediction_intervals(p.values,alphas[alpha])
  ordinal_prediction_sets <- opi_prediction_sets(p.values,alphas[alpha])
  
  for(class in 1:4)
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
cond_coverage_opi
cond_coverage_ops
cond_coverage_aps <- data.frame(rbind(aps_covg_neg,
                                      aps_covg_mil,
                                      aps_covg_mod,
                                      aps_covg_sev))
cond_coverage_lac <- data.frame(rbind(lac_covg_neg,
                                      lac_covg_mil,
                                      lac_covg_mod,
                                      lac_covg_sev))
cond_coverage_cdf <- data.frame(rbind(cdf_covg_neg,
                                      cdf_covg_mil,
                                      cdf_covg_mod,
                                      cdf_covg_sev))
aps_covg_neg <- c(0.9957744422239848,
                  0.9854834614948661,
                  0.9684389727755139,
                  0.9512777136963316,
                  0.9395030651278888)
aps_covg_mil <- c(0.988885034849625,
                  0.8826516777276306,
                  0.7484552326602397,
                  0.62330308257775,
                  0.5402431813432913)
aps_covg_mod <- c(0.9901317784979403,
                  0.884418944754871,
                  0.7516118474214574,
                  0.6433459244789762,
                  0.5694407470014649)
aps_covg_sev <- c(0.9633154404534454,
                  0.888072658461072,
                  0.802570032317735,
                  0.7268548418625115,
                  0.6878354013945277)
lac_covg_neg <- c(0.9966368025960292,
                  0.9850210449033124,
                  0.9675479921936464,
                  0.950506048775898,
                  0.935051017779028)
cdf_covg_neg <- c(0.99637552117294,
                  0.9685295335605324,
                  0.9396192268743654,
                  0.9335101717797961,
                  0.9330596298314933)

lac_covg_mil <- c(0.9896093164254245,
                  0.8835595872246769,
                  0.7451828617991133,
                  0.6168515063634167,
                  0.5009910992158976)

cdf_covg_mil <- c(0.9881801819795218,
                  0.9122502688589639,
                  0.8001428205234656,
                  0.653598382784491,
                  0.5419679935543181)

cdf_covg_sev <- c(1.0, 0.9996875, 0.963674281459094, 0.8363940495264706, 0.7210298540230009)
lac_covg_sev <- c(0.9702288842502945,
                  0.8821758195540386,
                  0.7922965180912996,
                  0.7178543052089326,
                  0.6598434498792998)
aps_covg_mod <- c(0.9901317784979403,
                  0.884418944754871,
                  0.7516118474214574,
                  0.6433459244789762,
                  0.5694407470014649)
lac_covg_mod <- c(0.9913785453491072,
                  0.8863679612883485,
                  0.7567346537628614,
                  0.6445976428957012,
                  0.5400122045523134)
cdf_covg_mod <- c(0.9886451954625407,
                  0.9042383371909816,
                  0.7958309313975924,
                  0.6757661325011368,
                  0.5771447776103454)
plt <- plot_ly(x=~alphas)
plt <- plt %>% add_lines(y=~aps_covg_neg,name="<b>Ordinal APS")
plt <- plt %>% add_lines(y=~lac_covg_neg,name="<b>Naive LAC")
plt <- plt %>% add_lines(y=~cdf_covg_neg,name="<b>Ordinal CDF")
plt <- plt %>% add_lines(y=~cond_coverage_opi[1,],name="<b>Ordinal Prediction Interval")
plt <- plt %>% add_lines(y=~cond_coverage_ops[1,],name="<b>Ordinal Prediction Set")
plt <- plt %>% layout(xaxis = list(title = "miscoverage"),
                      yaxis = list(title = "conditional coverage",range = c(0,1)),
                      title = "no stenosis")

plt <- plot_ly(x=~alphas)
plt <- plt %>% add_lines(y=~aps_covg_mod,name="<b>Ordinal APS")
plt <- plt %>% add_lines(y=~lac_covg_mod,name="<b>Naive LAC")
plt <- plt %>% add_lines(y=~cdf_covg_mod,name="<b>Ordinal CDF")
plt <- plt %>% add_lines(y=~cond_coverage_opi[3,],name="<b>Ordinal Prediction Interval")
plt <- plt %>% add_lines(y=~cond_coverage_ops[3,],name="<b>Ordinal Prediction Set")
plt <- plt %>% layout(xaxis = list(title = "miscoverage"),
                      yaxis = list(title = "conditional coverage",range = c(0,1)),
                      title = "moderate stenosis")

plt <- plot_ly(x=~alphas)
plt <- plt %>% add_lines(y=~aps_covg_sev,name="<b>Ordinal APS")
plt <- plt %>% add_lines(y=~lac_covg_sev,name="<b>Naive LAC")
plt <- plt %>% add_lines(y=~cdf_covg_sev,name="<b>Ordinal CDF")
plt <- plt %>% add_lines(y=~cond_coverage_opi[4,],name="<b>Ordinal Prediction Interval")
plt <- plt %>% add_lines(y=~cond_coverage_ops[4,],name="<b>Ordinal Prediction Set")
plt <- plt %>% layout(xaxis = list(title = "miscoverage"),
                      yaxis = list(title = "conditional coverage",range = c(0,1)),
                      title = "severe stenosis")


#==============================================================================#
#------------------------------------------------------------------------------#
#==============================================================================#

cdf_conditional_coverage <- NULL
aps_conditional_coverage <- NULL
lac_conditional_coverage <- NULL
opi_conditional_coverage <- NULL
ops_conditional_coverage <- NULL

for(alpha in 1:length(alphas))
{
  cdf_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_cdf[,alpha])
  aps_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_aps[,alpha])
  lac_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_lac[,alpha])
  opi_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_opi[,alpha])
  ops_conditional_coverage[alpha] <- max((1-alphas[alpha]) - cond_coverage_ops[,alpha])
}

conditional_coverage_plot <- plot_ly(x=~alphas)
conditional_coverage_plot <- conditional_coverage_plot %>% add_lines(y=~cdf_conditional_coverage,name="<b>Ordinal CDF")
conditional_coverage_plot <- conditional_coverage_plot %>% add_lines(y=~aps_conditional_coverage,name="<b>Ordinal APS")
conditional_coverage_plot <- conditional_coverage_plot %>% add_lines(y=~lac_conditional_coverage,name="<b>Naive LAC")
conditional_coverage_plot <- conditional_coverage_plot %>% add_lines(y=~opi_conditional_coverage,name="<b>Ordinal Prediction Interval")
conditional_coverage_plot <- conditional_coverage_plot %>% add_lines(y=~ops_conditional_coverage,name="<b>Ordinal Prediction Set")
conditional_coverage_plot <- conditional_coverage_plot %>% layout(xaxis = list(title = "miscoverage"),
                                                                yaxis = list(title = "conditional coverage",range = c(0,1)),
                                                                title = "Novel measure for conditional coverage")




