# 1. generate p-features from standard normal distribution
# 2. they are correlated with each other
# 3. based on the features simulate y

# beta_0 = -7, beta_1 = beta_2 = beta_3 = 4, beta_4 = -6 rad(2), beta_5 = 4/3
# beta_p = 0, p > 6

# Simulating X

#install.packages("ggplot2")


library(ggplot2)
library(MASS)

d <- 5 # dimensions
N <- 1000 # number of generated samples

#-----------------------Construction of covariance matrix----------------------#

covariance.matrix <- matrix(1/2, nrow = d, ncol = d)

index <- c(1:d)

covariance.matrix[4,index[!index %in% 4]] = 1/sqrt(2)
covariance.matrix[index[!index %in% 4],4] = 1/sqrt(2)

covariance.matrix[5,index[!index %in% 5]] = 0
covariance.matrix[index[!index %in% 5],5] = 0

for(i in 1:d)
{
  covariance.matrix[i,i] <- 1
}

# View(covariance.matrix)

X <- mvrnorm(N,rep(0,d),covariance.matrix)
X <- cbind(rep(1,nrow(X)),X)
beta_hat <- c(0,rep(4,3),-6*sqrt(2),4/3,rep(0,d-5))

# Generating the response

product.function <- X %*% beta_hat
df <- data.frame(cbind(X[,-1],product.function))

df$odds <- exp(product.function)/(1 + exp(product.function))
df$class <- ifelse(df$odds > 0.5,1,0)

# View(df)

names(df)[d+1] <- "logits"
names(df)[d+2] <- "pi_hat"
names(df)[d+3] <- "label"

df_1 <- df[,-c(d+1,d+2)]
write.csv(df_1,"/Users/subhrasishchakraborty/Desktop/df_1.csv")

getwd()
setwd("/Users/subhrasishchakraborty/Desktop")
results <- read.csv("rs.csv",header = TRUE)

View(results)

plt <- ggplot(results,aes(x = d)) + geom_line(aes(y = auc, color = beta_0))
plt <- plt + labs(x = "dimension", y = "AUC", color = "Legend")
