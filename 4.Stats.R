##### Purified enzyme stats

setwd("./PurifiedEnzymeData")

library(quantreg)
library(stringi)
library(multcomp)
library(lmtest)


### A) Quantile regression
## Effects of pH and temperature on the specific activity of APase
data <- read.table("Output/apase_kinetics.csv", sep = ",", row.names = "", header = T, check.names = F)
head(data)
data$pH <- as.factor(data$pH)
data$temp <- as.factor(data$temp)
str(data)

## Full model
r1 <- rq(SpecificActivity ~ pH*temp, data, tau = 0.5)
summary.rq(r1, iid = T, se = "iid")

## Model with only pH
r2 <- rq(SpecificActivity ~ pH, data, tau = 0.5)
## Model with only temperature
r3 <- rq(SpecificActivity ~ temp, data, tau = 0.5)
## Model with pH and temperature but without interaction
r4 <- rq(SpecificActivity ~ pH + temp, data, tau = 0.5)

## Wald test - significance of pH, temperature, and their interaction on the specific activity of APase
anova.rq(r1, r2, r3, r4, test = "Wald", iid = TRUE)
rm(list = ls())


### B) Log-likelihood ratios
## Effects of pH on the temperature sensitivity of APase
data <- read.table("Output/apase_kinetics.csv", sep = ",", row.names = "", header = T, check.names = F)
head(data)
data$pH <- as.factor(data$pH)
data$temp <- as.factor(data$temp)
data$revK <- 1/(as.numeric(as.character(data$temp))+273.15)
str(data)

## 1) Creating X matrix
# Creating an empty matrix X
full_Xmat <- matrix(ncol = 10, nrow = 120, data = 0)
## nrow is equal the number of observations
nrow(full_Xmat)
## ncol is equal the number of factors x2 (1 intercept and 1 slope each); the first 5 are the intercepts, and the other 5 are the slopes
ncol(full_Xmat)

phs <- levels(data$pH)
for (k in 1:120){
  if(data$pH[k] == "3.5"){                                                                
    i <- 1                                                                                
  } else if (data$pH[k] == "4.5") {
    i <- 2
  } else if (data$pH[k] == "5.5"){
    i <- 3
  } else if (data$pH[k] == "6.5"){
    i <- 4
  } else {
    i <- 5
  }
  ## the i value corresponds to the intercept column of the pH of the row
  full_Xmat[k, i] <- 1
  ## the [5+i] value corresponds to the slope column of the pH of the row
  full_Xmat[k, 5 + i] <- data$revK[k]
}
## intercepts have either 1 or 0 to indicate the right factor; the slopes have the temperatures (as 1/T)
head(full_Xmat)

## 2) Separating the intercept and the slope matrices
# a matrix only for the intercepts
intercept <- full_Xmat[,1:5]

for(i in 6:10){
  v <- c(3, 4, 5, 6, 7)
  nm <- paste0("pH", v[i-5], "5")
  assign(nm, as.matrix(full_Xmat[,i]))
}
ls()[grep("pH", ls())]

## 3) Creating Y matrix
# the Y matrix is the response variable
Y_mat <- as.matrix(data$lnactivity*(-8.314))
head(Y_mat)

rm(i, k, nm, v)

## 4) Creating comparison matrices
head(intercept)
phs <- ls()[grep("pH", ls())]
checklist <- data.frame()

i <- 1
j <- 0
for(p in phs){
  j <- 1 + j
  for(h in phs[-1:-j]){
    nm <- paste0("Xmat", i)
    assign(nm, cbind(intercept, (get(p) + get(h)), 
                     get(phs[phs != p & phs!= h][1]),
                     get(phs[phs != p & phs!= h][2]),
                     get(phs[phs != p & phs!= h][3])))
    checklist[i, "Matrix"] <- nm
    checklist[i, "Pair"] <- paste(p, "-", h)
    i <- i + 1
  }
}
ls()[grep("Xmat", ls())]
rm(i, j, nm, p, h)

## 5) Fitting the linear model to the full dataset
lm_fullXmat <- lm.fit(full_Xmat, Y_mat)
residuals(lm_fullXmat)

## 6) Fitting the linear model to the comparison matrices
xmats <- stri_sort(ls()[grep("Xmat\\d{1,2}", ls())], numeric = TRUE)
xmats

for(x in xmats){
  nm <- paste0('lm_', x)
  assign(nm, lm.fit(get(x), Y_mat))
}
stri_sort(ls()[grep("lm_Xmat", ls())], numeric = TRUE)
rm(x, xmats)

## 7) Finding the log-likelihood
# creating a vector for the bias correction for the standard deviation of the residuals
n_bias <- ((120 - 1)/120)

## 7.1) Finding the LL of the full model
# calculating the log-likelihood for the full model
LL_fullXmat <- sum(dnorm(residuals(lm_fullXmat), mean = 0, 
                         sd = sqrt(n_bias*(var(lm_fullXmat$residuals))), log = TRUE))
LL_fullXmat

## 7.2) Finding the LL of the comparison matrices
# creating a vector with the names of the linear models of the comparison matrices in order
lmxmats <- stri_sort(ls()[grep("lm_Xmat\\d{1,2}", ls())], numeric = TRUE)
lmxmats

for(l in lmxmats){
  ## using a loop to find the LL of each comparison matrix
  nm <- paste0("LL_", gsub("lm_(Xmat{1,2})", "\\1", l))
  assign(nm, sum(dnorm(residuals(get(l)), mean = 0, 
                       sd = sqrt(n_bias*(var(residuals(get(l))))), log = TRUE)))
}
stri_sort(ls()[grep("LL_Xmat", ls())], numeric = TRUE)
rm(lmxmats, l, nm)

## 8) Finding the log-likelihood ratio and p-values
checklist
LLres <- stri_sort(ls()[grep("LL_Xmat", ls())], numeric = TRUE)
LLres

for(LL in LLres){
  ## Comparing the LL of each comparison matrix (a simplified model) to the full model
  checklist[checklist$Matrix == gsub("LL_(Xmat{1,2})", "\\1", LL), 
            "LR"] <- -2*(get(LL) - LL_fullXmat)
  ## Using the ratio of the LL of these matrices and using the Chi-square distribution to find the p-value
  checklist[checklist$Matrix == gsub("LL_(Xmat{1,2})", "\\1", LL), 
            "pval"] <- (1 - pchisq(checklist[checklist$Matrix == gsub("LL_(Xmat{1,2})", "\\1", LL), 
                                             "LR"], df = 1))
}
checklist$sig <- checklist$pval <= 0.05
checklist[checklist$sig == TRUE, ]
rm(list = ls())