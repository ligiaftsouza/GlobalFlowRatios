##### APase, BGase, and NAGase dataset

setwd("./PurifiedEnzymeData")


### Creating a combined dataset for all enzymes
## BGase and NAGase data - Min et al. (2014)
data <- read.csv("Data/pH_enzyme.csv", sep = ",", check.names = F, header = T)
data$revK <- 1/(data$temp + 273.15)
str(data)
## APase data
apase <- read.csv("Output/apase_kinetics.csv", sep = ",", header = T, row.names = 1, check.names = F)
apase$revK <- 1/(apase$temp + 273.15)
apase$enzyme <- "APase"
apase <- apase[, c(8, 2, 1, 5, 6, 7)]
colnames(apase) <- colnames(data)
head(apase)
## Combined datasets
allenz <- rbind(data, apase)
write.csv(allenz, "Output/allenz.csv")
rm(list = ls())


### Calculating the activation energy
allenz <- read.csv("Output/allenz.csv", row.names = 1)
enzs <- unique(allenz$enzyme)
arrdata <- data.frame()
k <- 1   
for(e in enzs){
  tmp <- allenz[allenz$enzyme == e, ]
  if(e == "APase"){
    phs <- c(3.5, 4.5, 5.5, 6.5, 7.5)
  } else{
    phs <- c(4.5, 5.5, 6.5, 7.5, 8.5)
  }
  for (p in phs){
    lmph <- lm(lnactivity ~ revK, data = tmp[tmp$pH == p, ])
    out <- summary(lmph)
    arrdata[k, "enzyme"] <- tolower(e)
    arrdata[k, "ph"] <- p
    arrdata[k, "intercept"] <- lmph$coefficients[1]
    arrdata[k, "slope"] <- lmph$coefficients[2]
    ## Activation energy by multiplying the slope by R (the gas constant, which is 8.314 J-1 mol-1); the division by 1000 is to transform J mol-1 to kJ mol-1
    arrdata[k, "Ea"] <- (lmph$coefficients[2]*(-8.314))/1000
    ## getting the standard deviation in the same unit as the Ea
    arrdata[k, "sd"] <- abs(out$coefficients[2,2]*(-8.314))/1000
    k <- k + 1
  }
}
write.csv(arrdata, "Output/arrdata.csv")
rm(list = ls())

