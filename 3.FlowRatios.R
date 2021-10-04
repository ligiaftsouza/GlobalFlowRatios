##### Calculation of flow ratios

setwd("./PurifiedEnzymeData")


### Calculating the flow ratios
data <- read.csv("Output/arrdata.csv", check.names = F, row.names = "")
data$A <- exp(data$intercept)
data
## Calculating Vmax
kvmax <- data.frame()
enzymes <- unique(data$enzyme)
temps <- c(5, 15, 25)
k <- 1
for (e in enzymes){
  phs <- unique(data$ph[data$enzyme == e])
  for (p in phs){
    for (t in temps){
      ## exponential modification of the Arrhenius equation
      ex <- exp(data$Ea[data$enzyme == e & data$ph == p]/-0.008314/(273.15 + t))
      kvmax[k, "enzyme"] <- e
      kvmax[k, "pH"] <- p
      kvmax[k, "temp"] <- t
      kvmax[k, "kvmax"] <- (data$A[data$enzyme == e & data$ph == p])*ex
      k <- k + 1
    }
  }
}
flowratios <- data.frame()
enzpairs <- c("bgase:nagase", "bgase:apase", "nagase:apase")
names(enzpairs) <- c("C:N", "C:P", "N:P")
phs <- c(4.5, 5.5, 6.5, 7.5)
k <- 1
for (e in enzpairs){
  a <- gsub(pattern = "(\\w{5,6})\\W\\w{5,6}", replacement = "\\1", x = e)
  b <- gsub(pattern = "\\w{5,6}\\W(\\w{5,6})", replacement = "\\1", x = e)
  for (p in phs){
    for (t in temps){
      ratio <- (kvmax$kvmax[kvmax$enzyme == a & kvmax$pH == p & kvmax$temp == t])/
        (kvmax$kvmax[kvmax$enzyme == b & kvmax$pH == p & kvmax$temp == t])
      flowratios[k, "Pair"] <- names(enzpairs[enzpairs == e])
      flowratios[k, "ph"] <- p
      flowratios[k, 'tempC'] <- t
      if (e == "bgase:nagase"){
        flow <- ratio*6 + 8
      } else if (e == "bgase:apase"){
        flow <- ratio*6
      } else {
        flow <- ratio*1
      }
      flowratios[k, "ratio"] <- flow
      flowratios[k, "logratio"] <- log(flow)
      k <- k + 1
    }
  }
}
write.csv(flowratios, "Output/flowratios.csv")
rm(list = ls())

##### Calculating the CV of the flow ratios per enzyme pair and pH
flowratios <- read.csv("Output/flowratios.csv", row.names = 1)
fr_cv <- data.frame()
enzpairs <- c("bgase:nagase", "bgase:apase", "nagase:apase")
names(enzpairs) <- c("C:N", "C:P", "N:P")
phs <- unique(flowratios$ph)
k <- 1
for(e in enzpairs){
  for(p in phs){
    tmp <- flowratios[flowratios$Pair == names(enzpairs[enzpairs == e]) & 
                        flowratios$ph == p, ]
    fr_cv[k, "pair"] <- names(enzpairs[enzpairs == e])
    fr_cv[k, "ph"] <- p
    fr_cv[k, "avg_fr"] <- mean(tmp$ratio)
    fr_cv[k, "sd_fr"] <- sd(tmp$ratio)
    fr_cv[k, "cv_fr"] <- sd(tmp$ratio)/mean(tmp$ratio)
    k <- k + 1
  }
}
fr_cv$cvperc <- fr_cv$cv_fr*100
fr_cv
write.csv(fr_cv, "Output/flowratios_cv.csv")
rm(list = ls())
