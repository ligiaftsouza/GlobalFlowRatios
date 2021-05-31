##### APase kinetics calculations

setwd("./PurifiedEnzymeData")
library(inflection)

phs <- c("ph75", "ph65", "ph55", "ph45", "ph35")
index <- data.frame()
apase_well <- 0.0256/3.1
slopes <- data.frame()
i <- 1

for(p in phs){
  data <- read.csv(paste0("Data/", p,".csv"), header = T, check.names = F)
  for(c in colnames(data)[colnames(data) != "Time"]){
    smoothedcurve <- loess.smooth(x = data[, "Time"], y = data[,c], span = 0.25, evaluation = 31)
    smoothedcurve$x <- data[, "Time"]
    infpt <- ede(x = smoothedcurve$x, y = smoothedcurve$y, index = 0)[3]
    index[c, p] <- 0
    if(is.na(infpt) == TRUE){
      infpt <- ede(x = smoothedcurve$x, y = smoothedcurve$y, index = 1)[3]
      index[c, p] <- 1
    }
    write.csv(assign(paste0(p, "_", c), data[data$Time <= infpt, c("Time", c)]),
              file = paste0("Output/LinearAccumulation/loess", 
                            paste0(p, "_", c), ".csv"), quote = F, row.names = F)
    tmp <- get(paste0(p, "_", c))
    names(tmp)[2] <- "Curve"
    lm_curve <- lm(Curve ~ Time, data = tmp)
    slopes[i, "pH"] <- as.numeric(gsub(pattern = "^ph(\\d)(\\d)", replacement = "\\1.\\2", p))
    slopes[i, "temp"] <- as.numeric(gsub(pattern = "(\\d{1,2})_\\d{1}", replacement = "\\1", c))
    slopes[i, "Intercept"] <- lm_curve$coefficients[1]
    slopes[i, "Slope"] <- lm_curve$coefficients[2]
    ## finding the specific enzyme activity in umol h-1 mg-1 of enzyme in the well; dividing by 1000 to transform from nmol to umol
    slopes[i, "SpecificActivity"] <- slopes[i, "Slope"]/1000/apase_well
    slopes[i, "lnactivity"] <- log(slopes[i,"SpecificActivity"])
    rownames(slopes)[i] <- paste0(p, "_", c)
    i <- i + 1
  }
  
}
slopes
write.csv(slopes, "Output/apase_kinetics.csv", quote = F, row.names = T)
rm(list = ls())
