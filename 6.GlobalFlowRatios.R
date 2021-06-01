#####

setwd("./FlowRatiosProjections/")

library(raster)
library(rgdal)
library(leaflet)
library(viridisLite)


### Calculating global temperature sensitivity
data <- read.csv("../PurifiedEnzymeData/Output/arrdata.csv", row.names = "", header = T, check.names = F)           ## importing the Arrhenius calculations data with the slope (Ea) and intercept to
data <- data[data$ph != 3.5 & data$ph != 8.5, ]

ph_mask <- function(x){ 
  x[x > 7.5] <- NA
  x[x < 4.5] <- NA
  return(x)
}

phs <- dir("Data/soilph/")[grep(".tif$", dir("Data/soilph/"))] 
phs
enz <- c("apase", "nagase", "bgase")
for(e in enz){
  s_enz <- data[data$enzyme == e, ]
  ## Creating a slope spline function
  slopeSpline <- splinefun(x = s_enz$ph, y = s_enz$Ea)
  ## Creating an intercept spline function
  intSpline <- splinefun(x = s_enz$ph, y = s_enz$intercept)
  for(p in phs){
    writeRaster(ph_mask(calc(raster(paste0("Data/soilph/", p)), function(x) x/10)), 
                paste0("Output/ph/m", p), overwrite = TRUE)
    nm <- paste0(gsub("(\\w{2})\\w{3,4}", "\\1", e), gsub("ph(\\d{1,3})cm.tif", "\\1", p))
    ## Calculating the Ea globally
    writeRaster(calc(raster(paste0("Output/ph/m", p)), fun = function(x) slopeSpline(x)), 
                paste0("Output/spatialEa/", nm, "ea.tif"), overwrite = TRUE)
    ## Calculating ln(A) globally
    writeRaster(calc(raster(paste0("Output/ph/m", p)), fun = function(x) intSpline(x)), 
                paste0("Output/spatiallnA/", nm, "lnA.tif"), overwrite = TRUE)
  }
}
rm(list = ls())


### Calculating enzyme activity
depths <- c('000', "030", "100")
models <- dir("Data/projtemp/")
enz <- c("apase", "nagase", "bgase")
for(m in models){
  prjs <- dir(paste0("Data/projtemp/", m))[grep("4k.tif$", dir(paste0("Data/projtemp/", m)))]
  for(p in prjs){
    mpr <- raster(paste0("Data/projtemp/", m, "/", p))
    ssp <- gsub("\\w{1,3}\\d{0,1}(\\d{3})_\\d{1}k", "\\1", names(mpr))
    htK <- raster("Data/hist_temp/hmatk.tif")
    ptK <- mpr
    R <- 8.314/1000
    for(e in enz){
      for(d in depths){
        slp <- raster(paste0("Output/spatialEa/", gsub("(\\w{2})\\w{3,4}", "\\1", e), d, "ea.tif"))
        intercept <- raster(paste0("Output/spatiallnA/", gsub("(\\w{2})\\w{3,4}", "\\1", e), d, "lnA.tif"))
        ## calculating enzyme activity based on historical temperatures
        tmp <- overlay(slp, intercept, htK, fun = function(x, y, z){(-x/R)*(1/z) + y})
        writeRaster(tmp, filename = paste0("Output/activity/historical/", e, "/h", 
                                           gsub("(\\w{2})\\w{3,4}", "\\1", e), d, ".tif"),  
                    overwrite = TRUE)
        ## calculating enzyme activity based on temperature projections
        tmp <- overlay(slp, intercept, ptK, fun = function(x, y, z){(-x/R)*(1/z) + y})
        writeRaster(tmp, filename = paste0("Output/activity/projected/", e, "/", m, "/", 
                                           gsub("(\\w{2})\\w{3,4}", "\\1", e), d, tolower(m), ssp, ".tif"), 
                    overwrite = TRUE)
      }
    }
  }
}
rm(list = ls())


### Calculating average enzyme activity for all temperature projections
depths <- c("000", "030", "100")
ssp <- c(245, 585)
enz <- c("apase", "bgase", "nagase")
for(d in depths){
  for(e in enz){
    prjs <- list.files(path = paste0("Output/activity/projected/", e, "/"), pattern = d, 
                       include.dirs = TRUE, recursive = TRUE)
    for(s in ssp){
      models <- prjs[grep(s, prjs)]
        for(m in models){
          assign(paste0("m", which(models == m)), exp(raster(paste0("Output/activity/projected/", e, "/", m))))
        }
        files <- ls()[grep("m\\d{1,2}", ls())]
        avg <- calc(brick(mget(files)), fun = function(x) mean(x, na.rm = TRUE))
        writeRaster(avg, filename = paste0("Output/activity/projected/averages/avg_", 
                                           gsub("(\\w{2})\\w{3,4}", "\\1", e), d, s, ".tif"),
                    overwrite = TRUE)
        rm(list = ls()[grep("m\\d{1,2}", ls())])
    }
  }
}
rm(list = ls())


### Calculating the relative change in flow ratios
depths <- c("000", "030", "100")
es <- c("b", "n")
names(es) <- c("bgase", "nagase")
els <- c("c", "n", "p")
names(els) <- c("b", "n", "a")
ssp <- c(245, 585)
relchange <- data.frame()
k <- 1
for(d in depths){
  for(e in es){
    if(e == "b"){
      s_es <- c("n", "a")
      names(s_es) <- c("nagase", "apase")
    } else {
      s_es <- c("a")
      names(s_es) <- c("apase")
    }
    for(s in s_es){
      if(e == "b" & s == "n"){
        rmt <- function(x){x*6 + 8}
      } else if(e == "b" & s == "a"){
        rmt <- function(x){x*6}
      } else{
        rmt <- function(x){x*1}
      }
      ## calculating flow ratios beased on historical temperatures
      e1 <- exp(raster(paste0("Output/activity/historical/", names(es[es == e]), "/h", 
                          sub("(\\w{2})\\w{3,4}", "\\1", names(es[es == e])), d, ".tif")))
      e2 <- exp(raster(paste0("Output/activity/historical/", names(s_es[s_es == s]), "/h", 
                          sub("(\\w{2})\\w{3,4}", "\\1", names(s_es[s_es == s])), d, ".tif")))
      hratio <- rmt(e1/e2)
      writeRaster(hratio, filename = paste0("Output/ratios/historical/", 
                                            els[names(els) == e], els[names(els) == s], "/h",
                                           els[names(els) == e], els[names(els) == s], d, ".tif"),
                  overwrite = TRUE)
      for(sp in ssp){
        ## calculating flow ratios beased on projected temperatures
        e1 <- raster(paste0("Output/activity/projected/averages/avg_",sub("(\\w{2})\\w{3,4}", "\\1", names(es[es == e])),
                             d, sp, ".tif"))
        e2 <- raster(paste0("Output/activity/projected/averages/avg_", sub("(\\w{2})\\w{3,4}", "\\1", names(s_es[s_es == s])),
                       d, sp, ".tif"))
        pratio <- rmt(e1/e2) ## the rasters are not in log scale; they were transformed when getting the average
        writeRaster(pratio, filename = paste0("Output/ratios/projected/", 
                                              els[names(els) == e], els[names(els) == s], "/p",
                                              els[names(els) == e], els[names(els) == s], d, sp, ".tif"),
                    overwrite = TRUE)
        ## calculating the relative change in the flow ratios by the end of the century
        tmp <- (pratio - hratio)/hratio*100
        writeRaster(tmp, filename = paste0("Output/ratios/relchange/", 
                                              els[names(els) == e], els[names(els) == s], "/rc",
                                              els[names(els) == e], els[names(els) == s], d, sp, ".tif"),
                    overwrite = TRUE)
        tmp2 <- tmp[!is.na(tmp)]
        relchange[k, "ratio"] <- paste0(els[names(els) == e], els[names(els) == s])
        relchange[k, "ssp"] <- sp
        relchange[k, "depth"] <- d
        relchange[k, "min"] <- min(tmp2)
        relchange[k, "max"] <- max(tmp2)
        relchange[k, "mean"] <- mean(tmp2)
        relchange[k, "median"] <- median(tmp2)
        relchange[k, "positiveChangeArea"] <- length(tmp2[tmp2 > 0])/length(tmp2)*100
        relchange[k, "greaterThanMeanArea"]<- length(tmp2[tmp2 > mean(tmp2)])/length(tmp2)*100
        k <- k + 1
      }
    }
  }
}
write.csv(relchange, "Output/relativeflowratiochangestats.csv", row.names = F)
