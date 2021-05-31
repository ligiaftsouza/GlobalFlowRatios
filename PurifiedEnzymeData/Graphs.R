##### Purified enzymes plots

setwd("./PurifiedEnzymeData")


### Fig. 1
## Calculating the average per temperature and pH
data <- read.csv("Output/apase_kinetics.csv", check.names = F, row.names = 1)
head(data)
phs <- unique(data$pH)
temps <- unique(data$temp)
avg_apase <- data.frame()
k <- 1
for(p in phs){
  for(t in temps){
    avg_apase[k, "pH"] <- p
    avg_apase[k, "temp"] <- t
    avg_apase[k, "avg_SpAct"] <- mean(data$SpecificActivity[data$pH == p & data$temp == t], na.rm = T)
    avg_apase[k, "sd_SpAct"] <- sd(data$SpecificActivity[data$pH == p & data$temp == t], na.rm = T)
    k <- k + 1
  }
}
avg_apase

png("Plots/Fig1.png", width = 1600, height = 1300, pointsize = 60)
par(oma = c(0, 0, 0, 0), mar = c(3, 4, 0.5, 1))
plot(x = avg_apase$pH, y = avg_apase$avg_SpAct, type = "n", xlab = "", 
     ylab = "", ylim = c(0, 36), xaxt = "n", cex = 0.9, 
     cex.axis = 0.9)
axis(1, at = seq(3.5, 7.5, 1), cex = 0.9, cex.axis = 0.9)
mtext(side = 1, "pH", line = 2, cex = 0.9)
mtext(side = 2, (expression(paste("APase specific activity, ", mu, "mol h"^-1*"mg"["enzyme"]^-1*""))), 
      line = 1.8, cex = 0.9)
pchtemps <- c(21, 22, 24, 25)
names(pchtemps) <- c(5, 15, 25, 35)
ltytemps <- c(1, 2, 3, 5)
names(ltytemps) <- c(5, 15, 25, 35)
bgtemps <- c("white", "black", "white", "black")
names(bgtemps) <- c(5, 15, 25, 35)
for(t in temps){
    arrows(x0 = avg_apase$pH[avg_apase$temp == t], 
           x1 = avg_apase$pH[avg_apase$temp == t], 
           y0 = avg_apase$avg_SpAct[avg_apase$temp == t] - avg_apase$sd_SpAct[avg_apase$temp == t], 
           y1 = avg_apase$avg_SpAct[avg_apase$temp == t] + avg_apase$sd_SpAct[avg_apase$temp == t], 
           lty = 1, angle = 90, code = 3, length = 0.05, lwd = 2, 
           col = "black")
    lines(x = avg_apase$pH[avg_apase$temp == t], 
          y = avg_apase$avg_SpAct[avg_apase$temp == t], col = "black", 
          lwd = 2, lty = ltytemps[names(ltytemps) == t])
    points(x = avg_apase$pH[avg_apase$temp == t],
           y = avg_apase$avg_SpAct[avg_apase$temp == t],
           bg = bgtemps[names(bgtemps) == t], pch = pchtemps[names(pchtemps) == t], 
           cex = 1.6)
}
legend("bottomleft", col = "black", pch = rev(pchtemps), pt.bg = rev(bgtemps),
       c(expression(35~degree~C), expression(25~degree~C), expression(15~degree~C), expression("  5"~degree~"C")),
       lty = rev(ltytemps), lwd =2, cex = 0.6, pt.cex = 1.2, 
       ncol = 2, bty = "n", seg.len = 3)
text(x = 7.0, y = 33, c("pH: p < 0.001\ntemperature: p < 0.001\npH*temperature: p < 0.001"), cex = 0.4)
dev.off()
rm(list = ls())


### Fig. 2
data <- read.csv("Output/apase_kinetics.csv", check.names = F, row.names = 1)
head(data)
phs <- unique(data$pH)
temps <- unique(data$temp)
avg_apase <- data.frame()
k <- 1
for(p in phs){
  for(t in temps){
    avg_apase[k, "pH"] <- p
    avg_apase[k, "temp"] <- t
    avg_apase[k, "revK"] <- 1/(t + 273.15)
    avg_apase[k, "avg_lnSpAct"] <- mean(log(data$SpecificActivity[data$pH == p & data$temp == t]), na.rm = T)
    avg_apase[k, "sd_lnSpAct"] <- sd(log(data$SpecificActivity[data$pH == p & data$temp == t]), na.rm = T)
    k <- k + 1
  }
}

png("Plots/Fig2.png", width = 1600, height = 1300, pointsize = 60)
par(oma = c(0, 0, 0, 0), mar = c(3, 4, 0.5, 1))
plot(x = avg_apase$revK, y = avg_apase$avg_lnSpAct, 
     ylab = "", xlab = "", pch = c(1, 19, 2, 17, 15)[as.factor(avg_apase$ph)],
     ylim = c(0,4), xaxt = "n", yaxt = "n", cex = 1.2, cex.axis = 0.8,
     col = 'black', bg = "black", type = "n")
axis(1, at = round(unique(avg_apase$revK), digits = 5), cex.axis = 0.8, mgp = c(0, 0.7, 0), cex = 0.8)
axis(2, at = seq(0, 4, 1), mgp = c(0, 0.7, 0), cex.axis = 0.8, cex = 0.8)
mtext(side = 1, "1/K", line = 1.8, cex = 0.85)
mtext(side = 2, (expression(paste("ln(APase specific activity), ", mu,
                                  "mol h"^-1*"mg"["enzyme"]^-1*" "))), line = 2, cex = 0.85)
ltys <- c(6, 1, 3, 4, 5)
names(ltys) <- rev(phs)
pts <- c(15, 21, 19, 24, 17)
names(pts) <- rev(phs)
bgs <-  c("black", "white", "black", "white", "black")
names(bgs) <- rev(phs)
for(p in phs){                                                       ## using a loop to plot the error bars (using SEM but it can be changed)
  arrows(x0 = avg_apase$revK[avg_apase$pH == p], x1 = avg_apase$revK[avg_apase$pH == p], 
         y0 = avg_apase$avg_lnSpAct[avg_apase$pH == p] - avg_apase$sd_lnSpAct[avg_apase$pH == p], 
         y1 = avg_apase$avg_lnSpAct[avg_apase$pH == p] + avg_apase$sd_lnSpAct[avg_apase$pH == p], 
         lty = 1, angle = 90, code = 3, length = 0.05, col = "black", 
         lwd = 0.5)
  clip(min(avg_apase$revK), max(avg_apase$revK), 0, 5)                                  ## adding a clip to plot lines within the avg_apase points
  abline(lm(avg_apase$avg_lnSpAct[avg_apase$pH == p] ~avg_apase$revK[avg_apase$pH == p]), 
         col = "black", 
         lwd = 2, lty = ltys[names(ltys) == p])  
  do.call("clip", as.list(par("usr")))
  points(x = avg_apase$revK[avg_apase$pH == p], y = avg_apase$avg_lnSpAct[avg_apase$pH == p], 
         pch = pts[names(pts) == p], cex = 1.2,
         col = "black", 
         bg = bgs[names(bgs) == p])
}
legend("bottomleft", col = "black", 
       pch = pts, pt.bg = bgs, lty = ltys, lwd = 2,
       c("pH 3.5", "pH 4.5", "pH 5.5", "pH 6.5", "pH 7.5"), 
       cex = 0.75, pt.cex = 0.85, ncol = 2,  bty = "n")                                                             ## adding a legend to the plot
dev.off()     
rm(list = ls())


## Fig. 3
flowratios <- read.csv("Output/flowratios.csv", row.names = 1)
head(flowratios)
pts <- c(21, 22, 23)
names(pts) <- c("5", "15", "25")
ltys <- c(1, 2, 3)
names(ltys) <- c("5", "15", "25")
bgs <- c("black", "white", "black")
names(bgs) <- c("5", "15", "25")
temps <- c(5, 15, 25)
enzpairs <- c("bgase:apase", "nagase:apase", "bgase:nagase")
names(enzpairs) <- c("C:P", "N:P", "C:N")

png(paste0("Plots/Fig3.png"), width = 1000, height = 2000, pointsize = 60)
par(mfrow = c(1, 3))
layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1), heights = c(0.8, 0.8, 1))
layout.show(3)
for (e in enzpairs){
  ep <- names(enzpairs[enzpairs == e])
  nm <- gsub(":", "", ep)
  if(e == "bgase:apase"){
    par(mar = c(0, 4, 0.5, 0.5))
    ast <- c(4.5, 5.5, 6.5, 7.5)
    lb <- "(a)"
  } else if(e == "nagase:apase"){
    par(mar = c(0, 4, 0, 0.5))
    ast <- c(4.5, 5.5, 6.5, 7.5)
    lb <- "(b)"
  } else{
    par(mar = c(3.5, 4, 0, 0.5))
    ast <- c(4.5, 7.5)
    lb <- "(c)"
  }
  plot(x = flowratios$ph[flowratios$Pair == ep], 
       y = flowratios$ratio[flowratios$Pair == ep], 
       type = "n", col = "black", lwd = 3, 
       ylim = c(-1, ceiling(round(max(flowratios$ratio[flowratios$Pair == ep]), 1))*1.2),
       ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  if(e == "bgase:nagase"){
    axis(1, at = seq(4.5, 7.5, 1), labels = T, mgp = c(0, 0.8, 0), cex.axis = 0.8)
    mtext("pH", side = 1, cex = 0.7, line = 2.3)
  }
  axis(2, labels = T, mgp = c(0, 0.8, 0), cex.axis = 0.8)
  mtext(paste(ep, "\nestimated flow ratio"), side = 2, 
        mgp = c(0, 0.8, 0), cex = 0.7, line = 1.8)
  text(x = 4.5, y = max(flowratios$ratio[flowratios$Pair == ep])*1.2, lb, xpd = NA, font = 2)
  for (t in temps){
    lines(x = flowratios$ph[flowratios$Pair == ep & flowratios$tempC == t], 
          y = flowratios$ratio[flowratios$Pair == ep & flowratios$tempC == t], 
          col = "black", lwd = 2, lty = ltys[names(ltys) == t])
    points(x = flowratios$ph[flowratios$Pair == ep & flowratios$tempC == t], 
           y = flowratios$ratio[flowratios$Pair == ep & flowratios$tempC == t], 
           bg = bgs[names(bgs) == t], col = "black",
           pch = pts[names(pts) == t], cex = 1.6, lwd = 1.2)
  }
  for(a in ast){
    if(max(flowratios$ratio[flowratios$Pair == ep & flowratios$ph == a]) < ceiling(round(max(flowratios$ratio[flowratios$Pair == ep])))/2){
      ymax <- ceiling(round(max(flowratios$ratio[flowratios$Pair == ep]), 1))*0.4
    } else{
      ymax <- max(flowratios$ratio[flowratios$Pair == ep & flowratios$ph == a])*1.15
    }
    points(x = a, 
         y = ymax, 
         bg = bgs[names(bgs) == t], col = "black",
         pch = 8, cex = 0.7, lwd = 1.5)
    }
  legend("bottomright", legend =c(expression(paste(" 5", degree, "C")), 
                              expression(paste("15", degree, "C")),
                              expression(paste("25", degree, "C"))), 
         col = 'black', pch = pts, lty = ltys, lwd = , seg.len = 2, bty = "n", cex = 0.8,
         pt.bg = bgs, pt.cex = 1.4)
}
dev.off()
rm(list = ls())


### Fig. 4
data <- read.csv("Output/arrdata.csv", row.names = "", header = T, check.names = F) 
data <- data[data$ph != 3.5 & data$ph != 8.5, ]
enz <- c("APase", "BGase", "NAGase")
ltyenz <- c(3, 1, 2)
names(ltyenz) <- enz

png("Plots/Fig4.png", width = 1600, height = 1300, 
    pointsize = 60)
par(oma = c(0, 0, 0, 0), mar = c(3, 3, 0.5, 1))
x <- seq(3.5, 7.5, 0.1)
slopeSpline <- splinefun(x = data$ph, y = data$Ea)
plot(x, slopeSpline(x), type = "n", xlab = "", ylab = "", 
     ylim = c(0, 75), xlim = c(4.5, 7.5), xaxt = "n", yaxt = "n", cex = 0.9, xaxs = "i")
axis(1, at = seq(3.5, 7.5, 1), cex = 0.8, cex.axis = 0.8, mgp = c(0, 0.7, 0)) 
axis(2, at = seq(0, 75, 10), cex.axis = 0.8, mgp = c(0, 0.7, 0), cex = 0.8)
axis(3, labels = FALSE, tick = FALSE)
mtext(side = 1, "pH", line = 1.7, cex = 0.85)
mtext(expression(paste("E"["a"]*", kJ mol"^-1*"")), side = 2, line = 1.7, cex = 0.85)
# text(x = 3.5, y = 25, "(b)")
for(e in enz){
  slopeSpline <- splinefun(x = data$ph[data$enzyme == tolower(e)], 
                           y = data$Ea[data$enzyme == tolower(e)])
  slopeSplineSD <- splinefun(x = data$ph[data$enzyme == tolower(e)], 
                             y = (-data$sd[data$enzyme == tolower(e)]))
  x <- seq(3.5, 7.5, 0.1)
  x1 <- seq(3.5, 7.49, 0.01)
  y1 <- slopeSpline(x) + slopeSplineSD(x)
  y2 <- slopeSpline(x) - slopeSplineSD(x)
  polygon(c(x,rev(x)), c(y2,rev(y1)), col = adjustcolor("lightgray", alpha.f = 0.5), border = NA)
  lines(x = x1, y = slopeSpline(x1), lwd = 3, lty = ltyenz[names(ltyenz) == e])
}
box(which = "plot", lty = "solid")
legend("topright", legend = enz,
       lty = ltyenz, lwd = 3, bty = "n", cex = 0.7, seg.len = 1.6)
dev.off()
rm(list = ls())
