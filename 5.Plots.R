##### Purified enzymes plots

setwd("./PurifiedEnzymeData")

library(tidyverse)
library(cowplot)

### Fig. 1a and 1b
## APase kinetics
data <- read_csv("Output/apase_kinetics.csv")
head(data)

avg_apase <- data %>% group_by(pH, temp) %>% 
  summarise(avg_SpAct = mean(SpecificActivity),
            sd_SpAct = sd(SpecificActivity),
            avg_lnSpAct = mean(log(SpecificActivity)),
            sd_lnSpAct = sd(log(SpecificActivity)))

p1 <- ggplot(avg_apase, aes(x = pH, y = avg_SpAct, 
                            pch = factor(temp), linetype = factor(temp), fill = factor(temp)))+
  geom_errorbar(aes(ymin = avg_SpAct - sd_SpAct,
                    ymax = avg_SpAct + sd_SpAct), linetype = 1, width = 0.1)+
  geom_path()+
  geom_point(size = 3, stroke = 0.5)+
  scale_shape_manual(values = c("5" = 21, "15" = 22, "25" = 24, "35" = 25),
                     labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                                "25" = expression(25*degree*C), "35" = expression(35*degree*C))) +
  scale_fill_manual(values = c("5" = "white", "15" = "black", "25" = "white", "35" = "black"),
                    labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                               "25" = expression(25*degree*C), "35" = expression(35*degree*C)))+
  scale_linetype_manual(values = c("5" = 1, "15" = 2, "25" = 3, "35" = 5),
                        labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                                   "25" = expression(25*degree*C), "35" = expression(35*degree*C)))+
  scale_x_continuous("pH", breaks = c(3.5, 4.5, 5.5, 6.5, 7.5))+
  scale_y_continuous(expression(paste("APase specific activity, ", mu, "mol h"^-1*"mg"["enzyme"]^-1*"")),
                     limits = c(0, 35))+
  annotate(geom = "text", x = 6.7, y = 30, size = 3,
           label = c("pH: p < 0.001\ntemperature: p < 0.001\npH*temperature: p < 0.001"))+
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(),
        legend.key = element_rect(fill = NA), legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.2, 0.1), legend.spacing.x = unit(0.1, "cm"),
        axis.text = element_text(size = 12, colour = "black"), legend.key.height = unit(0.4, "cm"))+
  guides(fill = guide_legend(ncol = 2), 
         shape = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))

p2 <- ggplot(avg_apase, aes(x = 1/(temp + 273.15), y = avg_lnSpAct, 
                            pch = factor(pH), linetype = factor(pH), fill = factor(pH)))+
  geom_errorbar(aes(ymin = avg_lnSpAct - sd_lnSpAct,
                    ymax = avg_lnSpAct + sd_lnSpAct), linetype = 1, width = 0.000005)+
  geom_smooth(method = "lm", se = FALSE, color = "black", lwd = 0.5)+
  geom_point(size = 3, stroke = 0.5)+
  scale_shape_manual(values = c("3.5" = 15, "4.5" = 21, "5.5" = 19, "6.5" = 24, "7.5" = 17),
                     labels = c("3.5" = "pH 3.5", "4.5" = "pH 4.5", "5.5" = "pH 5.5", "6.5" = "pH 6.5", "7.5" = "pH 7.5")) +
  scale_fill_manual(values = c("3.5" = "black", "4.5" = "white", "5.5" = "black", "6.5" = "white", "7.5" = "black"),
                    labels = c("3.5" = "pH 3.5", "4.5" = "pH 4.5", "5.5" = "pH 5.5", "6.5" = "pH 6.5", "7.5" = "pH 7.5")) +
  scale_linetype_manual(values = c("3.5" = 6, "4.5" = 1, "5.5" = 3, "6.5" = 4, "7.5" = 5),
                        labels = c("3.5" = "pH 3.5", "4.5" = "pH 4.5", "5.5" = "pH 5.5", "6.5" = "pH 6.5", "7.5" = "pH 7.5")) +
  scale_x_continuous("1/T, K", breaks = round(unique(1/(avg_apase$temp + 273.15)), 5))+
  scale_y_continuous(expression(paste("ln(APase specific activity), ", mu, "mol h"^-1*"mg"["enzyme"]^-1*"")),
                     limits = c(0, 4))+
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(),
        legend.key = element_rect(fill = NA), legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.2, 0.1), legend.spacing.x = unit(0.1, "cm"),
        axis.text = element_text(size = 12, colour = "black"), legend.key.height = unit(0.4, "cm"))+
  guides(fill = guide_legend(ncol = 2), 
         shape = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))

save_plot("Plots/Fig1.png", plot_grid(p1, p2, ncol = 2, nrow = 1, labels = c("(a)", "(b)")), 
          dpi = 450, base_width = 12, base_height = 5)
rm(list = ls())

## Fig. 2
flowratios <- read_csv("Output/flowratios.csv")
head(flowratios)

svp2 <- ggplot(transform(flowratios, Pair = factor(Pair, levels = c("C:P", "N:P", "C:N"))), 
              aes(x = ph, y = ratio, 
                  shape = factor(tempC), fill = factor(tempC), linetype = factor(tempC)))+
  facet_grid(rows = vars(Pair), labeller = labeller(c("C:N" = "C:N estimated flow ratio",
                                                      "C:P" = "C:P estimated flow ratio", 
                                                      "N:P" = "N:P estimated flow ratio")),
             scales = "free_y", switch = "y")+
  geom_path() +
  geom_point(size = 3) +
  scale_shape_manual(values = c("5" = 21, "15" = 22, "25" = 23),
                     labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                                "25" = expression(25*degree*C)))+
  scale_fill_manual(values = c("5" = "black", "15" = "white", "25" = "black"),
                    labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                               "25" = expression(25*degree*C)))+
  scale_linetype_manual(values = c("5" = 1, "15" = 2, "25" = 3),
                        labels = c("5" = expression(5*degree*C), "15" = expression(15*degree*C),
                                   "25" = expression(25*degree*C))) +
  scale_y_continuous(limits = c(0, NA), c("Estimated flow ratio"))+
  scale_x_continuous("pH", limits = c(4.5, 7.5), breaks = c(4.5, 5.5, 6.5, 7.5))+
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(),
        legend.key = element_rect(fill = NA), legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.8, 0.05), legend.spacing.y = unit(0.05, "cm"),
        axis.text = element_text(size = 12, colour = "black"), 
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "white"), strip.placement = "outside",
        strip.text = element_text(size = 12), legend.key.height = unit(0.4, "cm"))+
  geom_text(x = 4.5, y = Inf, vjust = 1.5, aes(label = label), 
            inherit.aes = F, fontface = "bold", size = 4.5,
            data = data.frame(label = c("(a)", "(b)", "(c)"), Pair = as.factor(c("C:P", "N:P", "C:N"))))

save_plot("Plots/Fig2.png", svp2, dpi = 450, base_height = 8, base_width = 4)
rm(list = ls())

### Fig. 3
data <- read_csv("Output/arrdata.csv") 
data <- data %>% filter(ph != 3.5 & ph != 8.5)

for(enz in unique(data$enzyme)){
  dtf <- setNames(data.frame(spline(x = data$ph[data$enzyme == enz], 
                                    y = data$Ea[data$enzyme == enz], n = 31)), c("ph", "EaSpline"))
  dtf[, "EaSplineSd"] <- data.frame(spline(x = data$ph[data$enzyme == enz], 
                                           y = data$sd[data$enzyme == enz], n = 31))[, 2]
  dtf[, "enzyme"] <- enz
  assign(paste0(enz), dtf)
}
splines <- rbind(apase, bgase, nagase)

svp3 <- ggplot(splines, aes(x = ph, y = EaSpline, linetype = factor(enzyme)))+
  geom_ribbon(aes(ymin = EaSpline - EaSplineSd,
                  ymax = EaSpline + EaSplineSd), fill = "grey", alpha = 0.5,
              show.legend = F)+
  geom_line()+
  scale_linetype_manual(values = c("apase" = 3, "bgase" = 1, "nagase" = 2),
                        labels = c("apase" = "APase", "bgase" = "BGase", "nagase" = "NAGase"))+
  scale_x_continuous("pH", limits = c(4.5, 7.5), 
                     breaks = c(4.5, 5.5, 6.5, 7.5), expand = c(0.01, 0.01))+
  scale_y_continuous(expression(paste("E"["a"]*", kJ mol"^-1*"")), limits = c(-2, 75),
                     breaks = seq(0, 70, 10))+
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(),
        legend.key = element_rect(fill = NA), legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.85, 0.9), legend.spacing.y = unit(0.01, "cm"),
        legend.key.height = unit(0.4, "cm"),
        axis.text = element_text(size = 12, colour = "black"), 
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "white"), strip.placement = "outside",
        strip.text = element_text(size = 12))
save_plot("Plots/Fig3.png", svp3, dpi = 450, base_height = 4, base_width = 5)
rm(list = ls())


## Fig. S1
cvflows <- read_csv("Output/flowratios_cv.csv")
head(cvflows)

ggplot(transform(cvflows, pair = factor(pair, levels = c("C:P", "N:P", "C:N"))), 
       aes(x = ph, y = cv_fr*100, shape = pair, linetype = pair, fill = pair))+
  geom_path()+
  geom_point(size = 4, stroke = 0.5)+
  scale_shape_manual(values = c("C:P" = 22, "N:P" = 23, "C:N" = 21),
                     labels = c("C:P" = "C:P flow ratio", "N:P" = "N:P flow ratio", "C:N" = "C:N flow ratio"))+
  scale_fill_manual(values = c("C:P" = "white", "N:P" = "black", "C:N" = "white"),
                    labels = c("C:P" = "C:P flow ratio", "N:P" = "N:P flow ratio", "C:N" = "C:N flow ratio"))+
  scale_linetype_manual(values = c("C:P" = 2, "N:P" = 3, "C:N" = 1),
                        labels = c("C:P" = "C:P flow ratio", "N:P" = "N:P flow ratio", "C:N" = "C:N flow ratio"))+
  scale_x_continuous("pH", limits = c(4.5, 7.5), breaks = seq(4.5, 7.5, 1))+
  scale_y_continuous("Estimated CV with warming, %", limits = c(0, 70), breaks = seq(0, 70, 10))+
  theme(panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA), legend.title = element_blank(),
        legend.key = element_rect(fill = NA), legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.85, 0.9), legend.spacing.y = unit(0.01, "cm"),
        legend.key.height = unit(0.4, "cm"),
        axis.text = element_text(size = 12, colour = "black"), 
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "white"), strip.placement = "outside",
        strip.text = element_text(size = 12))
save_plot("Plots/FigS1.png", last_plot(), dpi = 450, base_height = 4, base_width = 5)
rm(list = ls())
