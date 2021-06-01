#####

setwd("./FlowRatiosProjections/")

library(rgdal)
library(ggplot2)
library(sp)
library(raster)
library(viridis)
library(ggpolypath)
library(cowplot)
library(scales)

world <- readOGR("Data/NaturalEarth/ne_110m_land.shp")
world_robin <- spTransform(world, CRS("+proj=robin"))
world_df <- fortify(world_robin)

bbox <- shapefile("Data/NaturalEarth/ne_110m_wgs84_bounding_box.shp") 
bbox_robin <- spTransform(bbox, CRS("+proj=robin"))
bbox_df<- fortify(bbox_robin)

graticule <- shapefile("Data/NaturalEarth/ne_110m_graticules_30.shp")
graticule_robin <- spTransform(graticule, CRS("+proj=robin"))
graticule_df <- fortify(graticule_robin)

plots_bg <- list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.y = element_blank(),
                       plot.title = element_text(size = 15),
                       plot.margin = unit(c(.5,-.5,.5,-.5), "lines"), 
                       legend.title = element_text(),
                       axis.title.x = element_blank(),
                       legend.margin=margin(t=-1, r=0, b=0, l=0, unit="lines"), 
                       legend.position="bottom", 
                       legend.direction = "horizontal",
                       legend.key = element_rect(colour = "black", size = 4),
                       legend.key.width = unit(2.0, "cm"), # 1.5
                       legend.key.height = unit(0.3, "cm"), #0.2
                       legend.text = element_text(margin = margin(t = .2, unit = "lines"))),
                 guides(colour = guide_colorbar(title.position = "top")))

depths <- c("000", "030", "100")
depths
ssp <- c(245, 585)
ratios <- c("cp", "np", "cn")
names(ratios) <- c("BG:AP", "NAG:AP", "BG:NAG")

dfratio <- read.csv("Output/relativeflowratiochangestats.csv")
dfratio

for(r in ratios){
  i <- 1
  avg <- round(mean(dfratio$mean[dfratio$ratio == r]), 0)
  sd <- round(sd(dfratio$mean[dfratio$ratio == r]), 0)
  mn <- round(min(dfratio$min[dfratio$ratio == r]), 0)
  mx <- round(max(dfratio$max[dfratio$ratio == r]), 0)
  mmn <- abs(round((abs(avg -2*sd) - abs(mn))/2, 0))
  mmx <- abs(round((abs(mx) - abs(avg+2*sd))/2, 0))
  files <- dir(paste0("Output/ratios/relchange/", r))
  for(s in ssp){
    for(d in depths){
      p <- paste0(d, s)
      pd <- files[grep(p, files)]
      rf <- raster(paste0("Output/ratios/relchange/", r, "/", pd))
      if(r == "cn"){
        lb <- c("Relative change in C:N flow ratio, %")
      } else if(r == "cp"){
        lb <- c("Relative change in C:P flow ratio, %")
      } else{
        lb <- c("Relative change in N:P flow ratio, %")
      }
      tmp <- rf
      memory.limit(size = 100000)
      tmp_p <- projectRaster(tmp, crs = "+proj=robin", over = T)
      tmp_df <- as.data.frame(tmp_p, xy = TRUE)
      colnames(tmp_df) <- c("x", "y", "data")
      gc()
      
      gg <- ggplot() + 
        geom_polygon(data = bbox_df, 
                     aes(x = long, y = lat), colour = "lightgray", 
                     fill= "lightgray", size = 0.15) +
        geom_path(data = graticule_df, aes(long, lat, group = group),
                  linetype = "dashed", color = "black", size = 0.1) +
        geom_polypath(data = world_df, aes(long, lat, group = group),
                      fill = "white", color = "white", size = 0.3) +  
        geom_raster(data = tmp_df, aes(x, y, fill = data)) +
        scale_fill_gradientn(na.value = "transparent", colors = rev(viridis(7)),
                             values = rescale(c(floor(mn), floor(mn + mmn), floor(avg - 2*sd),
                                                round(avg, 1),
                                                ceiling(avg + 2*sd), ceiling(mx - mmx), ceiling(mx))),
                             breaks = c(floor(mn), floor(avg - 2*sd),
                                        round(avg, 1),
                                        ceiling(avg + 2*sd), ceiling(mx)),
                             minor_breaks = c(floor(avg - sd), ceiling(avg + sd)),
                             limits = c(mn, mx),
                             aesthetics = "fill",
                             guide = guide_colourbar(even.steps = FALSE,
                                                     show.limits = FALSE,
                                                     title.position = "top",
                                                     title.hjust = 0.5, ticks = TRUE,
                                                     label.theme = element_text(size = 8,
                                                                                lineheight = 0.8),
                                                     label = TRUE,
                                                     frame.colour = "black")) +
        coord_equal() + 
        theme_classic(base_size = 12) + 
        plots_bg +
        labs(fill = lb)
      p <- gg + geom_polygon(data = world_df, aes(x = long, y = lat, group = group), 
                             color = "black", alpha = 0, size = 0.25) + 
        geom_polygon(data = bbox_df, aes(x = long, y = lat, group = group),
                     color = "black", alpha = 0, size = 0.25)
      
      assign(paste0("pl", i), p)
      i <- i + 1
    }
  }
  gc()
  memory.limit(size = 100000)
  prow <- plot_grid(pl1 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    pl2 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    pl3 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    pl4 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    pl5 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    pl6 + theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none", plot.title = element_blank()),
                    align = 'h',
                    labels = c("(a)", "(d)", 
                               "(b)", "(e)", 
                               "(c)", "(f)"),
                    byrow = F,
                    hjust = -3, 
                    nrow = 3,
                    ncol = 2)
  gc()
  lg <- get_legend(pl1)
  prow2 <- plot_grid(prow, lg, ncol = 1, rel_heights = c(1, .1))
  save_plot(paste0("Maps/Fig", 4 + which(ratios == r), ".png"), prow2, dpi=500, 
            base_width = 7, base_height = 4, nrow = 3, ncol = 2)
  gc()
}

