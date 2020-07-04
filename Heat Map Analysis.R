library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast");
library("cowplot")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData")

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

#These Parameters Are Based on MAP from Model Fitting

parameterspace <- expand.grid("betaAA" = seq(0,0, by =5), "betaHA" = seq(1,200, by =5))

parms = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.074716), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.010948457, theta = 0.008345866, alpha = 0.28247322)

heatmap <- list()

for(j in 1:2) {
  
  heatmap[[j]] = local({
    
    i = 0
    
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 5))
    
    for(i in 1:nrow(parameterspace)) {
      
      print(paste0("Scenario ", j," - ", round(i/nrow(parameterspace), digits = 2)))
      parms["betaAA"] <- parameterspace[i,1]
      parms["betaHA"] <- parameterspace[i,2] 
      
      parms["scen"] <- j
      
      out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
      
      scendata[i,] <- c("icombh" = max(out$I), "cum" = max(out$C), "scen" = parms[["scen"]], 
                        "tstart" = parms[["tstart"]], "t_dur" = parms[["t_dur"]])
    }
    
    colnames(scendata) <- c("peak", "cum", "scen", "tstart", "t_dur")
    
    p1 <- ggplot(scendata, aes(x = tstart, y = t_dur, fill= peak))  + geom_tile()  +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = "Intervention Trigger", y = "Intervention Duration", fill = "Peak I(t)", title = paste("Scenario", j)) + 
      scale_fill_viridis_c(direction = -1)
    
    p2<- ggplot(scendata, aes(x = tstart, y = t_dur, fill = cum))  + geom_tile() +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=15), legend.text=element_text(size=15),  axis.text=element_text(size=15),
            axis.title.y=element_text(size=15),axis.title.x = element_text(size=15),  plot.title = element_text(size = 20, vjust = 3, hjust = -0.2),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
            legend.key.width =  unit(0.5, "cm")) + 
      labs(x = "Intervention Trigger", y = "Intervention Duration", fill = "Cumulative\nIncidence", title = "") + 
      scale_fill_viridis_c(direction = -1, option = "magma") 
    
    combplot <- ggarrange(p1,p2, ncol = 2, nrow = 1, widths = c(1,1.05), align = "h")
    print(combplot)
    return(combplot)
  })
}

ggsave(combplot, filename = "Heat_5_scenarios_sensitivity.png", dpi = 300, type = "cairo", width = 10, height = 16, units = "in")
