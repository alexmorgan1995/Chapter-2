library("deSolve"); library("ggplot2"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast");
library("cowplot"); library("metR")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData")

# Model Functions ----------------------------------------------------------

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model ODEs
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa  
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

#These Parameters Are Based on MAP from Model Fitting

parmstet_pigs = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.074716), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.010948457, theta = 0.008345866, alpha = 0.28247322, tau = 0)
parmstet_broil = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.074716), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.010948457, theta = 0.008345866, alpha = 0.28247322, tau = 0)
parms_amppigs = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.074716), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.010948457, theta = 0.008345866, alpha = 0.28247322, tau = 0)

heatmap <- list()

for(j in 1:3) {
  parms = c(parmstet_pigs, parmstet_broil, parms_amppigs)[j]
  
  parameterspace <- expand.grid("betaAA" = seq(parms[["betaAA"]]*0.75, parms[["betaAA"]], by = (parms[["betaAA"]]*0.25)/50), 
                                "betaHA" = seq(0.00001*0.75, 0.00001, by = 0.00001*0.25/50))
  
  heatmap[[j]] = local({
    
    i = 0
    
    scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 6))
    
    for(i in 1:nrow(parameterspace)) {
      
      print(paste0("HeatMap - ", c("ICombH", "ResRate")[j]," | ", round(i/nrow(parameterspace), digits = 3)*100, "%"))
      parms["betaAA"] <- parameterspace[i,1]
      parms["betaHA"] <- parameterspace[i,2] 
      
      out <- data.frame(ode(y = init, func = amr, times = times, parms = parms))
      
      scendata[i,] <- c("icombh" = (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000 , 
                        "resrat" =  rounding(out[nrow(out),7])/(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])), 
                        "betaAA" = parms[["betaAA"]], 
                        "percbetaAA" = (parms[["betaAA"]]/ 0.074716)*100,
                        "betaHA" = parms[["betaHA"]],
                        "percbetaHA" = (parms[["betaHA"]]/ 0.00001)*100)
    }
    
    colnames(scendata) <- c("icombh", "resrat", "betaAA","percbetaAA", "betaHA", "percbetaHA")
    
    return(scendata)
  })
}


test <- heatmap[[1]]
#test$icombh[test$icombh <= 3.26] <- 3.26


breaks <- c(0, 3.26, seq(3.4, 4, by = 0.1))

heatmaptetpig <- ggplot(test, aes(percbetaAA, percbetaHA, z = icombh)) +
  metR::geom_contour_fill(breaks = breaks, color = "black", size = 0.1)  + 
  geom_contour(color = "red", size = 1, breaks = 3.26, alpha = 0.8) +
  metR::geom_text_contour(col = "white", fontface = "bold", size = 5, breaks = breaks) +
  scale_fill_viridis_b(breaks = breaks, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.7,1, by = 0.3/8))) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  theme(legend.position = "right", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 18, vjust = 3, hjust = 0.5, face = "bold"),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(1, "cm"),
        legend.key.width =  unit(0.5, "cm")) + 
  labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", 
       title = paste("Tetracycline Usage in Fattening Pigs"))


ggsave(heatmaptetpig, filename = "HeatMap_tetpigs.png", dpi = 300, type = "cairo", width = 8, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft Figures")
