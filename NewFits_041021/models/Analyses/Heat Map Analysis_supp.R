library("deSolve"); library("ggplot2"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("metR"); 
library("grid"); library("gridExtra"); library("rootSolve"); library("metR")

rm(list=ls())
setwd("/Users/amorgan/Documents/PhD_Work/Chapter-2/NewFits_041021/data/new/full")

# Model Functions ----------------------------------------------------------

#Model ODEs
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      (0.5*zeta)*Sa*(1-alpha) - (0.5*zeta)*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + (0.5*zeta)*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + (0.5*zeta)*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    
    CumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    CumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), CumS, CumR))
  })
}

#Importing in the Datasets
import <- function(id) {
  data <- data.frame(matrix(ncol = 6, nrow = 0))
  for(i in 1:length(grep(paste0("post_", id), list.files(), value = TRUE))) {
    test  <- cbind(read.csv(paste0("ABC_post_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

# Identify the MAP for the Parameter Sets ---------------------------------
#Import of Posterior Distributions

data <- list(import("tetbroil"), import("ampbroil"), import("tetpigs"), import("amppigs"))
lapply(1:4, function(x) data[[x]]$group = factor(data[[x]]$group, levels = unique(data[[x]]$group)))

#Obtain the MAPs for each dataset

MAP <- rbind(c(colMeans(data[[2]][which(data[[2]]$group == tail(unique(data[[2]]$group),1)),][,1:6])),
             c(colMeans(data[[1]][which(data[[1]]$group == tail(unique(data[[1]]$group),1)),][,1:6])),
             c(colMeans(data[[4]][which(data[[4]]$group == tail(unique(data[[4]]$group),1)),][,1:6])),
             c(colMeans(data[[3]][which(data[[3]]$group == tail(unique(data[[3]]$group),1)),][,1:6])))

colnames(MAP) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
rownames(MAP) <- c("ampbroil", "tetbroil","amppigs", "tetpigs")

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

#These Parameters Are Based on MAP from Model Fitting
parms_ampbroil = c(ra = 0, rh = (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP["ampbroil","betaAA"], 
                   betaAH = 0.00001, betaHH = 0.00001, betaHA = MAP["ampbroil","betaHA"], phi = MAP["ampbroil","phi"] , 
                   kappa = MAP["ampbroil","kappa"], alpha = MAP["ampbroil","alpha"] , tau = 0, zeta = MAP["ampbroil","zeta"])

parms_tetbroil = c(ra = 0, rh = (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP["tetbroil","betaAA"], 
                   betaAH = 0.00001, betaHH = 0.00001, betaHA = MAP["tetbroil","betaHA"], phi = MAP["tetbroil","phi"], 
                   kappa = MAP["tetbroil","kappa"], alpha = MAP["tetbroil","alpha"] , tau = 0, 
                   zeta = MAP["tetbroil","zeta"])

parms_amppigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["amppigs","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = MAP["amppigs","betaHA"], phi = MAP["amppigs","phi"], 
                  kappa = MAP["amppigs","kappa"], alpha = MAP["amppigs","alpha"], tau = 0, zeta = MAP["amppigs","zeta"])

parms_tetpigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["tetpigs","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = MAP["tetpigs","betaHA"], phi = MAP["tetpigs","phi"], 
                  kappa = MAP["tetpigs","kappa"], alpha = MAP["tetpigs","alpha"], tau = 0, zeta = MAP["tetpigs","zeta"])

heatmap <- list()

for(j in 1:4) {
  parms = list(parms_ampbroil, parms_tetbroil, parms_amppigs,parms_tetpigs)[[j]]
  
  heatmap[[j]] = local({
    
    parametertest <- list()
    
    parms1 <- parms
    i = 0
    
    for(z in 1:3) {
      
      if(z == 1) {
        parameterspace <- expand.grid("betaAA" = seq(parms[["betaAA"]]*0, parms[["betaAA"]], by = parms[["betaAA"]]/25), 
                                      "betaHA" = seq(parms[["betaHA"]]*0.75, parms[["betaHA"]], by = parms[["betaHA"]]*0.25/25))  
      }
      if(z == 2) {
        parameterspace <- expand.grid("zeta" = seq(parms[["zeta"]]*0, parms[["zeta"]], by = parms[["zeta"]]/25), 
                                      "betaHA" = seq(parms[["betaHA"]]*0.75, parms[["betaHA"]], by = parms[["betaHA"]]*0.25/25))
      }
      if(z == 3) {
        parameterspace <- expand.grid("both" = seq(0, 1, by = 1/25), 
                                      "betaHA" = seq(parms[["betaHA"]]*0.75, parms[["betaHA"]], by = parms[["betaHA"]]*0.25/25))
      }
      
      scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 10))
      
      for(i in 1:nrow(parameterspace)) {
        
        #print(paste0("HeatMap - ", c("tet_pigs", "tet_broil", "amp_pigs")[j]," | ", round(i/nrow(parameterspace), digits = 3)*100, "%"))
        parms1["betaHA"] <- parameterspace[i,2] 
        
        if(z == 1) {
          parms1["betaAA"] <- parameterspace[i,1]
        }
        if(z == 2) {
          parms1["zeta"] <- parameterspace[i,1]
        }
        if(z == 3) {
          parms1["betaAA"] <- parms["betaAA"]*parameterspace[i,1]
          parms1["zeta"] <- parms["zeta"]*parameterspace[i,1]
        }
        
        out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms1)
        scendata[i,] <- c("icombh" = ((out[[2]] + out[[3]])*(446000000))/100000, 
                          "prevH" =  out[[1]][["Ish"]] + out[[1]][["Irh"]], 
                          "resrat" =  out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]), 
                          "betaAA" = parms[["betaAA"]], 
                          "percbetaAA" = (parms1[["betaAA"]]/ parms[["betaAA"]])*100,
                          "betaHA" = parms1[["betaHA"]],
                          "percbetaHA" = (parms1[["betaHA"]]/ parms[["betaHA"]])*100,
                          "zeta" = parms1[["zeta"]],
                          "perczeta" = (parms1[["zeta"]]/ parms[["zeta"]])*100,
                          "percdecrease_both" = parameterspace[i,1]*100)
        if(z == 1 | 2) {
          scendata[["percdecrease_both"]] == 1 
        }
        
        print(paste0(c("ampbroil", "tetbroil", "amppigs", "tetpigs")[j], " - " ,
                     c("beta 1", "zeta 2", "both 3")[z], " - ",round(i/nrow(parameterspace), digit = 2)*100, "%"))
      }
      colnames(scendata) <- c("icombh","prevH", "resrat", "betaAA","percbetaAA", "betaHA", "percbetaHA", "zeta", "perczeta", "percdecrease")
      parametertest[[z]] <- scendata
    }
    return(parametertest)
  })
}


#Lone Plots

scentest1 <- heatmap[[1]][[3]]
scentest2 <- heatmap[[2]][[3]]
scentest3 <- heatmap[[3]][[3]]
scentest4 <- heatmap[[4]][[3]]

breaks1 <- c(0, seq(0.593, max(scentest1$icombh)+0.05, by = 0.02))

plot1 <- ggplot(scentest1, aes(percdecrease, percbetaHA, z = icombh)) + metR::geom_contour_fill(breaks = breaks1, color = "black", size = 0.1)  + 
  labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), fill = expression(paste("Daily \nIncidence")), title = "Ampicillin Usage in Broiler Poultry") +
  geom_contour(color = "red", size = 1, breaks = 0.593, alpha = 0.8) +
  geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks1, label.placer = label_placer_fraction(frac = 0.5),
                          stroke = 0.05, stroke.color = "black",) +
  scale_fill_viridis_b(breaks = breaks1, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
  scale_y_continuous(expand = c(0,0), limits = c(75, 100)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + theme_bw() +
  theme(legend.position = "right", legend.title = element_text(size=12), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
        legend.key.width =  unit(0.5, "cm"))

breaks2 <- c(0,  seq(0.593, max(scentest2$icombh)+0.05, by = 0.015))

plot2 <- ggplot(scentest2, aes(percdecrease, percbetaHA, z = icombh)) + metR::geom_contour_fill(breaks = breaks2, color = "black", size = 0.1)  + 
  labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), fill = expression(paste("Daily \nIncidence")), 
       title = "Tetracycline Usage in Broiler Poultry") +
  geom_contour(color = "red", size = 1, breaks = 0.593, alpha = 0.8) +
  geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks2, label.placer = label_placer_fraction(frac = 0.5),
                          stroke = 0.05, stroke.color = "black",) +
  scale_fill_viridis_b(breaks = breaks2, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
  scale_y_continuous(expand = c(0,0), limits = c(75, 100)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + theme_bw() +
  theme(legend.position = "right", legend.title = element_text(size=12), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
        legend.key.width =  unit(0.5, "cm"))

breaks3 <- c(0,  seq(0.593, max(scentest3$icombh)+0.05, by = 0.025))

plot3 <- ggplot(scentest3, aes(percdecrease, percbetaHA, z = icombh)) + metR::geom_contour_fill(breaks = breaks3, color = "black", size = 0.1)  + 
  labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), fill = expression(paste("Daily \nIncidence")), title = "Ampicillin Usage in Fattening Pigs") +
  geom_contour(color = "red", size = 1, breaks = 0.593, alpha = 0.8) +
  geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks3, label.placer = label_placer_fraction(frac = 0.5),
                          stroke = 0.05, stroke.color = "black",) +
  scale_fill_viridis_b(breaks = breaks3, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
  scale_y_continuous(expand = c(0,0), limits = c(75, 100)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + theme_bw() +
  theme(legend.position = "right", legend.title = element_text(size=12), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
        legend.key.width =  unit(0.5, "cm"))

breaks4 <- c(0,  seq(0.593, max(scentest4$icombh)+0.05, by = 0.03))

plot4 <- ggplot(scentest4, aes(percdecrease, percbetaHA, z = icombh)) + metR::geom_contour_fill(breaks = breaks4, color = "black", size = 0.1)  + 
  labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), 
       fill = expression(paste("Daily \nIncidence")), title = "Tetracycline Usage in Fattening Pigs") +
  geom_contour(color = "red", size = 1, breaks = 0.593, alpha = 0.8) +
  geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks3, label.placer = label_placer_fraction(frac = 0.5),
                          stroke = 0.05, stroke.color = "black",) +
  scale_fill_viridis_b(breaks = breaks4, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
  scale_y_continuous(expand = c(0,0), limits = c(75, 100)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + theme_bw() +
  theme(legend.position = "right", legend.title = element_text(size=12), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.7, "cm"),
        legend.key.width =  unit(0.5, "cm"))

plot1 <- plot1 + theme(legend.position = "bottom", legend.key.height =unit(0.55, "cm"), legend.key.width =  unit(1.5, "cm"))
plot2 <- plot2 + theme(legend.position = "bottom", legend.key.height =unit(0.55, "cm"), legend.key.width =  unit(1.5, "cm"))
plot3 <- plot3 + theme(legend.position = "bottom", legend.key.height =unit(0.55, "cm"), legend.key.width =  unit(1.5, "cm"))
plot4 <- plot4 + theme(legend.position = "bottom", legend.key.height =unit(0.55, "cm"), legend.key.width =  unit(1.5, "cm"))
combplot_solo4by4 <- ggarrange(plot1, plot2, plot3,  plot4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), font.label = list(size = 20), vjust = 1.2)  +
  bgcolor("white")

ggsave(combplot_solo4by4, filename = "HeatMapcomb_solo_4x4_supplementary.png", dpi = 500, width = 10, height = 10, units = "in",
       path = "/Users/amorgan/Documents/PhD_Work/Chapter-2/NewFits_041021/figures")
