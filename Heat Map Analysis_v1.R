library("deSolve"); library("ggplot2"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast");
library("cowplot"); library("metR"); library("grid"); library("gridExtra")

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
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

#Importing in the Datasets
import <- function(id) {
  data <- data.frame(matrix(ncol = 6, nrow = 0))
  for(i in 1:5) {
    test  <- cbind(read.csv(paste0("results_ABC_SMC_gen_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

# Identify the MAP for the Parameter Sets ---------------------------------
#Import of Posterior Distributions
data <- do.call(rbind, list(import(tet), import(amp), import(broil)))

MAPtet <- c("phi" <- mean(data$phi[which(data$group == "data5" & data$fit == "tet")]),
            "theta" <- mean(data$theta[which(data$group == "data5" & data$fit == "tet")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data5" & data$fit == "tet")]),
            "alpha" <- mean(data$alpha[which(data$group == "data5" & data$fit == "tet")]),
            "zeta" <- mean(data$zeta[which(data$group == "data5" & data$fit == "tet")]))

MAPamp <- c("phi" <- mean(data$phi[which(data$group == "data5" & data$fit == "amp")]),
            "theta" <- mean(data$theta[which(data$group == "data5" & data$fit == "amp")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data5" & data$fit == "amp")]),
            "alpha" <- mean(data$alpha[which(data$group == "data5" & data$fit == "amp")]),
            "zeta" <- mean(data$zeta[which(data$group == "data5" & data$fit == "amp")]))

MAPbroil <- c("phi" <- mean(data$phi[which(data$group == "data5" & data$fit == "broil")]),
              "theta" <- mean(data$theta[which(data$group == "data5" & data$fit == "broil")]),
              "betaAA" <- mean(data$betaAA[which(data$group == "data5" & data$fit == "broil")]),
              "alpha" <- mean(data$alpha[which(data$group == "data5" & data$fit == "broil")]),
              "zeta" <- mean(data$zeta[which(data$group == "data5" & data$fit == "broil")]))

MAP <- rbind(MAPtet, MAPamp, MAPbroil); colnames(MAP) <- c("phi", "theta", "betaAA", "alpha", "zeta")

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

#These Parameters Are Based on MAP from Model Fitting

parmstet_pigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPtet","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPtet","phi"] , 
                  theta = MAP["MAPtet","theta"], alpha = MAP["MAPtet","alpha"] , tau = 0, zeta = MAP["MAPtet","zeta"])

parmstet_broil = c(ra = 0, rh = (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP["MAPbroil","betaAA"], 
                   betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPbroil","phi"], 
                   theta = MAP["MAPbroil","theta"], alpha = MAP["MAPbroil","alpha"] , tau = 0, 
                   zeta = MAP["MAPbroil","zeta"])

parms_amppigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPamp","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPamp","phi"], 
                  theta = MAP["MAPamp","theta"], alpha = MAP["MAPamp","alpha"], tau = 0, zeta = MAP["MAPamp","zeta"])

heatmap <- list()

for(j in 1:3) {
  parms = list(parmstet_pigs, parmstet_broil, parms_amppigs)[[j]]
  
  heatmap[[j]] = local({
    
    parametertest <- list()
    
    parms1 <- parms
    i = 0
    
    for(z in 1:3) {
      
      if(z == 1) {
        parameterspace <- expand.grid("betaAA" = seq(parms[["betaAA"]]*0.75, parms[["betaAA"]], by = (parms[["betaAA"]]*0.25)/25), 
                                      "betaHA" = seq(0.00001*0.75, 0.00001, by = 0.00001*0.25/25))  
      }
      if(z == 2) {
        parameterspace <- expand.grid("zeta" = seq(parms[["zeta"]]*0.75, parms[["zeta"]], by = (parms[["zeta"]]*0.25)/25), 
                                      "betaHA" = seq(0.00001*0.75, 0.00001, by = 0.00001*0.25/25))
      }
      
      if(z == 3) {
        parameterspace <- expand.grid("both" = seq(0.75, 1, by = 0.01), 
                                      "betaHA" = seq(0.00001*0.75, 0.00001, by = 0.00001*0.25/25))
      }
      
      scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 9))
      
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
        
        out <- data.frame(ode(y = init, func = amr, times = times, parms = parms1))
        scendata[i,] <- c("icombh" = (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000 , 
                          "resrat" =  rounding(out[nrow(out),7])/(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])), 
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
        
        print(paste0(c("parmstet_pigs 1", "parmstet_broil 2", "parms_amppigs 3")[j], " - " ,
                     c("beta 1", "zeta 2", "both 3")[z], " - ",round(i/nrow(parameterspace), digit = 2)*100, "%"))
      }
      colnames(scendata) <- c("icombh", "resrat", "betaAA","percbetaAA", "betaHA", "percbetaHA", "zeta", "perczeta", "percdecrease")
      parametertest[[z]] <- scendata
    }
    return(parametertest)
  })
}

# Plotting ----------------------------------------------------------------

plotheat <- list()

for(i in 1:3) {
  
  plottemp <- list()
  
  for(j in 1:3) {
    scentest <- heatmap[[i]][[j]]
    
    breaks <- c(0, 3.26, seq(3.4, max(scentest$icombh)+0.1, by = 0.1))
    
    if(j == 1) {
      plot <- ggplot(scentest, aes(percbetaAA, percbetaHA, z = icombh))
      if(i == 1) {
        plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", 
                            title = paste(""))
      }
      if(i == 2) {
        plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", 
                            title = paste(""))
      }
      if(i == 3) {
        plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", 
                            title = paste(""))
      }
    } 
    if(j == 2) {
      plot <- ggplot(scentest, aes(perczeta, percbetaHA, z = icombh)) + 
        labs(x = bquote("% of Baseline"~zeta), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", title = "")
    } 
    if(j == 3) {
      plot <- ggplot(scentest, aes(percdecrease, percbetaHA, z = icombh)) + 
        labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", title = "")
    } 
    
    plot <- plot + metR::geom_contour_fill(breaks = breaks, color = "black", size = 0.1)  + 
      geom_contour(color = "red", size = 1, breaks = 3.26, alpha = 0.8) +
      metR::geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks, label.placement = label_placement_fraction(frac = 0.5),
                              stroke = 0.05, stroke.color = "black",) +
      scale_fill_viridis_b(breaks = breaks, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
            plot.title = element_text(size = 18, vjust = 3, hjust = 0.1, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.75, "cm"),
            legend.key.width =  unit(2, "cm"))
    plottemp[[j]] <- plot
  }
  
  combplot <- ggarrange(plottemp[[1]], plottemp[[2]], plottemp[[3]], ncol = 3, nrow = 1, common.legend = TRUE,
                        legend = "bottom", labels = c("A", "B", "C")[i], font.label = list(size = 25), vjust = 1.2)
  title1=text_grob(c("Tetracycline Sales in Fattening Pigs", "Tetracycline Sales in Broiler Poultry",
                     "Ampicillin Sales in Fattening Pigs")[i], size = 18, face = "bold", vjust = 2)
  
  combplot <- grid.arrange(combplot, top = title1)
  plotheat[[i]] <- combplot
}

# Collating the Plots -----------------------------------------------------

combplot <- ggarrange(plotheat[[1]], plotheat[[2]], plotheat[[3]], ncol = 1, nrow = 3)
ggsave(combplot, filename = "HeatMapcomb_poster.png", dpi = 300, type = "cairo", width = 11, height = 13, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

ggsave(combplot, filename = "HeatMapcomb.png", dpi = 300, type = "cairo", width = 12, height = 14, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")