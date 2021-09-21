library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")

# Model Functions ----------------------------------------------------------

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model ODEs
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
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
  for(i in 1:10) {
    test  <- cbind(read.csv(paste0("results_ABC_SMC_gen_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

# Data Import -------------------------------------------------------------
#Import of Posterior Distributions
data <- do.call(rbind, list(import(tet), import(amp), import(broil)))

data$group <- factor(data$group, levels = c("data1", "data2","data3","data4","data5","data6","data7","data8","data9","data10"))

#Import of Fitting Data
dataamp <- read.csv("resistanceprofAnim_amp.csv")
dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))
dataamp <- dataamp[!dataamp$N < 10,]

datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetra$lower <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[1]]))
datatetra$upper <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[2]]))
datatetra <- datatetra[!datatetra$N < 10,]

databroil <- read.csv("salm_broilers_2018.csv")
databroil$lower <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[1]]))
databroil$upper <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[2]]))
databroil <- databroil[!databroil$N < 10,]

#Obtain the MAPs for each dataset

MAPtet <- c("phi" <- mean(data$phi[which(data$group == "data10" & data$fit == "tet")]),
            "kappa" <- mean(data$kappa[which(data$group == "data10" & data$fit == "tet")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "tet")]),
            "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "tet")]),
            "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "tet")]))

MAPamp <- c("phi" <- mean(data$phi[which(data$group == "data10" & data$fit == "amp")]),
            "kappa" <- mean(data$kappa[which(data$group == "data10" & data$fit == "amp")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "amp")]),
            "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "amp")]),
            "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "amp")]))

MAPbroil <- c("phi" <- mean(data$phi[which(data$group == "data10" & data$fit == "broil")]),
              "kappa" <- mean(data$kappa[which(data$group == "data10" & data$fit == "broil")]),
              "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "broil")]),
              "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "broil")]),
              "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "broil")]))

MAP <- rbind(MAPtet, MAPamp, MAPbroil); colnames(MAP) <- c("phi", "kappa", "betaAA", "alpha", "zeta")

# ABC-SMC Posterior -------------------------------------------------------

plotlist1 <- list()

for(j in 1:length(unique(data$fit))) {
  plotlist2 <- list()
  
  for (i in 1:length(colnames(MAP))) { # Loop over loop.vector
    
    dens <- density(data[,colnames(MAP)[i]][data$fit == unique(data$fit)[j] & data$group == "data10"])
    dataplot <- melt(data[data$fit == unique(data$fit)[j],], id.vars = "group", measure.vars = colnames(MAP)[i])
    
    plotlist2[[i]] <- local({
      i = i
      p1 <- ggplot(data[data$fit == unique(data$fit)[j],], aes(x=get(colnames(MAP)[i]), fill=group)) + geom_density(alpha=.5) + theme_bw()  +
        scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5",
                                       "Generation 6", "Generation 7", "Generation 8", "Generation 9", "Generation 10"))+
        theme(legend.text=element_text(size=10), axis.text.x=element_text(size=10),axis.ticks.y=element_blank(), axis.text.y=element_blank(),
              axis.title.y=element_text(size=10), axis.title.x= element_text(size=10), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
              plot.title = element_text(size = 12, vjust = 3, hjust = 0.5, face = "bold"))
      if(colnames(MAP)[i] == "phi") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.05), expand = c(0, 0), name = expression(paste("Rate of Resistance Reversion (", phi, ")"))) +
        labs(fill = NULL, title = c("Tetracycline Sales in Fattening Pigs", "Ampicillin Sales in Fattening Pigs", "Tetracycline Sales in Boiler Poultry")[j]) + 
          scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "kappa") {
        p1 <- p1 + scale_x_continuous(limits = c(0,4),expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Recovery (", kappa, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, 1), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "betaAA") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.2),expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")))+
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      if(colnames(MAP)[i] == "alpha") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.8),expand = c(0, 0), name = expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      if(colnames(MAP)[i] == "zeta") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.15),expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      
      p1 <- p1 + geom_vline(xintercept = mean(dataplot$value[dataplot$group == "data10"]), size  = 1.2, color = "red", alpha = 0.5)
      
      return(p1)
      
    })
    
    print(paste0("Plot Parameter: ",colnames(MAP)[i], " | Data: ", unique(data$fit)[j] ))
  }
  plotlist1[[unique(data$fit)[j]]] <- plotlist2
}

abc <- ggarrange(plotlist1[[1]][[1]], plotlist1[[2]][[1]], plotlist1[[3]][[1]],
                 plotlist1[[1]][[2]], plotlist1[[2]][[2]], plotlist1[[3]][[2]],
                 plotlist1[[1]][[3]], plotlist1[[2]][[3]], plotlist1[[3]][[3]],
                 plotlist1[[1]][[4]], plotlist1[[2]][[4]], plotlist1[[3]][[4]],
                 plotlist1[[1]][[5]], plotlist1[[2]][[5]], plotlist1[[3]][[5]],
                 nrow = 5, ncol =3, 
                 labels =  c("A","B","C",
                             "","","",
                             "","","",
                             "","","",
                             "","",""),
                 font.label = c(size = 20), common.legend = TRUE, legend = "bottom",
                 align = "hv", vjust = 1.05)

ggsave(abc, filename = "ABC_SMC_Post_v1.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# Model Fit with Data - Ribbon -----------------------------------------------------

start_time <- Sys.time()

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
ribbon_final <- data.frame()
times <- seq(0, 200000, by = 100)

for(j in 1:nrow(MAP)) {
  output1 <- data.frame()
  output_ribbon <- data.frame()
  ribbondata <- data[data$group == "data10" & data$fit == unique(data$fit)[j],] 
  
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=3))
    temp_ribbon <- data.frame(matrix(NA, nrow = 0, ncol=4))
    
    if(rownames(MAP)[j] == "MAPbroil") {
      parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])
    } 
    else {
      parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])
    }
    
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[i]
    temp[1,2] <- rounding(out[nrow(out),5]) 
    temp[1,3] <- rounding(out[nrow(out),6]) 
    temp[1,4] <- rounding(out[nrow(out),7])
    temp[1,5] <- temp[1,3] + temp[1,4]
    temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
    temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[1,8] <- rownames(MAP)[j]
    
    for(z in 1:1000) {
      temp_ribbon_tau <- data.frame(matrix(NA, nrow = 1, ncol=4))
      
      if(rownames(MAP)[j] == "MAPbroil") {
        parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])
      } 
      else {
        parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi = ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])
      }
      
      out <- ode(y = init, func = amr, times = times, parms = parms2)
      temp_ribbon_tau[1,1] <- parmtau[i]
      temp_ribbon_tau[1,2] <- z
      temp_ribbon_tau[1,3] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
      temp_ribbon_tau[1,4] <- rownames(MAP)[j]
      
      temp_ribbon <- rbind.data.frame(temp_ribbon, temp_ribbon_tau)
      
      print(paste0(temp_ribbon_tau[1,4], ", tau: ", temp_ribbon_tau[1,1], ", ", (z/1000)*100, "%"))
      
    }
    output1 <- rbind.data.frame(output1, temp)
    output_ribbon <- rbind.data.frame(output_ribbon, temp_ribbon)
  }
  icombhdata <- rbind(icombhdata, output1)
  ribbon_final <- rbind(ribbon_final, output_ribbon)
}

colnames(icombhdata)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
icombhdata[,2:5] <- icombhdata[,2:5]*100000 #Scaling the prevalence (per 100,000)

colnames(ribbon_final)[1:4] <- c("tau","particle","IResRatA", "group")

HDI_ribbon <- data.frame()

for(j in 1:length(unique(ribbon_final$group))) {
  for(i in 1:length(unique(ribbon_final$tau))) {
    HDI_ribbon <- rbind(HDI_ribbon, 
                        data.frame("tau" = unique(ribbon_final$tau)[i],
                                   "lowHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[2]],
                                   "highHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[3]],
                                   "scen" = unique(ribbon_final$group)[j]))
  }
}

end_time <- Sys.time(); end_time - start_time

amp <- ggplot(dataamp, aes(x = pig_amp_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.03),expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "MAPamp",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  labs(x ="Fattening Pig Ampicillin Sales (g/PCU)", y = "Ampicillin-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPamp",], aes(x = tau, y= IResRatA), col = "red", size = 1.1) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

tet <- ggplot(datatetra, aes(x = pig_tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.034), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Fattening Pig Carriage") +
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "MAPtet",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPtet",], aes(x = tau, y= IResRatA), col = "red", size = 1.1) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

broil <- ggplot(databroil, aes(x = tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.022), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Broiler Poultry Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Broiler Poultry Carriage") +
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "MAPbroil",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.1)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

fits <- ggarrange(tet, amp, broil, nrow = 1, ncol = 3, align = "h", labels = c("A","B", "C"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(fits, filename = "Model_Fits_v1.png", dpi = 300, type = "cairo", width = 14, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# Tau and ICombH Base Plot - Requires you to run previous section ----------------------------------------------------------

icombhlist <- list()

#tet, amp, broil
averagesales <- c(0.0122887, 0.01156391, 0.006666697)

for(i in 1: length(unique(icombhdata$group))){
  i = i 
  
  icombhlist[[i]] <- local({
    
    plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[i],],
                     id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 
    
    p1 <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
      geom_vline(xintercept = averagesales[i], alpha = 0.3, size = 2) + 
      geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
      scale_y_continuous(limits = c(0,5), expand = c(0, 0))  + 
      geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[i]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
                position = "stack", angle = 45) +
      theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm')) + 
      scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) 
    
    if(unique(icombhdata$group)[i] == "MAPtet") {
      p1 <- p1 + labs(x ="Tetracycline Usage in Fattening Pig (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    if(unique(icombhdata$group)[i] == "MAPamp") {
      p1 <- p1 + labs(x ="Ampicillin Usage in Fattening Pig (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    if(unique(icombhdata$group)[i] == "MAPbroil") {
      p1 <- p1 + labs(x ="Tetracycline Usage in Broiler Poultry (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    
    return(p1)
  })
}

icombh <- ggarrange(icombhlist[[1]], icombhlist[[2]], icombhlist[[3]], nrow = 3, ncol = 1, align = "v", labels = c("A","B", "C"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(icombh, filename = "Icombh.png", dpi = 300, type = "cairo", width = 8, height = 11, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

ggsave(icombh, filename = "Icombh_poster.png", dpi = 300, type = "cairo", width = 10, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# HDI for parameter estimates ---------------------------------------------

#Data is defined at the beginning 

HDI_parms <- list()

for(i in 1:length(unique(data$fit))) {
  temp_data <- data[data$group == "data10" & data$fit == c("tet", "amp", "broil")[i],]
  HDI_parms[[i]] <- data.frame("parm" = c("phi", "kappa", "betaAA", "alpha", "zeta"),
                               "mean_parm" = unlist(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                                temp_data[["alpha"]], temp_data[["zeta"]]), mean)),
                               "lowHDI" = sapply(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                             temp_data[["alpha"]], temp_data[["zeta"]]), hdi), "[[", 2),
                               "highHDI" = sapply(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                              temp_data[["alpha"]], temp_data[["zeta"]]), hdi), "[[", 3))
}; HDI_parms

# Distances for Fits ------------------------------------------------------

summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, parms, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol=4)
  for (i in 1:length(tau_range)) {
    parms1 <- parms
    parms1[["tau"]] <- tau_range[[i]]
    temp <- matrix(NA, nrow = 1, ncol=4)
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms1)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput)))
}


datatetra <- read.csv("resistanceprofAnim_v1.csv"); datatetra <- datatetra[!datatetra$N < 10,]
colnames(datatetra)[12] <- "usage"; datatetra$usage <- datatetra$usage/1000 
databroil <- read.csv("salm_broilers_2018.csv"); databroil <- databroil[!databroil$N < 10,]
colnames(databroil)[9] <- "usage"; databroil$usage <- databroil$usage/1000 
dataamp <- read.csv("resistanceprofAnim_amp.csv"); dataamp <- dataamp[!dataamp$N < 10,]
colnames(dataamp)[9] <- "usage"; dataamp$usage <- dataamp$usage/1000 

parmstet_pigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPtet","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPtet","phi"] , 
                  kappa = MAP["MAPtet","kappa"], alpha = MAP["MAPtet","alpha"] , tau = 0, zeta = MAP["MAPtet","zeta"])

parmstet_broil = c(ra = 0, rh = (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP["MAPbroil","betaAA"], 
                   betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPbroil","phi"], 
                   kappa = MAP["MAPbroil","kappa"], alpha = MAP["MAPbroil","alpha"] , tau = 0, 
                   zeta = MAP["MAPbroil","zeta"])

parms_amppigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPamp","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPamp","phi"], 
                  kappa = MAP["MAPamp","kappa"], alpha = MAP["MAPamp","alpha"], tau = 0, zeta = MAP["MAPamp","zeta"])

for(i in 1:3) {
  casestudy <- list(datatetra, databroil, dataamp)[[i]]
  parms <- list(parmstet_pigs, parmstet_broil, parms_amppigs)[[i]]
  #print(parms)
  test<- computeDistanceABC_ALEX(sum.stats = summarystatprev, 
                          distanceABC = sum_square_diff_dist, 
                          fitmodel = amr, 
                          tau_range = casestudy$usage, 
                          parms = parms,
                          init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                          times = seq(0, 2000, by = 100), 
                          data = casestudy)
  print(test)
}

#Identify the Exact Values of the Fit
tauoutput <- matrix(nrow = 0, ncol=4)

averagesales <- c(0.0122887, 0.01156391, 0.006666697)

exploredtau <- c(datatetra$usage, averagesales[1], 0)

init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times = seq(0, 2000, by = 100)

for (i in 1:length(exploredtau)) {
  parms1 <- parmstet_pigs
  parms1[["tau"]] <- exploredtau[[i]]
  temp <- matrix(NA, nrow = 1, ncol=4)
  out <- ode(y = init.state, func = amr, times = times, parms = parms1)
  temp[1,1] <- exploredtau[i]
  temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
  temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
  temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
  tauoutput <- rbind(tauoutput, temp)
}

tauoutput <- data.frame(tauoutput)
colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  

tauoutput$ICombH[tauoutput$tau == 0]/tauoutput$ICombH[tauoutput$tau == averagesales[1]]

# Plotting Prior Distribution ---------------------------------------------

prior.data <- setNames(data.frame(runif(1000, min = 0, max = 0.2), runif(1000, min = 0, max = 0.04), runif(1000, min = 0, max = 2),
                                  "p_alpha" <- rbeta(1000, 1.5, 8.5), runif(1000, min = 0, max = 0.15)), 
                       c("p_betaAA", "p_phi", "p_kappa", "p_alpha", "p_zeta"))

prior.plot.list <- list()

for (i in 1:length(colnames(prior.data))) { # Loop over loop.vector
  
  data <- data.frame("data" = prior.data[,colnames(prior.data)[i]])
  p1 <- ggplot(data, aes(x=data)) + geom_histogram(fill="darkgrey", bins = 10, size = 0.7) + theme_bw() +
    theme(legend.text=element_text(size=14), axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
          plot.title = element_text(size = 14, vjust = 3, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  
  max.y <- max(ggplot_build(p1)$data[[1]][1])*1.2
  
    
  if(colnames(prior.data)[i] == "p_phi") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Resistance Reversion (", phi, ")"))) +
        labs(fill = NULL, title = "") + 
        scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency") 
    }
    if(colnames(prior.data)[i] == "p_kappa") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Recovery (", kappa, ")"))) +
        labs(fill = NULL, title = "") + 
        scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency") 
    }
    if(colnames(prior.data)[i] == "p_betaAA") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")))+
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency")
    }
    if(colnames(prior.data)[i] == "p_alpha") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")"))) +
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y),expand = c(0, 0), name = "Frequency")
    }
    if(colnames(prior.data)[i] == "p_zeta") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) +
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y),expand = c(0, 0), name = "Frequency")
    }
  prior.plot.list[[i]] <- p1
}

comb.prior.plot <- ggarrange(prior.plot.list[[1]],prior.plot.list[[2]],prior.plot.list[[3]],prior.plot.list[[4]],prior.plot.list[[5]],
                             ncol = 1,nrow = 5)

ggsave(comb.prior.plot, filename = "prior_parms.png", dpi = 300, type = "cairo", width = 8, height = 13, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# Supplementary Plots (kappa Responsible for Co-Existence) -----------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)


for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[1,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = MAP[1,"phi"], kappa = 0, alpha = 0, tau = parmtau[i],
             zeta = MAP[1,"zeta"])
  
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  temp[1,8] <- rownames(MAP)[1]
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
output1[,2:5] <- output1[,2:5]*100000 #Scaling the prevalence (per 100,000)

plotdata <- melt(output1,
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

averagesales <- 0.0122887

p1 <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = averagesales, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,7), expand = c(0, 0))  + 
  geom_text(label= c(round(output1$IResRat,digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Tetracycline Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)")  


ggsave(p1, filename = "Icombh_nokappa.png", dpi = 300, type = "cairo", width = 7, height = 4, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")
