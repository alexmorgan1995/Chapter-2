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
            "theta" <- mean(data$theta[which(data$group == "data10" & data$fit == "tet")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "tet")]),
            "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "tet")]),
            "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "tet")]))

MAPamp <- c("phi" <- mean(data$phi[which(data$group == "data10" & data$fit == "amp")]),
            "theta" <- mean(data$theta[which(data$group == "data10" & data$fit == "amp")]),
            "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "amp")]),
            "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "amp")]),
            "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "amp")]))

MAPbroil <- c("phi" <- mean(data$phi[which(data$group == "data10" & data$fit == "broil")]),
              "theta" <- mean(data$theta[which(data$group == "data10" & data$fit == "broil")]),
              "betaAA" <- mean(data$betaAA[which(data$group == "data10" & data$fit == "broil")]),
              "alpha" <- mean(data$alpha[which(data$group == "data10" & data$fit == "broil")]),
              "zeta" <- mean(data$zeta[which(data$group == "data10" & data$fit == "broil")]))

MAP <- rbind(MAPtet, MAPamp, MAPbroil); colnames(MAP) <- c("phi", "theta", "betaAA", "alpha", "zeta")

# Model Fit with Data -----------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:nrow(MAP)) {
  output1 <- data.frame()
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
    
    if(rownames(MAP)[j] == "MAPbroil") {
      parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])
    } 
    
    else {
      parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i],
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
    print(temp[1,3])
    output1 <- rbind.data.frame(output1, temp)
  }
  icombhdata <- rbind(icombhdata, output1)
}

colnames(icombhdata)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
icombhdata[,2:5] <- icombhdata[,2:5]*100000 #Scaling the prevalence (per 100,000)

amp <- ggplot(dataamp, aes(x = pig_amp_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.03),expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Ampicillin Sales (g/PCU)", y = "Ampicillin-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPtet",], aes(x = tau, y= IResRatA), col = "red", size = 1.02) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

tet <- ggplot(datatetra, aes(x = pig_tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.034), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPamp",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

broil <- ggplot(databroil, aes(x = tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.022), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Broiler Poultry Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Broiler Poultry Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

fits <- ggarrange(tet, amp, broil, nrow = 1, ncol = 3, align = "h", labels = c("A","B", "C"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(fits, filename = "Model_Fits.png", dpi = 300, type = "cairo", width = 14, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")


# Ribbon ------------------------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 3, nrow = 0))
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
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])
    } 
    else {
      parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])
    }
    
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[i]
    temp[1,2] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[1,3] <- rownames(MAP)[j]
    
    for(z in 1:1000) {
      temp_ribbon_tau <- data.frame(matrix(NA, nrow = 1, ncol=4))
      
      if(rownames(MAP)[j] == "MAPbroil") {
        parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], theta = ribbondata[z,"theta"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])
      } 
      else {
        parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi = ribbondata[z,"phi"], theta = ribbondata[z,"theta"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
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

colnames(icombhdata)[1:3] <- c("tau","IResRatA", "group")
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