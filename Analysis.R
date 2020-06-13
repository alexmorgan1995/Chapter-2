library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

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

# Data Import -------------------------------------------------------------
#Import of Posterior Distributions
data <- do.call(rbind, list(import(tet), import(amp), import(broil)))

#Import of Fitting Data
dataamp <- read.csv("resistanceprofAnim_amp.csv")
dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetra$lower <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[1]]))
datatetra$upper <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[2]]))

databroil <- read.csv("salm_broilers_2018.csv")
databroil$lower <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[1]]))
databroil$upper <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[2]]))

#Obtain the MAPs for each dataset

MAPtet <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "tet")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "tet")]),
            "betaAA" <- map_estimate(data$d_betaAA[which(data$group == "data5" & data$fit == "tet")]),
            "alpha" <- map_estimate(data$d_alpha[which(data$group == "data5" & data$fit == "tet")]))

MAPamp <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "amp")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "amp")]),
            "betaAA" <- map_estimate(data$d_betaAA[which(data$group == "data5" & data$fit == "amp")]),
            "alpha" <- map_estimate(data$d_alpha[which(data$group == "data5" & data$fit == "amp")]))

MAPbroil <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "broil")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "broil")]),
            "betaAA" <- map_estimate(data$d_betaAA[which(data$group == "data5" & data$fit == "broil")]),
            "alpha" <- map_estimate(data$d_alpha[which(data$group == "data5" & data$fit == "broil")]))

MAP <- rbind(MAPtet, MAPamp, MAPbroil); colnames(MAP) <- c("phi", "theta", "betaAA", "alpha")

# ABC-SMC Posterior -------------------------------------------------------

testphi <- melt(rbind(data1, data2, data3, data4, data5), id.vars = "group",measure.vars = "phi")
testtheta <- melt(rbind(data1, data2, data3, data4,data5), id.vars = "group",measure.vars = "theta")
testbetaAA <- melt(rbind(data1, data2, data3, data4, data5), id.vars = "group",measure.vars = "d_betaAA")
testalpha <- melt(rbind(data1, data2, data3, data4, data5), id.vars = "group",measure.vars = "d_alpha")

p1 <- ggplot(testphi, aes(x=value, fill=group)) + geom_density(alpha=.5) + 
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")"))) + 
  scale_y_continuous(limits = c(0,350), expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p2 <- ggplot(testtheta, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", theta, ")"))) + 
  scale_y_continuous(limits = c(0,45), expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p3<- ggplot(testbetaAA, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")"))) + 
  scale_y_continuous(limits = c(0,40),expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p4 <- ggplot(testalpha, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")"))) + 
  scale_y_continuous(limits = c(0,12),expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# Model Fit with Data -----------------------------------------------------

parmtau <- seq(0,0.035, by = 0.001)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:nrow(MAP)) {
  output1 <- data.frame()
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
    parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
               betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i])
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
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPtet",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)

tet <- ggplot(datatetra, aes(x = pig_tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPamp",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)

broil <- ggplot(databroil, aes(x = tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)

ggarrange(amp, tet, broil, nrow = 1, ncol = 3) 

# Tau and ICombH ----------------------------------------------------------


