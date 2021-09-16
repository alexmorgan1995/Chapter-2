library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("betareg")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")
#setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model
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

#### Data Import ####
#Tetracycline Usage in Fattening Pigs 
datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetra$mgpcuuseage <- datatetra$mgpcuuseage / 1000
datatetra$pig_tetra_sales <- datatetra$pig_tetra_sales / 1000
datatetra <- datatetra[!datatetra$N < 10,]
datatetra$lower <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[1]]))
datatetra$upper <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[2]]))


ggplot()  + geom_point(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim)) +
  geom_text(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

#Ampicillin Usage in Fattening Pigs 

dataamp <- read.csv("resistanceprofAnim_amp.csv")

dataamp$mgpcuuseage <- dataamp$mgpcuuseage / 1000
dataamp$pig_amp_sales <- dataamp$pig_amp_sales / 1000
dataamp <- dataamp[!dataamp$N < 10,]
dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

#Tetracycline Usage in Broiler Poultry 

databroil <- read.csv("salm_broilers_2018.csv")

databroil$mgpcuuseage <- databroil$mgpcuuseage / 1000
databroil$tetra_sales <- databroil$tetra_sales / 1000
databroil <- databroil[!databroil$N < 10,]
databroil$lower <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[1]]))
databroil$upper <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[2]]))

# Beta Regression ------------------------------------------------------
#We are looking at two variables here - Pig_tetra_sales and ResPropAnim 

#Get summary statistics of your variables 
summary(datatetra)
summary(dataamp)
summary(databroil)

#If we wanted to conduct a poisson regression or aany other type of linea regression - at this point we would have to check for
#Colinearity, homoscdasticty and make sure the data conforms to the assumptions for each model 

#Transform the Data to allow for x = (0,1)
datatetra$ResPropAnim_trans <- ((datatetra$ResPropAnim*(nrow(datatetra)-1)) + 0.5)/nrow(datatetra)
dataamp$ResPropAnim_trans <- ((dataamp$ResPropAnim*(nrow(dataamp)-1)) + 0.5)/nrow(dataamp)
databroil$ResPropAnim_trans <- ((databroil$ResPropAnim*(nrow(databroil)-1)) + 0.5)/nrow(databroil)

#Start the beta-regression

#Tet fattening Pigs
output_tet <- betareg(ResPropAnim_trans ~ pig_tetra_sales, data = datatetra)
summary(output_tet)
exp(coef(output_tet))[2]

#Amp fattening Pigs
output_amp <- betareg(ResPropAnim_trans ~ pig_amp_sales, data = dataamp)
summary(output_amp)
exp(coef(output_amp))[2]

#Tet broiler Pigs
output_broil <- betareg(ResPropAnim_trans ~ tetra_sales , data = databroil)
summary(output_broil)
exp(coef(output_broil))[2]


ggplot(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim))  + geom_point() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE) + 
  geom_line(aes(y = predict(output_tet, datatetra), x = pig_tetra_sales), colour = "red", size = 1.2)
  
ggplot(dataamp, aes(x = pig_amp_sales, y= ResPropAnim))  + geom_point() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.03)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE) +
  geom_line(aes(y = predict(output_amp, dataamp), x = pig_amp_sales), colour = "red", size = 1.2)

ggplot(data = databroil, aes(x = tetra_sales, y= ResPropAnim))  + geom_point() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.0225)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")  +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE) + 
  geom_line(aes(y = predict(output_broil, databroil), x = tetra_sales), colour = "red", size = 1.2)
