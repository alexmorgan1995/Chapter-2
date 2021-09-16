library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("betareg")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Data")
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
datafluor<- read.csv("test_fluoroquine.csv")
datafluor$scaled_sales_fluor <- datafluor$scaled_sales_fluor / 1000
datafluor <- datafluor[!datafluor$num_test_fluor  < 10,]

datafluor$lower <- unlist(lapply(1:nrow(datafluor), function(i) prop.test(datafluor$isol_res_fluor[i],datafluor$num_test_fluor[i])[[6]][[1]]))
datafluor$upper <- unlist(lapply(1:nrow(datafluor), function(i) prop.test(datafluor$isol_res_fluor[i],datafluor$num_test_fluor[i])[[6]][[2]]))

ggplot()  + geom_point(data = datafluor, aes(x = scaled_sales_fluor, y= propres_fluor )) +
  geom_text(data = datafluor, aes(x = scaled_sales_fluor, y= propres_fluor , label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.005)) + scale_y_continuous(expand = c(0, 0), limits = c(0,0.35)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

#Ampicillin Usage in Fattening Pigs 

datasulfa<- read.csv("test_sulfa.csv")
datasulfa$scaled_sales_sulfa <- datasulfa$scaled_sales_sulfa / 1000
datasulfa <- datasulfa[!datasulfa$num_test_sulfa  < 10,]

datasulfa$lower <- unlist(lapply(1:nrow(datasulfa), function(i) prop.test(datasulfa$isol_res_sulfa[i],datasulfa$num_test_sulfa[i])[[6]][[1]]))
datasulfa$upper <- unlist(lapply(1:nrow(datasulfa), function(i) prop.test(datasulfa$isol_res_sulfa[i],datasulfa$num_test_sulfa[i])[[6]][[2]]))

ggplot()  + geom_point(data = datasulfa, aes(x = scaled_sales_sulfa, y= propres_sulfa )) +
  geom_text(data = datasulfa, aes(x = scaled_sales_sulfa, y= propres_sulfa , label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.01)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

# Beta Regression ------------------------------------------------------
#We are looking at two variables here - Pig_tetra_sales and ResPropAnim 

#Get summary statistics of your variables 
summary(datafluor)
summary(datasulfa)

#If we wanted to conduct a poisson regression or aany other type of linea regression - at this point we would have to check for
#Colinearity, homoscdasticty and make sure the data conforms to the assumptions for each model 

#Transform the Data to allow for x = (0,1)
datafluor$ResPropAnim_trans <- ((datafluor$propres_fluor*(nrow(datafluor)-1)) + 0.5)/nrow(datafluor)
datasulfa$ResPropAnim_trans <- ((datasulfa$propres_sulfa*(nrow(datasulfa)-1)) + 0.5)/nrow(datasulfa)

#Start the beta-regression

#Tet fattening Pigs
output_fluor <- betareg(ResPropAnim_trans ~ scaled_sales_fluor, data = datafluor)
summary(output_fluor)
exp(coef(output_fluor))[2]

#Amp fattening Pigs
output_sulfa <- betareg(ResPropAnim_trans ~ scaled_sales_sulfa, data = datasulfa)
summary(output_sulfa)
exp(coef(output_sulfa))[2]

p1 <- ggplot(data = datafluor, aes(x = scaled_sales_fluor, y= propres_fluor))  + geom_point() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.0045)) + scale_y_continuous(expand = c(0, 0), limits = c(0,0.5)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE) + 
  geom_line(aes(y = predict(output_fluor, datafluor), x = scaled_sales_fluor), colour = "red", size = 1.2)+
  labs(title = "Fluoroquinolone Usage in Fattening Pigs")
  
p2 <- ggplot(datasulfa, aes(x = scaled_sales_sulfa, y= propres_sulfa))  + geom_point() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.009)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE) +
  geom_line(aes(y = predict(output_sulfa, datasulfa), x = scaled_sales_sulfa), colour = "red", size = 1.2)+
  labs(title = "Sulfamethoxazole Usage in Fattening Pigs")

ggarrange(p1,p2, nrow = 2, ncol = 1)