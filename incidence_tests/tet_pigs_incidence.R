library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")

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
    
    dCumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    dCumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh,dCumS,dCumR)))
  })
}

#### Data Import ####
datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetrahum <- read.csv("resistanceprofHum_v1.csv")
datatetra$mgpcuuseage <- datatetra$mgpcuuseage / 1000; datatetrahum$mgpcuuseage <- datatetrahum$mgpcuuseage / 1000
datatetra$pig_tetra_sales <- datatetra$pig_tetra_sales / 1000; datatetrahum$pig_tetra_sales <- datatetrahum$pig_tetra_sales / 1000
datatetrahum$ResPropHum <- datatetrahum$ResPropHum/ 100 
datatetra <- datatetra[!datatetra$N < 10,]

mean(datatetra$ResPropHum ,na.rm = TRUE)
mean(datatetra$pig_tetra_sales ,na.rm = TRUE)

ggplot()  + geom_point(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim)) +
  geom_text(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

ggplot()  + geom_point(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum)) +
  geom_text(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Human Carriage")

#### Test Data ####

data10 <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit/results_ABC_SMC_gen_tet_10.csv", header = TRUE)

map_tet <- map_estimate(data10)

#### Testing the Model #### 

parmtau <- seq(0, 0.03, by=0.001)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0, CumS = 0, CumR = 0)
output1 <- data.frame()
times <- seq(0, 200000, by = 50)

parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_tet[map_tet$Parameter == "betaAA", 2], betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00001), phi = map_tet[map_tet$Parameter == "phi", 2], kappa = map_tet[map_tet$Parameter == "kappa", 2], 
           alpha = map_tet[map_tet$Parameter == "alpha", 2], tau = 0.001, zeta = map_tet[map_tet$Parameter == "zeta", 2])
out <- ode(y = init, func = amr, times = times, parms = parms2)


for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol= 11))
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_tet[map_tet$Parameter == "betaAA", 2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = map_tet[map_tet$Parameter == "phi", 2], kappa = map_tet[map_tet$Parameter == "kappa", 2], 
             alpha = map_tet[map_tet$Parameter == "alpha", 2], tau = parmtau[i], zeta = map_tet[map_tet$Parameter == "zeta", 2])
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5])
  temp[1,3] <- rounding(out[nrow(out),6])*100000
  temp[1,4] <- rounding(out[nrow(out),7])*100000
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  temp[1,8] <- tail(diff(as.matrix(out[,"CumS"])), 1)
  temp[1,9] <- tail(diff(as.matrix(out[,"CumR"])), 1)
  temp[1,10] <- temp[1,8] + temp[1,9] 
  temp[1,11] <- round((temp[1,9] / temp[1,10]), digits = 3) 
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:11] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA","IncS","IncR", "comb_inc", "comb_res")

output2 <- output1
output2$IncS <- (output1$IncS*446000000)/100000
output2$IncR <- (output1$IncR*446000000)/100000
output2$comb_inc <- output2$IncS + output2$IncR

plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)", textfont = list(size = 30)) %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,6), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (g/PCU)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                                               xshift =3))

plot_ly(output2, x= ~tau, y = ~IncS, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~IncR, name = "Antibiotic-Resistant Infection (Human)", textfont = list(size = 30)) %>% 
  layout(yaxis = list(title = "Incidence of EU Human Salmonellosis (ICombH) per 100,000", range = c(0,2.5), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (g/PCU)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T, barmode = "stack", 
         annotations = list(x = ~tau, y = ~comb_inc, text = ~comb_res, yanchor = "bottom", showarrow = FALSE, textangle = 310, xshift =3))
