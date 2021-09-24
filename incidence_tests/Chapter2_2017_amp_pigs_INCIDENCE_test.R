library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data")

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

        
    tCumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    tCumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), tCumS = tCumS, tCumR = tCumR))
  })
}

#### Data Import ####
dataamp <- read.csv("resistanceprofAnim_amp.csv")

dataamp$mgpcuuseage <- dataamp$mgpcuuseage / 1000
dataamp$pig_amp_sales <- dataamp$pig_amp_sales / 1000
dataamp <- dataamp[!dataamp$N < 10,]

dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

mean(dataamp$hum_res_amp ,na.rm = TRUE)
mean(dataamp$pig_amp_sales ,na.rm = TRUE)

ggplot(dataamp, aes(x = pig_amp_sales, y= ResPropAnim))  + geom_point() +
  geom_text(aes(x = pig_amp_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE)



#### Test Data ####
data1 <- cbind(read.csv("INC_RESULTS_1.csv", header = TRUE), "group" = "data1")
data2 <- cbind(read.csv("INC_RESULTS_2.csv", header = TRUE), "group" = "data2")
data3 <- cbind(read.csv("INC_RESULTS_3.csv", header = TRUE), "group" = "data3")
data4 <- cbind(read.csv("INC_RESULTS_4.csv", header = TRUE), "group" = "data4") 
data5 <- cbind(read.csv("INC_RESULTS_5.csv", header = TRUE), "group" = "data5") 
data6 <- cbind(read.csv("INC_RESULTS_6.csv", header = TRUE), "group" = "data6")
data7 <- cbind(read.csv("INC_RESULTS_7.csv", header = TRUE), "group" = "data7")
data8 <- cbind(read.csv("INC_RESULTS_8.csv", header = TRUE), "group" = "data8")
data9 <- cbind(read.csv("INC_RESULTS_9.csv", header = TRUE), "group" = "data9")

map_phi <- map_estimate(data9[,"phi"], precision = 20) 
map_kappa <- map_estimate(data9[,"kappa"], precision = 20) 
map_betaAA <- map_estimate(data9[,"betaAA"], precision = 20) 
map_alpha <- map_estimate(data9[,"alpha"], precision = 20) 
map_zeta <- map_estimate(data9[,"zeta"], precision = 20) 

#### Testing the Model #### 

parmtau <- c(seq(0,0.035,by=0.001), 0.01156391)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00001), phi = map_phi, kappa = map_kappa, alpha = map_alpha, tau = 0.01156391, zeta = map_zeta)
out <- runsteady(y = init, func = amr, parms = parms2, times = c(0, Inf))




for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = map_phi, kappa = map_kappa, alpha = map_alpha, tau = parmtau[i], zeta = map_zeta)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:7] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

output2 <- output1
output2[,2:5] <- output2[,2:5]*100000 #Scaling the prevalence (per 100,000)


plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)", textfont = list(size = 30)) %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,6), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (g/PCU)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                                               xshift =3))
