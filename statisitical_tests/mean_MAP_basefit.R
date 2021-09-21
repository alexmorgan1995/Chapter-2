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

# Testing the Model -------------------------------------------------------

thetaparm_post_tetra <- read.csv(paste0("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit/",
                             list.files("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit", pattern = "tet_10")))
thetaparm_post_amp <- read.csv(paste0("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit/",
                                  list.files("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit", pattern = "amp_9")))
thetaparm_post_broil <- read.csv(paste0("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit/",
                                  list.files("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit/Remote_Fit", pattern = "broil_10")))
#Functions
summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

#Parameters

init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times = seq(0, 2000, by = 50)

tau_list <- list("tetra" = datatetra$pig_tetra_sales ,
                 "amp" = dataamp$pig_amp_sales ,
                 "broil" = databroil$tetra_sales)

map_list <- list("tetra" = map_estimate(thetaparm_post_tetra),
                 "amp" = map_estimate(thetaparm_post_amp),
                 "broil" = map_estimate(thetaparm_post_broil))

#map_list <- list("tetra" = data.frame("Parameter" = colnames(thetaparm_post_tetra), "mean" = colMeans(thetaparm_post_tetra)),
#                 "amp" = data.frame("Parameter" = colnames(thetaparm_post_amp), "mean" = colMeans(thetaparm_post_amp)),
#                 "broil" = data.frame("Parameter" = colnames(thetaparm_post_broil), "mean" = colMeans(thetaparm_post_broil)))


#Running the Model

test_stat <- list()

for(j in 1:3) {
  test_stat[[j]] <- local({
    tauoutput = matrix(nrow = 0, ncol=4)
    data <- list(datatetra, dataamp, databroil)[[j]]
    tau_range = data[[ncol(data)-4]]
    thetaparm = map_list[[j]]
    
    for (i in 1:length(tau_range)) {
      temp <- matrix(NA, nrow = 1, ncol=4)
      parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua =  240^-1, uh = 28835^-1, 
                 betaAA = thetaparm[thetaparm$Parameter == "betaAA",2], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = 0.00001, phi = thetaparm[thetaparm$Parameter == "phi",2], tau = tau_range[i], kappa = thetaparm[thetaparm$Parameter == "kappa",2], 
                 alpha = thetaparm[thetaparm$Parameter == "alpha",2], zeta = thetaparm[thetaparm$Parameter == "zeta",2])
      out <- ode(y = init.state, func = amr, times = times, parms = parms2)
      temp[1,1] <- tau_range[i]
      temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
      temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
      temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
      tauoutput <- rbind(tauoutput, temp)
    }
    tauoutput <- data.frame(tauoutput)
    colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")
    
    tauoutput_list <- list("dyn" = data.frame("country" = data$Country, "tau" = tau_range, "ResPropAnim" = tauoutput$ResPropAnim),
                           "stat" = data.frame("country" = data$Country, "tau" = tau_range, "ResPropAnim" = predict(list(output_tet, output_amp, output_broil)[[j]], data)))
    tauoutput_list[["dyn_sumsquare"]] <- sum_square_diff_dist(list(summarystatprev), data, tauoutput_list[[1]])
    tauoutput_list[["stat_sumsquare"]] <- sum_square_diff_dist(list(summarystatprev), data, tauoutput_list[[2]])
  
    p1 <- ggplot(tauoutput_list[["dyn"]], aes(x = tau, y = ResPropAnim)) + geom_line(col = "red") + 
      labs(title = paste0(c("Tet Fattening Pigs", "Amp Fattening Pigs", "Tet Broiler Poultry")[j],": Dynamic Model Prediction"),
    x = "Livestock Antibiotic Usage") +  
      geom_point(aes(x = data[[ncol(data)-4]], y = data[,3]), inherit.aes = FALSE) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + 
      annotate("text", x =  max(tau_range)/2, y = 0.9, label = paste0("Sum of Squares: " ,round(tauoutput_list[["dyn_sumsquare"]], 4)), size = 5)
    
    p2 <- ggplot(tauoutput_list[["stat"]], aes(x = tau, y = ResPropAnim)) + geom_line(col = "red") + 
      labs(title = paste0(c("Tet Fattening Pigs", "Amp Fattening Pigs", "Tet Broiler Poultry")[j],": Statistical Model Prediction"),
           x = "Livestock Antibiotic Usage") + 
      geom_point(aes(x = data[[ncol(data)-4]], y = data[,3]), inherit.aes = FALSE) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + 
      annotate("text", x = max(tau_range)/2, y = 0.9, label = paste0("Sum of Squares: " ,round(tauoutput_list[["stat_sumsquare"]], 4)), size = 5)
    
    tauoutput_list[["comb_plot"]] <- ggarrange(p1,p2, ncol = 1, nrow = 2)
    
    return(tauoutput_list)
  })
}

test_stat

