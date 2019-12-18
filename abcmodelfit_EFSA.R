setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Mathematical Models")
rm(list=ls())
library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2"); library("ggrepel")
library("tidyr"); library("")

#### Model Functions ####
#Foodborne Disease Model with Integrated Lambda and Zeta Parameters
amr <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - alpha*(betaAH*Irh*Sa) - alpha*(betaAA*Ira*Sa) - ua*Sa  
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa
    dIra = alpha*betaAH*Irh*Sa + alpha*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - alpha*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - alpha*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = alpha*(betaHH*Irh*Sh) + alpha*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Model Testbed - Basic Model Output ####
#Initial Parameter Set and Times for the Model
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

times <- seq(0,10000,by=1)

parms = c(ra = 60^-1 , rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.02), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.05, tau = 0.03, theta = 1.2, alpha = 0.9)

out <- ode(y = init, func = amr, times = times, parms = parms) # Solve the ODE

#Transforming the deSolve output for GGplot and plotting

outdata <- data.frame(out); outdata$IComb <- outdata$Ish + outdata$Irh #Outdata is still the original unscaled function 
outdata1 <- outdata; outdata1[5:8] <- outdata[5:8]*100000 # Scale the Human Compartments for (per 100,000 pop)

meltedout <- melt(outdata, id = "time", variable.name = "Compartment", value.name = "Value")
meltedout1 <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#For the Animal Population
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(as.numeric(which(meltedout$Compartment == "Sa" & meltedout$time == 0)):
                                   as.numeric(which(meltedout$Compartment == "Sh" & meltedout$time == 0)) -1),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#For the un-altered human population
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(as.numeric(which((meltedout$Compartment == "Sh" & meltedout$time == 0))):
                                   length(meltedout$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#For the scaled human population - Just for Ih and Irh
ggplot(data = meltedout1, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout1[(as.numeric(which((meltedout1$Compartment == "Ish" & meltedout1$time == 0))):
                                    length(meltedout1$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Basic ICombH/Tau Plot ####
parmtau <- seq(0,0.1,by=0.005)

init <- c(Sa=0.99, Isa=0.01, Ira=0, Sh=1, Ish=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 100000, by = 100)



#parms2 = c(ra = 60^-1 , rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.04), betaAH = 0.00001, betaHH = 0.00001, 
#           betaHA = (0.00001), phi = 0.04, tau = parmtau[i], theta = 0.7, alpha = 0.9)
        
for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 60^-1, rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.05), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.04, theta = 0.4, alpha = 0.9, tau = parmtau[i])
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")

output2 <- output1
output2[,2:5] <- output2[,2:5]*100000 #Scaling the prevalence (per 100,000)

plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)", textfont = list(size = 30)) %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,3.5), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (Tau)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                                               xshift =3))

#plot(parmtau, output1$IResRat)

#### Import Data for Model Fitting ####

datatetra <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter 2 Write Up/resistanceprof.csv")
datatetra$mgpcuuseage <- datatetra$mgpcuuseage / 1000 
cidata <- data.frame()
datatetra <- datatetra[!datatetra$N < 5,]

for (i in 1:length(datatetra$Country)) {
  t <- prop.test(datatetra$Positive.Sample[i], datatetra$N[i])
  cidata[i,1] <- t$conf.int[1]
  cidata[i,2] <- t$conf.int[2]
}

datatetra$lowci <- cidata[,1]; datatetra$highci <- cidata[,2]
init <- c(Sa=0.99, Isa=0.01, Ira=0, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200000, by = 10)
usagedata <- rbind(data.frame("Parameter" = "Observed", "Usage" = c(datatetra$mgpcuuseage)),
                   data.frame("Parameter" = "Model", "Usage" = seq(0,0.1,by=0.01)))

parms = c(ra = 60^-1 , rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.05), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.05, tau = parmtau[i], theta = 1.2, alpha = 0.9)

for (j in 1:length(unique(usagedata$Parameter))) {
  tempusage <- data.frame()
  parmtau <- usagedata[,2][usagedata == as.character(unique(usagedata[,1])[j])]
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
    parms2 = c(ra = 60^-1 , rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.05), betaAH = 0.00001, betaHH = 0.00001, 
              betaHA = (0.00001), phi = 0.046, tau = parmtau[i], theta = 0.8, alpha = 0.9)
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[i]
    temp[1,2] <- rounding(out[nrow(out),3]) 
    temp[1,3] <- rounding(out[nrow(out),4]) 
    temp[1,4] <- temp[1,2] + temp[1,3]
    temp[1,5] <- temp[1,3]/temp[1,4]
    print(temp[1,5])
    tempusage <- rbind.data.frame(tempusage, temp)
  }
  colnames(tempusage)[1:5] <- c("tau", "SensAnim","ResAnim","ICombA","ResPropAnim")
  assign(paste("usage", as.character(unique(usagedata[,1])[j]), sep=""), tempusage)
}

usageObserved$country <- c(as.character(datatetra$Country)); usageObserved$isolates <- c(datatetra$N)

#Plotting Model Predictions - Resistance
ggplot(data = usageObserved[!usageObserved$isolates < 5,], aes(x = tau, y = ResPropAnim)) + geom_point() + 
  geom_text_repel(data = usageObserved[!usageObserved$isolates < 5,], 
                  mapping=aes(x=usageObserved$tau[!usageObserved$isolates < 5], 
                              y= usageObserved$ResPropAnim[!usageObserved$isolates < 5], label=usageObserved$country[!usageObserved$isolates < 5]), 
                  size=4, box.padding = unit(1.5, "lines")) +
  geom_line(data = usageModel, aes(x = usageModel$tau, usageModel$ResPropAnim, col = "red")) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

ggplot(data = datatetra[!datatetra$N < 5,], aes(x = mgpcuuseage, y = ResPropAnim)) + geom_point() +
  geom_text_repel(data = datatetra[!datatetra$N < 5,], 
                  mapping=aes(x=datatetra$mgpcuuseage[!datatetra$N < 5], 
                              y= datatetra$ResPropAnim[!datatetra$N < 5], label=datatetra$Country[!datatetra$N < 5]), 
                  size=4, box.padding = unit(1.5, "lines")) +
  geom_line(data = usageModel, aes(x = usageModel$tau, usageModel$ResPropAnim, col = "red")) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) 

#Plotting CIs on graph

ggplot(data = datatetra[!datatetra$N < 5,], aes(x = mgpcuuseage, y = ResPropAnim)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = datatetra$lowci[!datatetra$N < 5], ymax = datatetra$highci[!datatetra$N < 5])) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0))

#### ABC Model Fitting #### 
summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

sum_square_diff_dist(list(summarystatprev), datatetra, usageObserved)

#This I guess is very similar to the for loop I programmed before - which will output a distance for any given parameter set or initial conditions. 

times <-  seq(0, 5000, by = 100)
init <- c(Sa=0.99, Isa=0.01, Ira=0, Sh=1, Ish=0, Irh=0)
thetaparm = c(ra = 60^-1 , rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.05), betaAH = 0.00001, betaHH = 0.00001, 
                     betaHA = (0.00001), phi = 0.04, theta = 0.4, alpha = 0.9)

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol=3)
  tau_range <- append(tau_range, 0.032)
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol=3)
    parms2 = c(ra = thetaparm[["ra"]], rh =  thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], 
               betaAA = thetaparm[["betaAA"]], betaAH = thetaparm[["betaAH"]], betaHH = thetaparm[["betaHH"]], 
               betaHA = thetaparm[["betaHA"]], phi = thetaparm[["phi"]], tau = tau_range[i], theta = thetaparm[["theta"]], 
               alpha = thetaparm[["alpha"]])
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim")  
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == 0.032,]), 
         abs(tauoutput$ICombH[tauoutput$tau == 0.032] - 3.26),
         abs(tauoutput$ResPropAnim[tauoutput$tau == 0.032] - 0.409)))
}

test_dist <- computeDistanceABC_ALEX(sum.stats = summarystatprev, 
                                     distanceABC = sum_square_diff_dist,  
                                     fitmodel = amr, 
                                     tau_range = datatetra$mgpcuuseage, 
                                     thetaparm = thetaparm, 
                                     init.state = init, 
                                     times = times,
                                     data = datatetra)

# use the ABC rejection algorithm to find population 1 in the ABC-SMC algorithm
start_time <- Sys.time()

times <-  seq(0, 10000, by = 100)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

hist(runif(1000, min = 0, max = 0.1))
hist(runif(1000, min = 0, max = 2))
hist(runif(1000, min = 0.5, max = 1))
hist(runif(1000, min = 0, max = 0.1))

ABC_algorithm <- function(N, epsilon, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  dump1 <- matrix(nrow = 0 , ncol = 11)
  i <- 0
  while(i < N) {
    d_phi <- runif(1, min = 0.02, max = 0.09)
    d_theta <- runif(1, min = 0, max = 2)
    d_alpha <- runif(1, min = 0.7, max = 1)
    d_betaAA <- runif(1, min = 0, max = 0.15)
    thetaparm <- c(ra = 60^-1, rh =  (7^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
      betaHA = (0.00001), phi = d_phi, theta = d_theta, alpha = d_alpha)
    dist <- computeDistanceABC_ALEX(sum.stats = summarystatprev, 
                                    distanceABC = sum_square_diff_dist, 
                                    fitmodel = amr, 
                                    tau_range = datatetra$mgpcuuseage, 
                                    thetaparm = thetaparm, 
                                    init.state = init, 
                                    times = times,
                                    data = datatetra)
    if((dist[1] <= epsilon) && (!is.na(dist)) && (dist[2] <= epsilon[2]) && (dist[3] <= epsilon[3])) {
      dump1 <- rbind(dump1, thetaparm)
    }
    i <- dim(dump1)[1]
    print(c(i, dist, thetaparm[c(5,9:11)]))
  }
  return(dump1)
}

data_ABC <- ABC_algorithm(N = 1000, 
                          epsilon = c(0.85, 0.489, 0.0614),
                          sum.stats = list(summarystatprev), 
                          distanceABC = sum_square_diff_dist, 
                          fitmodel = amr, 
                          tau_range = datatetra$mgpcuuseage, 
                          init.state = init, 
                          times = times, 
                          data = datatetra)

end_time <- Sys.time()
end_time - start_time

saved <- data_ABC

hist(data_ABC[,5], xlim = c(0,0.1))
hist(data_ABC[,9], xlim = c(0,0.1))
hist(data_ABC[,10], xlim = c(0,2))
hist(data_ABC[,11], xlim = c(0.5,1))
