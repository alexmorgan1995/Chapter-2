setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Mathematical Models")

rm(list=ls()); library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("tidyr")

#### Model Functions ####
#Foodborne Disease Model with Integrated Lambda and Zeta Parameters
amr <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSa = ua + ra*(Ia + Ira) + zeta*tau*Ia - (betaAA*Ia*Sa) - (betaAH*Ih*Sa) - lambda*(betaAH*Irh*Sa) - lambda*(betaAA*Ira*Sa) - ua*Sa  
    dIa = betaAA*Ia*Sa + betaAH*Ih*Sa + phi*Ira - zeta*tau*Ia - tau*theta*Ia - ra*Ia - ua*Ia
    dIra = lambda*betaAH*Irh*Sa + lambda*betaAA*Ira*Sa + tau*theta*Ia - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ih+Irh) - (betaHH*Ih*Sh) - lambda*(betaHH*Irh*Sh) - (betaHA*Ia*Sh) - lambda*(betaHA*Ira*Sh) - uh*Sh 
    dIh = betaHH*Ih*Sh + betaHA*Ia*Sh - rh*Ih - uh*Ih 
    dIrh = lambda*(betaHH*Irh*Sh) + lambda*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIa,dIra,dSh,dIh,dIrh)))
  })
}

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Model Testbed - Basic Model Output ####
#Initial Parameter Set and Times for the Model
init <- c(Sa=0.98, Ia=0.01, Ira=0.01, Sh=1, Ih=0, Irh=0)

times <- seq(0,10000,by=1)

parms = c(ra = 0, rh =  7^-1, uh = 28835^-1, ua = 42^-1, betaAA = 0.0415, betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = 0.00001, phi = 0.02, tau = 0.02, theta = 0.5, zeta = 1, lambda = 1)

out <- ode(y = init, func = amr, times = times, parms = parms) # Solve the ODE

#Transforming the deSolve output for GGplot and plotting

outdata <- data.frame(out)
outdata$IComb <- outdata$Ih + outdata$Irh #Outdata is still the original unscaled function 
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
  geom_line(data=(y = meltedout[(as.numeric(which((meltedout$Compartment == "Ih" & meltedout$time == 0))):
                                   length(meltedout$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#For the scaled human population - Just for Ih and Irh
ggplot(data = meltedout1, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(as.numeric(which((meltedout1$Compartment == "Ih" & meltedout1$time == 0))):
                                   length(meltedout1$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Sensitivity Analysis ####

start_time <- Sys.time()

parms = fast_parameters(minimum = c(0, 70^-1, 420^-1, 288350^-1, 0.00415, 0.000001, 0.000001, 0.000001, 0, 0, 0), 
                        maximum = c(420^-1, 0.7^-1, 4.2^-1, 2883.5^-1, 0.415, 0.0001, 0.0001, 0.0001, 0.2, 0.2, 2), 
                             factor=11, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                  "phi", "tau", "theta"))

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0,200000, by = 10) 
output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms1 = c(ra = parms$ra[i], rh = parms$rh[i], ua = parms$ua[i], uh = parms$uh[i], betaAA = parms$betaAA[i],
             betaAH = parms$betaAH[i], betaHH = parms$betaHH[i], betaHA = parms$betaHA[i], phi=parms$phi[i],
             tau=parms$tau[i], theta=parms$theta[i], zeta = 1, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  temp[1,1] <- rounding(out[nrow(out),5]) 
  temp[1,2] <- rounding(out[nrow(out),6]) 
  temp[1,3] <- rounding(out[nrow(out),7])
  temp[1,4] <- temp[1,2] + temp[1,3]
  temp[1,5] <- temp[1,3]/temp[1,4]
  temp[1,6] <- signif(((rounding(out[nrow(out),6]) + rounding(out[nrow(out)-1,6]) + rounding(out[nrow(out)-2,6]))/3), digits = 6)
  if(temp[1,6] == temp[1,2]) {temp[1,7] <- "YES"}
  print(temp[1,7])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:7] <- c("SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat", "ValueatEqui","Equi?")
end_time <- Sys.time(); end_time - start_time

#### Sensitivity Analysis ####
#For ICombH
sensit <- output$ICombH #Creating Variable for the output variable of interest
sens <- sensitivity(x=sensit, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                               "phi", "tau", "theta"))
df.equilibrium <- NULL; df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi","tau", "theta"), value=sens)

p <- ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23"); p

#For ResRat (Humans)

sensit1 <- output$IResRat #Creating Variable for the output variable of interest
sensit1[is.nan(sensit1)] <- 0
sens1 <-sensitivity(x=sensit1, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                  "phi", "tau", "theta"))
df.equilibrium1 <- NULL; df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi","tau", "theta"), value=sens1)

p1 <- ggplot(df.equilibrium11, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23"); p1

#Combined Bar Plot with Both Sensitivity Analysis
df.equilibrium$state <- "Overall Foodborne Disease"
df.equilibrium1$state <- "Resistance Ratio"

newdf <- rbind(df.equilibrium, df.equilibrium1)

ggplot(data = newdf, aes(x = parameter, y = value, fill = state)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_x_discrete("Parameter", waiver(), c(expression(paste(beta[AA])), expression(paste(beta[AH])),
                                            expression(paste(beta[HA])), expression(paste(beta[HH])),
                                            expression(paste(phi)), expression(paste(italic("r")[A])),
                                            expression(paste(italic("r")[H])), expression(paste(tau)),
                                            expression(paste(theta)), expression(paste(mu[A])),
                                            expression(paste(mu[H])))) +
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0.00009)) + 
  scale_fill_manual(values=c("dodgerblue2", "orangered1")) +
  theme(axis.text=element_text(size = 12, colour = "black"),
        axis.line.x = element_line(color="black", size = 0.7),
        axis.title=element_text(size=12), 
        legend.position=c(0.8, 0.9), legend.title = element_blank()) +
  ylab("Partial Variance") 

#### Basic ICombH/Tau Plot ####
parmtau <- seq(0,0.1,by=0.005)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 0 , rh =  (7^-1), ua = 42^-1, uh = 28835^-1, betaAA = (0.0415), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.02, tau = parmtau[i], theta = 0.5, zeta = 1, lambda = 0.7)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
icombrat <- output1$ICombH[output1$tau == 0.5]/output1$ICombH[output1$tau == 0]
output1$IResRat <- signif(output1$IResRat, digits = 3)

output2 <- output1
output2$InfHumans <- output2$InfHumans*100000; 
output2$ResInfHumans <- output2$ResInfHumans*100000
output2$ICombH <- output2$ICombH*100000

table <- data.frame(x = c("Baseline Curtailment", "70% BetaAA", "70% BetaHA", "70% BetaAA/BetaHA"),
                    y = c(3.055010, 1.7933500, 2.138530, 1.2553500))
table$x <- factor(table$x, levels = c("Baseline Curtailment", "70% BetaAA", "70% BetaHA", "70% BetaAA/BetaHA"))

plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)") %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,3.5), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (Tau)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", 
         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                            xshift =3))

#### Comb IRH and IH Measure Testing - Effect of Treatment - NEW PLOT ####

#Testing for the Effect of Changing Treatment on the Combined Measure

#parmtau <- seq(0,0.5,by=0.03)
parmtau <- seq(0,0.1,by=0.0005)
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 200000, by = 100)

#Have to alter a number of the outputnames to generate all the values for the figures
outputua13 <- data.frame()

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 0, rh =  (7^-1), ua = (42^-1)*1.3, uh = 28835^-1, betaAA = (0.0415), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = (0.02), tau = parmtau[i], theta = (0.5))
  out <- ode(y = init, func = amrold, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,3])
  outputua13 <- rbind.data.frame(outputua13, temp)
}

colnames(outputua13)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")

outputua13$InfHumans <- outputua13$InfHumans*100000
outputua13$ResInfHumans <- outputua13$ResInfHumans*100000
outputua13$ICombH <- outputua13$ICombH*100000

#### Plotting ####
outputAA1$state <- "Baseline" #BetaHA = 0.00001; BetaAA = 0.1
outputAA9$state <- "BetaHA = 100%; BetaAA = 90%" #BetaHA = 0.00001; BetaAA = 0.09 
outputAA8$state <- "BetaHA = 100%; BetaAA = 80%" #BetaHA = 0.00001; BetaAA = 0.08
outputAA7$state <- "BetaHA = 100%; BetaAA = 70%" #BetaHA = 0.00001; BetaAA = 0.07
new.data3 <- do.call("rbind", list(outputAA1, outputAA9, outputAA8, outputAA7))

ggplot(data = new.data3, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.00005), expand = c(0,0)) +
  scale_colour_manual(values = c("black","steelblue1","slateblue1","slategray"), 
                      labels = c("Baseline   ", 
                                 expression(paste(beta[AA]," = 70%")), 
                                 expression(paste(beta[AA]," = 80%")), 
                                 expression(paste(beta[AA]," = 90%"))))+
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +

ggplot() + 
  geom_line(data=outputAA1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputAA1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputAA9, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputAA9, aes(tau, replace(ICombH, ICombH>1.26, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputAA8, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputAA8, aes(tau, replace(ICombH, ICombH>1.26, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputAA7, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputAA7, aes(tau, replace(ICombH, ICombH>1.28, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste(" Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                " per 100,000"))) +
  theme(text = element_text(size=13.5),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputAA1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputAA9, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputAA8, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputAA7, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##

outputHA1$state <- "Baseline Parameters" #BetaHA = 0.00001; BetaAA = 0.1
outputHA9$state <- "Beta" #BetaHA = 0.000009; BetaAA = 0.1
outputHA8$state <- "Eight" #BetaHA = 0.000008; BetaAA = 0.1
outputHA7$state <- "Seven" #BetaHA = 0.000007; BetaAA = 0.1
new.data4 <- do.call("rbind", list(outputHA1, outputHA9, outputHA8, outputHA7))

ggplot(data = new.data1, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.00005), expand = c(0,0)) +
  scale_colour_manual(values = c("black","steelblue1","slateblue1","slategray"), 
                      labels = c("Baseline   ", 
                                 expression(paste(beta[HA]," = 70%")), 
                                 expression(paste(beta[HA]," = 80%")), 
                                 expression(paste(beta[HA]," = 90%"))))+
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_hline(yintercept = 1.595806e-05, linetype = 2, size = 1) + 
  annotate("text", 0.42, 1.595806e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputHA1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputHA1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputHA9, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputHA9, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputHA8, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputHA8, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputHA7, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputHA7, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                          " per 100,000")))  +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputHA1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputHA9, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputHA8, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputHA7, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##

outputphi12$state <- "Twelve" #BetaHA = 0.000009; BetaAA = 0.1
outputphi11$state <- "Eleven" #BetaHA = 0.000009; BetaAA = 0.0
outputphi1$state <- "One" #BetaHA = 0.000009; BetaAA = 0.08
outputphi09$state <- "Nine" #BetaHA = 0.000009; BetaAA = 0.07
outputphi08$state <- "Eight" #BetaHA = 0.000009; BetaAA = 0.07
new.data5 <- do.call("rbind", list(outputphi12, outputphi11, outputphi1, outputphi09, outputphi08))

ggplot(data = new.data5, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                          " per 100,000"))) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  scale_colour_manual(values = c("steelblue1","slateblue","black","firebrick","red"), 
                      labels = c(expression(paste(phi," = 80%")), 
                                 expression(paste(phi," = 90%")),
                                 "Baseline", 
                                 expression(paste(phi," = 110%")), 
                                 expression(paste(phi," = 120%"))))+
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.165, 1.388e-05, vjust = -0.6, label = "Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputphi12, aes(tau,ICombH), linetype= 3, size = 1, colour = "red") + 
  geom_line(data=outputphi12, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "red") +
  geom_line(data=outputphi11, aes(tau,ICombH), linetype= 3, size = 1, colour = "firebrick") + 
  geom_line(data=outputphi11, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "firebrick") +
  geom_line(data=outputphi1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputphi1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputphi09, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputphi09, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputphi08, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputphi08, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                " per 100,000"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputphi12, aes(tau,IResRat), linetype= 1, size = 1, colour = "red") +
  geom_line(data=outputphi11, aes(tau,IResRat), linetype= 1, size = 1, colour = "firebrick") +
  geom_line(data=outputphi1, aes(tau,IResRat), linetype= 1, size = 1, colour = "black") + 
  geom_line(data=outputphi09, aes(tau,IResRat), linetype= 1, size = 1, colour = "slateblue1") +
  geom_line(data=outputphi08, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##

outputrh1$state <- "Baseline" #BetaHA = 0.000008; BetaAA = 0.1
outputrh11$state <- "Eleven" #BetaHA = 0.000008; BetaAA = 0.09
outputrh12$state <- "Twelve" #BetaHA = 0.000008; BetaAA = 0.08
outputrh13$state <- "Thirteen" #BetaHA = 0.000008; BetaAA = 0.07
new.data6 <- do.call("rbind", list(outputrh1, outputrh11, outputrh12, outputrh13))

ggplot(data = new.data6, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.00005), expand = c(0,0)) +
  scale_colour_manual(values = c("black","steelblue1","slateblue1","slategray"), 
                      labels = c("Baseline   ", 
                                 expression(paste(italic("r")[H]," = 130%")),
                                 expression(paste(italic("r")[H]," = 120%")),
                                 expression(paste(italic("r")[H]," = 110%")))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_hline(yintercept = 1.595806e-05, linetype = 2, size = 1) + 
  annotate("text", 0.42, 1.595806e-05, vjust = -0.6, label = "Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputrh1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputrh1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputrh11, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputrh11, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputrh12, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputrh12, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputrh13, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputrh13, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                          " per 100,000"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputrh1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputrh11, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputrh12, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputrh13, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

## Slaughter rate (ra)

outputua1$state <- "Baseline" #BetaHA = 0.000008; BetaAA = 0.1
outputua11$state <- "Eleven" #BetaHA = 0.000008; BetaAA = 0.09
outputua12$state <- "Twelve" #BetaHA = 0.000008; BetaAA = 0.08
outputua13$state <- "Thirteen" #BetaHA = 0.000008; BetaAA = 0.07
new.data9 <- do.call("rbind", list(outputua1, outputua11, outputua12, outputua13))

ggplot(data = new.data9, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.00005), expand = c(0,0)) +
  scale_colour_manual(values = c("black","steelblue1","slateblue1","slategray"), 
                      labels = c("Baseline   ", 
                                 expression(paste(mu[A]," = 130%")),
                                 expression(paste(mu[A]," = 120%")),
                                 expression(paste(mu[A]," = 110%")))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_hline(yintercept = 1.595806e-05, linetype = 2, size = 1) + 
  annotate("text", 0.42, 1.595806e-05, vjust = -0.6, label = "Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputua1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputua1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputua11, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputua11, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputua12, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputua12, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputua13, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputua13, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                          " per 100,000"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputua1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputua11, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputua12, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputua13, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))


##Theta 

outputtheta12$state <- "Twelve" #BetaHA = 0.000009; BetaAA = 0.1
outputtheta11$state <- "Eleven" #BetaHA = 0.000009; BetaAA = 0.0
outputtheta1$state <- "One" #BetaHA = 0.000009; BetaAA = 0.08
outputtheta09$state <- "Nine" #BetaHA = 0.000009; BetaAA = 0.07
outputtheta08$state <- "Eight" #BetaHA = 0.000009; BetaAA = 0.07
new.data7 <- do.call("rbind", list(outputtheta12, outputtheta11, outputtheta1, outputtheta09, outputtheta08))

ggplot(data = new.data7, aes(x = tau, y = ICombH, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.00005), expand = c(0,0)) +
  scale_colour_manual(values = c("steelblue1","slateblue","black","firebrick","red"), 
                      labels = c(expression(paste(theta," = 80%")), 
                                 expression(paste(theta," = 90%")),
                                 "Baseline", 
                                 expression(paste(theta," = 110%")), 
                                 expression(paste(theta," = 120%"))))+
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm")) +
  geom_hline(yintercept = 1.595806e-05, linetype = 2, size = 1) + 
  annotate("text", 0.42, 1.595806e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputtheta12, aes(tau,ICombH), linetype= 3, size = 1, colour = "red") + 
  geom_line(data=outputtheta12, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "red") +
  geom_line(data=outputtheta11, aes(tau,ICombH), linetype= 3, size = 1, colour = "firebrick") + 
  geom_line(data=outputtheta11, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "firebrick") +
  geom_line(data=outputtheta1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputtheta1, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "black") +
  geom_line(data=outputtheta09, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputtheta09, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputtheta08, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputtheta08, aes(tau, replace(ICombH, ICombH>1.245, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.245, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.245, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans "," (I"["CombH"],")", 
                                                                                " per 100,000"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputtheta12, aes(tau,IResRat), linetype= 1, size = 1, colour = "red") +
  geom_line(data=outputtheta11, aes(tau,IResRat), linetype= 1, size = 1, colour = "firebrick") +
  geom_line(data=outputtheta1, aes(tau,IResRat), linetype= 1, size = 1, colour = "black") + 
  geom_line(data=outputtheta09, aes(tau,IResRat), linetype= 1, size = 1, colour = "slateblue1") +
  geom_line(data=outputtheta08, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Resistant Human Infection (ResProp)"))) +
  theme(legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=13.5),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Scenario Modelling ####

#Baseline - 125% Increase in Recovery
zetabase <- 1
#Baseline - 62.5% Increase in Recovery - BetaAA = 0.075
zeta5 <- 0.5
#Baseline - 25% Increase in Recovery - BetaAA = 0.0665
zeta25 <- 0.25
#Zero Effect - No Increase in Recovery - BetaAA = 0.0565
zetazero <- 0

parmtau <- seq(0,0.1,by=0.005)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 0 , rh =  (7^-1), ua = 42^-1, uh = 28835^-1, betaAA = (0.0415), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.02, tau = parmtau[i], theta = 0.5, zeta = 1, lambda = 0.7)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
icombrat <- output1$ICombH[output1$tau == 0.5]/output1$ICombH[output1$tau == 0]
output1$IResRat <- signif(output1$IResRat, digits = 3)

output2 <- output1
output2$InfHumans <- output2$InfHumans*100000
output2$ResInfHumans <- output2$ResInfHumans*100000
output2$ICombH <- output2$ICombH*100000

table <- data.frame(x = c("Baseline Curtailment", "70% BetaAA", "70% BetaHA", "70% BetaAA/BetaHA"),
                    y = c(3.055010, 1.7933500, 2.138530, 1.2553500))
table$x <- factor(table$x, levels = c("Baseline Curtailment", "70% BetaAA", "70% BetaHA", "70% BetaAA/BetaHA"))

plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)") %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,3.5), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (Tau)"),
        legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", 
         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                            xshift =3))

#Contour plot
#Ranges for Parameter Testing

betaHArange <- seq(0.000005,0.00001, by=0.0000001)
betaAArange <- seq(0.02075,0.0415, by = 0.000415)
betaaaperc <- seq(50, 101, by = 1)
betahaperc <- seq(50,101, by = 1)
#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(betaHArange, betaAArange)
colnames(combparm1)[1:2] <- c("betaHA","betaAA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 100000, by = 500)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 0, rh =  (7^-1), ua = 42^-1, uh = 28835^-1, betaAA = combparm1[i,2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = combparm1[i,1], phi = 0.02, tau = 0, theta = 0.5, zeta = 1, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  print(out[nrow(out),7])
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- rounding(out[nrow(out),5]) 
  temp[1,4] <- rounding(out[nrow(out),6])
  temp[1,5] <- rounding(out[nrow(out),7])
  temp[1,6] <- temp[1,5]/(temp[1,4] + temp[1,5])
  temp[1,7] <- temp[1,4] + temp[1,5]
  print(temp[1,7])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

#temp[1,6] <- out[nrow(out),6] + out[nrow(out),7]
colnames(surfaceoutput1)[1:7] <- c("betaHA","betaAA","SuscHum","InfSensHum", "InfResHum", "PropRes","Comb Inf Res/Sens Humans")
plot(surfaceoutput1$InfSensHum, surfaceoutput1$InfResHum, xlim = c(-0.001,1), ylim = c(-0.001,1))

surfaceoutputplot <- surfaceoutput1[,c(1,2,7)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaAA", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$betaHA
mat$betaHA <- NULL
mat1 <- data.matrix(mat)
mat2 <- mat1
#mat2[mat2 < 100] <- 0.356e-05
mat3 <- mat1*100000
mat4 <- mat2*100000

p100 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour", 
               contours = list(start = 1.245, end = 3, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Animal-to-Animal Transmission (% of Baseline Value)"),
         yaxis = list(title = "Animal-to-Human Transmission (% of Baseline Value)"))
p100

#### Better Plots for Zeta and Lambda Testing ####

parmtau <- seq(0,0.1,by=0.001)
#parmtau <- seq(0,0.05,by=0.005)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
outputzeta15 <- data.frame()
times <- seq(0, 200000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 0, rh =  (7^-1), ua = 42^-1, uh = 28835^-1, betaAA = (0.0415), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.02, tau = parmtau[i], theta = 0.5, zeta = 1.5, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,3])
  outputzeta15 <- rbind.data.frame(outputzeta15, temp)
}

colnames(outputzeta15)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")

#icombrat <- outputAA1$ICombH[outputAA1$tau == 0.5]/output1$ICombH[outputAA1$tau == 0]
#outputAA1$IResRat <- signif(outputAA1$IResRat, digits = 3)

outputzeta15$InfHumans <- outputzeta15$InfHumans*100000
outputzeta15$ResInfHumans <- outputzeta15$ResInfHumans*100000
outputzeta15$ICombH <- outputzeta15$ICombH*100000

#
outputlambda1$state <- "LambdaBase" #BetaHA = 0.00001; BetaAA = 0.1
outputlambda09$state <- "Lambda09" #BetaHA = 0.00001; BetaAA = 0.09 
outputlambda08$state <- "Lambda08" #BetaHA = 0.00001; BetaAA = 0.08
outputlambda07$state <- "Lambda07" #BetaHA = 0.00001; BetaAA = 0.07
outputlambda06$state <- "Lambda06" #BetaHA = 0.00001; BetaAA = 0.07
outputlambda05$state <- "Lambda05" #BetaHA = 0.00001; BetaAA = 0.07
new.data1 <- do.call("rbind", list(outputlambda1, outputlambda09, outputlambda08, outputlambda07, outputlambda06,
                                   outputlambda05))
#
outputzeta1$state <- "ZetaBase" #BetaHA = 0.00001; BetaAA = 0.1
outputzeta05$state <- "Zeta05" #BetaHA = 0.00001; BetaAA = 0.09 
outputzeta025$state <- "Zeta025" #BetaHA = 0.00001; BetaAA = 0.08
outputzeta0$state <- "Zeta0" #BetaHA = 0.00001; BetaAA = 0.07
outputzeta15$state <- "Zeta15" #BetaHA = 0.00001; BetaAA = 0.07
new.data <- do.call("rbind", list(outputzeta1, outputzeta05, outputzeta025, outputzeta0, outputzeta15))
#

# new lambda with 05, 04
ggplot(data = new.data1, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1.2) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = "Proportion of Resistant Human Infection (ResProp)") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_colour_manual(values = c("gold","salmon","firebrick2","indianred4","grey40","black"), 
                      labels = c(expression(paste(lambda," = 50%")),
                                 expression(paste(lambda," = 60%")),
                                 expression(paste(lambda," = 70%")), 
                                 expression(paste(lambda," = 80%")), 
                                 expression(paste(lambda," = 90%")),
                                 "Baseline   "))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

#"Proportion of Infected Humans (ICombH) per 100,000"
ggplot() + 
  geom_line(data=outputlambda1, aes(tau,ICombH), size = 1.1, colour = "black") + 
  geom_line(data=outputlambda09, aes(tau,ICombH), size = 1.1, colour = "grey40") + 
  geom_line(data=outputlambda08, aes(tau,ICombH), size = 1.1, colour = "indianred4") + 
  geom_line(data=outputlambda07, aes(tau,ICombH), size = 1.1, colour = "firebrick2") +
  geom_line(data=outputlambda06, aes(tau,ICombH), size = 1.1, colour = "salmon") +
  geom_line(data=outputlambda05, aes(tau,ICombH), size = 1.1, colour = "gold") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans ","(I"["CombH"],")", " per 100,000"))) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

#

ggplot(data = new.data, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1.2) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = "Proportion of Resistant Human Infection (ResProp)") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_colour_manual(values = c("gold","salmon","indianred4","grey50","black"), 
                      labels = c(expression(paste(zeta," = 0%")),
                                 expression(paste(zeta," = 25%")),
                                 expression(paste(zeta," = 50%")), 
                                 expression(paste(zeta," = 150%")),
                                 "Baseline   "))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))

ggplot() + 
  geom_line(data=outputzeta1, aes(tau,ICombH), size = 1.1, colour = "black") + 
  geom_line(data=outputzeta05, aes(tau,ICombH), size = 1.1, colour = "indianred4") + 
  geom_line(data=outputzeta025, aes(tau,ICombH), size = 1.1, colour = "salmon") + 
  geom_line(data=outputzeta0, aes(tau,ICombH), size = 1.1, colour = "gold") +
  geom_line(data=outputzeta15, aes(tau,ICombH), size = 1.1, colour = "grey50") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3.5), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), y = expression(paste("Proportion of Infected Humans ","(I"["CombH"],")", " per 100,000"))) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))


#### Contour plot Testing ####
#Ranges for Parameter Testing

betaHArange <- seq(0.000005,0.00001, by=0.0000001)
betaAArange <- seq(0.02075,0.0415, by = 0.000415)
betaaaperc <- seq(50, 101, by = 1)
betahaperc <- seq(50,101, by = 1)
#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(betaHArange, betaAArange)
colnames(combparm1)[1:2] <- c("betaHA","betaAA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 100000, by = 500)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 0, rh =  (7^-1), ua = 56^-1, uh = 28835^-1, betaAA = combparm1[i,2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = combparm1[i,1], phi = 0.02, tau = 0, theta = 0.5, zeta = zetabase, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  print(out[nrow(out),7])
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- rounding(out[nrow(out),5]) 
  temp[1,4] <- rounding(out[nrow(out),6])
  temp[1,5] <- rounding(out[nrow(out),7])
  temp[1,6] <- temp[1,5]/(temp[1,4] + temp[1,5])
  temp[1,7] <- temp[1,4] + temp[1,5]
  print(temp[1,7])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

#temp[1,6] <- out[nrow(out),6] + out[nrow(out),7]
colnames(surfaceoutput1)[1:7] <- c("betaHA","betaAA","SuscHum","InfSensHum", "InfResHum", "PropRes","Comb Inf Res/Sens Humans")
plot(surfaceoutput1$InfSensHum, surfaceoutput1$InfResHum, xlim = c(-0.001,1), ylim = c(-0.001,1))

surfaceoutputplot <- surfaceoutput1[,c(1,2,7)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaAA", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$betaHA
mat$betaHA <- NULL
mat1 <- data.matrix(mat)
mat2 <- mat1
#mat2[mat2 < 100] <- 0.356e-05
mat3 <- mat1*100000
mat4 <- mat2*100000

p100 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
                contours = list(start = 1.2446, end = 4, size = 0.2, showlabels = TRUE,
                                labelfont = list(size = 12, color = 'white')),
                colorbar = list(title = "I<sub>CombH</sub>", 
                                thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Animal-to-Animal Transmission (% of Baseline Value)"),
         yaxis = list(title = "Animal-to-Human Transmission (% of Baseline Value)"))
p100

p15 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
               contours = list(start = 0.356, end = 3.05501, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Animal-to-Animal Transmission (%)"),
         yaxis = list(title = "Animal-to-Human Transmission (%)"))
p15

p16 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
               contours = list(start = 2.0788, end = 3.05501, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "BetaAA (%)"),
         yaxis = list(title = "BetaHA (%)"))
p16


p17 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
               contours = list(start = 2.53, end = 3.05501, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "BetaAA (%)"),
         yaxis = list(title = "BetaHA (%)"))
p17

p18 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
               contours = list(start = 3.1, end = 3.1, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "BetaAA (%)"),
         yaxis = list(title = "BetaHA (%)"))
p18
#### Surface Plots for Tau vs BAA RElationship - Supplementary ####

taurange <- seq(0,0.1, by = 0.0005)
betaAAsupprange <- seq(0,0.08, by = 0.001)
#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(taurange, betaAAsupprange)
colnames(combparm1)[1:2] <- c("tau","betaAA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 100000, by = 500)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 0, rh =  (7^-1), ua = 84^-1, uh = 28835^-1, betaAA = combparm1[i,2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.00001, phi = 0.02, tau = combparm1[i,1], theta = 0.5, zeta = 1, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  print(out[nrow(out),7])
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- rounding(out[nrow(out),5]) 
  temp[1,4] <- rounding(out[nrow(out),6])
  temp[1,5] <- rounding(out[nrow(out),7])
  temp[1,6] <- temp[1,5]/(temp[1,4] + temp[1,5])
  temp[1,7] <- temp[1,4] + temp[1,5]
  print(temp[1,7])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

#temp[1,6] <- out[nrow(out),6] + out[nrow(out),7]
colnames(surfaceoutput1)[1:7] <- c("tau","betaAA","SuscHum","InfSensHum", "InfResHum", "PropRes","Comb Inf Res/Sens Humans")

surfaceoutputploticomb <- surfaceoutput1[,c(1,2,7)]
surfaceoutputploticomb[3] <- surfaceoutputploticomb[3] * 100000
surfaceoutputplotres <- surfaceoutput1[,c(1,2,6)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputploticomb$new <- NULL

mat <- spread(surfaceoutputploticomb, key = "betaAA", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)
mat2 <- t(mat1)

plot_ly(x = taurange, y = betaAAsupprange , z = mat2, type = "contour",
        contours = list(start = 0, end = 6, size = 0.6),
        colorbar = list(title = "I<sub>CombH</sub> ", 
                        thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Livestock Antibiotic Usage (Tau)"),
         yaxis = list(title = "Animal-to-Human Transmission (BetaAA)"))

surfaceoutputplotres$new <- NULL
mat4 <- spread(surfaceoutputplotres, key = "betaAA", value = "PropRes")
row.names(mat4) <- mat4$tau
mat4$tau <- NULL
mat5 <- data.matrix(mat4)
mat6 <- t(mat5)

plot_ly(x = taurange, y = betaAAsupprange , z = mat6, type = "contour",
        contours = list(start = 0, end = 1, size = 0.1),
        colorbar = list(title = "ResProp", 
                        thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Livestock Antibiotic Usage (Tau)"),
         yaxis = list(title = "Animal-to-Human Transmission (BetaAA)"))

#### Relationship between Livestock Tau and Livestock Resistance ####

parmtau <- seq(0,0.1,by=0.001)
#parmtau <- seq(0,0.05,by=0.005)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
outputanres <- data.frame()
times <- seq(0, 200000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms2 = c(ra = 0, rh =  (7^-1), ua = 42^-1, uh = 28835^-1, betaAA = (0.0415), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.02, tau = parmtau[i], theta = 0.5, zeta = 1, lambda = 1)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  temp[1,7] <- rounding(out[nrow(out),4])/(rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])) 
  print(temp[1,7])
  outputanres <- rbind.data.frame(outputanres, temp)
}

colnames(outputanres)[1:7] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat", "AniResRat")

outputanres$InfHumans <- outputanres$InfHumans*100000
outputanres$ResInfHumans <- outputanres$ResInfHumans*100000
outputanres$ICombH <- outputanres$ICombH*100000

ggplot() + 
  geom_line(data=outputanres, aes(tau,AniResRat), linetype= 1, size = 2, colour = "black") +
  scale_x_continuous(limits = c(0,0.1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Livestock Antibiotic Usage (",tau,")")), 
       y = expression(paste("Proportion of Resistant Animal Infection (ResProp)"))) +
  theme(legend.title = element_blank(), text = element_text(size=13),plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

