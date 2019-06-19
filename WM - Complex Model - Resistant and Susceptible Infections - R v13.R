setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Mathematical Models")

rm(list=ls())

library("deSolve")
library("fast")
library("sensitivity")
library("ggplot2")
library("plotly")
library("tidyr")
library("nlmeODE")
library("phaseR")
library("reshape2")

#### Model Functions + Output ####
amr <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSa = ua + ra*(Ia + Ira) + tau*Ia - (betaAA*Ia*Sa) - (betaAH*Ih*Sa) - (betaAH*Irh*Sa) - (betaAA*Ira*Sa) - ua*Sa  
    dIa = betaAA*Ia*Sa + betaAH*Ih*Sa + phi*Ira - tau*Ia - tau*theta*Ia - ra*Ia - ua*Ia
    dIra = betaAH*Irh*Sa + betaAA*Ira*Sa + tau*theta*Ia - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ih+Irh) - (betaHH*Ih*Sh) - (betaHH*Irh*Sh) - (betaHA*Ia*Sh) - (betaHA*Ira*Sh) - uh*Sh 
    dIh = betaHH*Ih*Sh + betaHA*Ia*Sh - rh*Ih - uh*Ih 
    dIrh = betaHH*Irh*Sh + betaHA*Ira*Sh - rh*Irh - uh*Irh 
    return(list(c(dSa,dIa,dIra,dSh,dIh,dIrh)))
  })
}

rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Model Testbed - Basic Model Output ####

init <- c(Sa=0.98, Ia=0.01, Ira=0.01, Sh=1, Ih=0, Irh=0)
times1 <- seq(0,1000,by=0.1)

#Need to Specify Model Parameters
parms = c(ra = 52^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.10, betaAH = 0.000001, betaHH = 0.000001, 
          betaHA = 0.00001, phi = 0.1, tau = 0.1, theta = 0.5)
parms = c(ra = 52^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.07, betaAH = 0.000001, betaHH = 0.000001, 
          betaHA = 0.00001, phi = 0.1, tau = 0.1, theta = 0.5)
parms = c(ra = 22^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.15, betaAH = 0.000001, betaHH = 0.000001, 
          betaHA = 0.00001, phi = 0.1, tau = 0.1, theta = 0.5)
parms = c(ra = 25^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.09, betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = 0.00001, phi = 0.05, tau = 0.05, theta = 0.5)

out <- ode(y = init, func = amr, times = times1, parms = parms)

#GGPLOT FIGURE
outdata <- data.frame(out)
outdata$combIa <- as.numeric(outdata$Ia) + as.numeric(outdata$Ira)
outdata$combIh <- as.numeric(outdata$Ih) + as.numeric(outdata$Irh)

#Animal Population 
par(mfrow=c(1,2))
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sa, colour = "Sa"), size = 1.1) +
  geom_line(aes(y = Ia, colour = "Ia"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Ira, colour = "Ira"), size = 1.1) +
  geom_line(aes(y = combIa, colour = "combIa"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sa" = "chartreuse3", "Ia" = "red", "Ira" = "blue", "combIa" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombA"],")")),
                                expression(paste("Resistant Infection ","(I"["RA"],")")),
                                expression(paste("Sensitive Infection ","(I"["A"],")")),
                                expression(paste("Susceptible ","(S"["A"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Human Population 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sh, colour = "Sh"), size = 1.1) +
  geom_line(aes(y = Ih, colour = "Ih"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Irh, colour = "Irh"), size = 1.1) +
  geom_line(aes(y = combIh, colour = "combIh"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sh" = "chartreuse3", "Ih" = "red", "Irh" = "blue", "combIh" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombH"],")")),
                                expression(paste("Resistant Infection ","(I"["RH"],")")),
                                expression(paste("Sensitive Infection ","(I"["H"],")")),
                                expression(paste("Susceptible ","(S"["H"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#High Resolution Human Population 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Ih, colour = "Ih"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Irh, colour = "Irh"), size = 1.1) +
  geom_line(aes(y = combIh, colour = "combIh"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Ih" = "red", "Irh" = "blue", "combIh" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombH"],")")),
                                expression(paste("Resistant Infection ","(I"["RH"],")")),
                                expression(paste("Sensitive Infection ","(I"["H"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000025), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

Icomb <- as.numeric(out[nrow(out),6]) + as.numeric(out[nrow(out),7])  
print(Icomb)
#plot(out[,"time"], out[,"Ia"]+ out[,"Ira"], type="l", lwd=2, col="black", ylab="IComb*", xlab="Time")

#### For 100000

init <- c(Sa=0.98, Ia=0.01, Ira=0.01, Sh=1, Ih=0, Irh=0)
times1 <- seq(0,1000,by=0.1)

#Need to Specify Model Parameters

parms = c(ra = 25^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.09, betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = 0.00001, phi = 0.05, tau = 0.05, theta = 0.5)

out <- ode(y = init, func = amr, times = times1, parms = parms)

#GGPLOT FIGURE
outdata <- data.frame(out)
outdata$Ih <- outdata$Ih*100000
outdata$Irh <- outdata$Irh*100000
outdata$combIa <- as.numeric(outdata$Ia) + as.numeric(outdata$Ira) 
outdata$combIh <- as.numeric(outdata$Ih) + as.numeric(outdata$Irh)

#Animal Population 
par(mfrow=c(1,2))
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sa, colour = "Sa"), size = 1.1) +
  geom_line(aes(y = Ia, colour = "Ia"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Ira, colour = "Ira"), size = 1.1) +
  geom_line(aes(y = combIa, colour = "combIa"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sa" = "chartreuse3", "Ia" = "red", "Ira" = "blue", "combIa" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombA"],")")),
                                expression(paste("Resistant Infection ","(I"["RA"],")")),
                                expression(paste("Sensitive Infection ","(I"["A"],")")),
                                expression(paste("Susceptible ","(S"["A"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Human Population 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sh, colour = "Sh"), size = 1.1) +
  geom_line(aes(y = Ih, colour = "Ih"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Irh, colour = "Irh"), size = 1.1) +
  geom_line(aes(y = combIh, colour = "combIh"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sh" = "chartreuse3", "Ih" = "red", "Irh" = "blue", "combIh" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombH"],")")),
                                expression(paste("Resistant Infection ","(I"["RH"],")")),
                                expression(paste("Sensitive Infection ","(I"["H"],")")),
                                expression(paste("Susceptible ","(S"["H"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#High Resolution Human Population 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Ih, colour = "Ih"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Irh, colour = "Irh"), size = 1.1) +
  geom_line(aes(y = combIh, colour = "combIh"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Ih" = "red", "Irh" = "blue", "combIh" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombH"],")")),
                                expression(paste("Resistant Infection ","(I"["RH"],")")),
                                expression(paste("Sensitive Infection ","(I"["H"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Function for Parameter Combinations - From FAST ####
start_time <- Sys.time()

parms = fast_parameters(minimum = c(2.5^-1,0.6^-1,24^-1,2883.5^-1,0,0,0,0,0,0,0), maximum = c(250^-1,60^-1,2400^-1,288350^-1,0.9,0.0001,0,0.001,0.1,0.1,5), 
                        factor=11, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "tau", "theta"))

#### Creating Model Output for Parameter Combinations from FAST Parameters #### 

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0,200000, by = 10) 

output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=8))
  parms1 = c(ra = parms$ra[i], rh = parms$rh[i] , ua = parms$ua[i], uh = parms$uh[i], betaAA = parms$betaAA[i],
             betaAH = parms$betaAH[i], betaHH = parms$betaHH[i], betaHA = parms$betaHA[i], phi=parms$phi[i],
             tau=parms$tau[i], theta=parms$theta[i])
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  temp[1,1] <- rounding(out[nrow(out),5]) 
  temp[1,2] <- rounding(out[nrow(out),6]) 
  temp[1,3] <- rounding(out[nrow(out),7])
  temp[1,4] <- temp[1,2] + temp[1,3]
  temp[1,5] <- temp[1,3]/temp[1,4]
  temp[1,6] <- signif(((rounding(out[nrow(out),6]) + rounding(out[nrow(out)-1,6]) + rounding(out[nrow(out)-2,6]))/3), digits = 6)
  if(temp[1,6] == temp[1,2]) {temp[1,7] <- "YES"}
  temp[1,8] <- rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]) 
  print(temp[1,3])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:8] <- c("SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat", "ValueatEqui","Equi?","ICombA")

end_time <- Sys.time()
end_time - start_time


#### Sensitivity Analysis ####
#For iCombH

sensit <- output$ICombH #Creating Variable for the output variable of interest
sens<-sensitivity(x=sensit, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                               "phi", "tau", "theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi","tau", "theta"), value=sens)

p <- ggplot(df.equilibrium, aes(parameter, value))
p + geom_bar(stat="identity", fill="grey23")

#For ResRat

test <- output[complete.cases(output),] # to get rid of NA found in the dataframe - make it analysable int he sensitivity analysis

sensit1 <- test$IResRat #Creating Variable for the output variable of interest
sens1 <-sensitivity(x=sensit1, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                  "phi", "tau", "theta"))

df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi","tau", "theta"), value=sens1)

#Combined Bar Plot
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
  scale_y_continuous(limits = c(0,0.25), expand = c(0,0.00009)) + 
  scale_fill_manual(values=c("dodgerblue2", "orangered1")) +
  theme(axis.text=element_text(size = 12, colour = "black"),
        axis.line.x = element_line(color="black", size = 0.7),
        axis.title=element_text(size=12), 
        legend.position=c(0.15, 0.9), legend.title = element_blank()) +
  ylab("Partial Variance") 

#### Comb IRH and IH Measure Testing - Effect of Treatment ####

#Testing for the Effect of Changing Treatment on the Combined Measure

parmtau <- seq(0,0.5,by=0.02)
parmtau <- seq(0,0.2,by=0.01)
#parmtau <- seq(0,0.05,by=0.005)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 25^-1, rh =  (6^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.09)*0.7, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001)*0.7, phi = 0.05, tau = parmtau[i], theta = 0.5)
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

plot_ly(
  data=table,
  x = ~x,
  y = ~y,
  type = "bar",
  marker = list(color = c('rgb(128,128,128)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(222,45,38,0.8)'))) %>%
  layout(yaxis = list(title = "Proportion of Infected Humans (per 100,000)", 
                      titlefont = list(size = 17) ,tickfont = list(size = 16),range = c(-0.01,3.5)),
         xaxis = list(title = "Curtailment Scenarios (Tau = 0)", tickfont = list(size = 16),  titlefont = list(size = 17)),
         margin = list(l = 30, r = 30, b = 30, t = 30, pad = 5))

#plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Sensitive Infection (Human)") %>%
#  add_trace(y= ~ResInfHumans, name = "Resistant Infection (Human)") %>% 
#  layout(yaxis = list(title = "Proportion of Infected Humans per 100,000", range = c(0,3.5), showline = TRUE),
#         xaxis = list(title = "Livestock Antibiotic Usage"),
#        legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
#         barmode = "stack", 
#         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
#                            xshift =3))

plot_ly(output1, x= ~tau, y = ~InfHumans, type = "bar", name = "Sensitive Infection") %>%
  add_trace(y= ~ResInfHumans, name = "Resistant Infection") %>% 
  layout(yaxis = list(title = "Proportion of Humans Infected", exponentformat= "E", range = c(0,3.5E-5), 
                      titlefont = list(size = 18) ,tickfont = list(size = 16)),
         xaxis = list(title = "Tau (Antibiotic Usage)", range = c(-0.01,0.2),
                      titlefont = list(size = 18), tickfont = list(size = 16)),
         margin = list(l = 30, r = 30, b = 30, t = 30, pad = 5),
         legend = list(orientation = "h", x = 0.70, y=1, font= list(size = 16)), showlegend = T,
         barmode = "stack",
         annotations = list(
           list(x = 0.18, y = 1.45e-05, text = "Baseline Parameters", yanchor = "bottom",
                showarrow = FALSE, xshift= -30, yshift = 3, font = list(size = 16)),
           list(x = 0.18, y = output1$ICombH[6], text = "Intervention w/ Antibiotic Usage", yanchor = "bottom",
                showarrow = FALSE, xshift= -75, yshift = 3, font = list(size = 16, color = "red")),
           list(ax = 0.05, ay = 1.45e-05, x = 0.05, y = output1$ICombH[6],
                axref = "x", ayref = "y", xref = "x", yref = "y", showarrow= TRUE, arrowhead=1, 
                arrowsize = 1.2, arrowwidth=2.5, arrowcolor= "red"))) %>%
  add_segments(x=-0.01, xend = 0.200001, y=1.388e-05, yend = 1.388e-05, line = list(color = "black", dash = "dot", width = 2),
               showlegend = FALSE) %>%
  add_segments(x=-0.01, xend = 0.200001, y=output1$ICombH[6], yend = output1$ICombH[6], 
               line = list(color = "red", dash = "dot", width = 3), showlegend = FALSE)

#Have put 6 since that is where Tau is equal to 0.05

#### Comb IRH and IH Measure Testing - Effect of Treatment - Just below 0.1####

#Testing for the Effect of Changing Treatment on the Combined Measure

parmtau <- seq(0,0.05,by=0.002)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 25^-1, rh =  (6^-1)*1.3, ua = 240^-1, uh = 28835^-1, betaAA = (0.09), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001)*0.7, phi = (0.05), tau = parmtau[i], theta = 0.5)
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

#plot_ly(output1, x= ~tau, y = ~InfHumans, type = "bar", name = "Sensitive Infection (Human)") %>%
#  add_trace(y= ~ResInfHumans, name = "Resistant Infection (Human)") %>% 
#  layout(yaxis = list(title = "Proportion Infected", exponentformat= "E", range = c(0,5E-5), showline = TRUE),
#         xaxis = list(title = "Tau (Antibiotic Usage)"),
#        legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
#         barmode = "stack", 
#         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
#                            xshift =3))

plot_ly(output1, x= ~tau, y = ~InfHumans, type = "bar", name = "Sensitive Infection") %>%
  add_trace(y= ~ResInfHumans, name = "Resistant Infection") %>% 
  layout(yaxis = list(title = "Proportion of Humans Infected", exponentformat= "E", range = c(0,3E-5), 
                      titlefont = list(size = 18) ,tickfont = list(size = 18)),
         xaxis = list(title = "Tau (Antibiotic Usage)", range = c(-0.001,0.05),
                      titlefont = list(size = 18), tickfont = list(size = 18)),
         margin = list(l = 30, r = 60, b = 30, t = 30, pad = 5),
         legend = list(orientation = "h", x = 0.70, y=1, font= list(size = 18)), showlegend = F,
         barmode = "stack",
         annotations = list(
           list(x = 0.045, y = 1.388e-05, text = "Baseline Parameters", yanchor = "bottom",
                showarrow = FALSE, xshift= -30, yshift = 3, font = list(size = 18)),
           list(ax = 0.05, ay = 1.45e-05, x = 0.05, y = output1$ICombH[26],
                axref = "x", ayref = "y", xref = "x", yref = "y", showarrow= TRUE, arrowhead=1, 
                arrowsize = 1.2, arrowwidth=2.5, arrowcolor= "red"))) %>%
  add_segments(x=-0.01, xend = 0.051, y=1.388e-05, yend = 1.388e-05, line = list(color = "black", dash = "dot", width = 2),
               showlegend = FALSE) %>%
  add_segments(x=-0.01, xend = 0.051, y=output1$ICombH[26], yend = output1$ICombH[26], 
               line = list(color = "red", dash = "dot", width = 3), showlegend = FALSE)

#### Comb IRH and IH Measure Testing - Effect of Treatment - NEW PLOT ####

#Testing for the Effect of Changing Treatment on the Combined Measure

#parmtau <- seq(0,0.5,by=0.03)
parmtau <- seq(0,0.2,by=0.0002)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

#Have to alter a number of the outputnames to generate all the values for the figures
outputtheta08 <- data.frame()

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 25^-1, rh =  (6^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.09), betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.00001, phi = (0.05), tau = parmtau[i], theta = (0.5)*0.8)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,3])
  outputtheta08 <- rbind.data.frame(outputtheta08, temp)
}

colnames(outputtheta08)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")

#### Plotting ####
outputAA1$state <- "Baseline" #BetaHA = 0.00001; BetaAA = 0.1
outputAA9$state <- "BetaHA = 100%; BetaAA = 90%" #BetaHA = 0.00001; BetaAA = 0.09 
outputAA8$state <- "BetaHA = 100%; BetaAA = 80%" #BetaHA = 0.00001; BetaAA = 0.08
outputAA7$state <- "BetaHA = 100%; BetaAA = 70%" #BetaHA = 0.00001; BetaAA = 0.07
new.data <- do.call("rbind", list(outputAA1, outputAA9, outputAA8, outputAA7))

ggplot(data = new.data, aes(x = tau, y = IResRat, colour = state)) +
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
  geom_hline(yintercept = 1.595806e-05, linetype = 2, size = 1) + 
  annotate("text", 0.42, 1.595806e-05, vjust = -0.6, label = "Current Antibiotic Usage")

ggplot() + 
  geom_line(data=outputAA1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputAA1, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "black") +
  geom_line(data=outputAA9, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputAA9, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputAA8, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputAA8, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputAA7, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputAA7, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.38806e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  theme(axis.line.x = element_line(color="black", size = 1), text = element_text(size=14),
        legend.position=c(0.8, 0.85), legend.title = element_blank(),
        legend.spacing.x = unit(0.7, 'cm'), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputAA1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputAA9, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputAA8, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputAA7, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Resistance Ratio (ResRat)"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##

outputHA1$state <- "Baseline Parameters" #BetaHA = 0.00001; BetaAA = 0.1
outputHA9$state <- "Beta" #BetaHA = 0.000009; BetaAA = 0.1
outputHA8$state <- "Eight" #BetaHA = 0.000008; BetaAA = 0.1
outputHA7$state <- "Seven" #BetaHA = 0.000007; BetaAA = 0.1
new.data1 <- do.call("rbind", list(outputHA1, outputHA9, outputHA8, outputHA7))

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
  geom_line(data=outputHA1, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "black") +
  geom_line(data=outputHA9, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputHA9, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputHA8, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputHA8, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputHA7, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputHA7, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.388e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputHA1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputHA9, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputHA8, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputHA7, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Resistance Ratio (ResRat)"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))
##

outputphi12$state <- "Twelve" #BetaHA = 0.000009; BetaAA = 0.1
outputphi11$state <- "Eleven" #BetaHA = 0.000009; BetaAA = 0.0
outputphi1$state <- "One" #BetaHA = 0.000009; BetaAA = 0.08
outputphi09$state <- "Nine" #BetaHA = 0.000009; BetaAA = 0.07
outputphi08$state <- "Eight" #BetaHA = 0.000009; BetaAA = 0.07
new.data2 <- do.call("rbind", list(outputphi12, outputphi11, outputphi1, outputphi09, outputphi08))

ggplot(data = new.data2, aes(x = tau, y = IResRat, colour = state)) +
  geom_line(size = 1) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
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
  geom_line(data=outputphi12, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "red") +
  geom_line(data=outputphi11, aes(tau,ICombH), linetype= 3, size = 1, colour = "firebrick") + 
  geom_line(data=outputphi11, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "firebrick") +
  geom_line(data=outputphi1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputphi1, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "black") +
  geom_line(data=outputphi09, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputphi09, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputphi08, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputphi08, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.388e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputphi12, aes(tau,IResRat), linetype= 1, size = 1, colour = "red") +
  geom_line(data=outputphi11, aes(tau,IResRat), linetype= 1, size = 1, colour = "firebrick") +
  geom_line(data=outputphi1, aes(tau,IResRat), linetype= 1, size = 1, colour = "black") + 
  geom_line(data=outputphi09, aes(tau,IResRat), linetype= 1, size = 1, colour = "slateblue1") +
  geom_line(data=outputphi08, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Resistance Ratio (ResRat)"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##

outputrh1$state <- "Baseline" #BetaHA = 0.000008; BetaAA = 0.1
outputrh11$state <- "Eleven" #BetaHA = 0.000008; BetaAA = 0.09
outputrh12$state <- "Twelve" #BetaHA = 0.000008; BetaAA = 0.08
outputrh13$state <- "Thirteen" #BetaHA = 0.000008; BetaAA = 0.07
new.data3 <- do.call("rbind", list(outputrh1, outputrh11, outputrh12, outputrh13))

ggplot(data = new.data1, aes(x = tau, y = IResRat, colour = state)) +
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
  geom_line(data=outputrh1, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "black") +
  geom_line(data=outputrh11, aes(tau,ICombH), linetype= 3, size = 1, colour = "slategray") + 
  geom_line(data=outputrh11, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slategray") +
  geom_line(data=outputrh12, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputrh12, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputrh13, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputrh13, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.388e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputrh1, aes(tau,IResRat), linetype= 1, size = 5, colour = "black") +
  geom_line(data=outputrh11, aes(tau,IResRat), linetype= 1, size = 3, colour = "slategray") +
  geom_line(data=outputrh12, aes(tau,IResRat), linetype= 1, size = 2, colour = "slateblue1") + 
  geom_line(data=outputrh13, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") + 
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Resistance Ratio (ResRat)"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

##Theta 

outputtheta12$state <- "Twelve" #BetaHA = 0.000009; BetaAA = 0.1
outputtheta11$state <- "Eleven" #BetaHA = 0.000009; BetaAA = 0.0
outputtheta1$state <- "One" #BetaHA = 0.000009; BetaAA = 0.08
outputtheta09$state <- "Nine" #BetaHA = 0.000009; BetaAA = 0.07
outputtheta08$state <- "Eight" #BetaHA = 0.000009; BetaAA = 0.07
new.data5 <- do.call("rbind", list(outputtheta12, outputtheta11, outputtheta1, outputtheta09, outputtheta08))

ggplot(data = new.data5, aes(x = tau, y = ICombH, colour = state)) +
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
  geom_line(data=outputtheta12, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "red") +
  geom_line(data=outputtheta11, aes(tau,ICombH), linetype= 3, size = 1, colour = "firebrick") + 
  geom_line(data=outputtheta11, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "firebrick") +
  geom_line(data=outputtheta1, aes(tau,ICombH), linetype= 3, size = 1, colour = "black") + 
  geom_line(data=outputtheta1, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "black") +
  geom_line(data=outputtheta09, aes(tau,ICombH), linetype= 3, size = 1, colour = "slateblue1") + 
  geom_line(data=outputtheta09, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "slateblue1") +
  geom_line(data=outputtheta08, aes(tau,ICombH), linetype= 3, size = 1, colour = "steelblue1") + 
  geom_line(data=outputtheta08, aes(tau, replace(ICombH, ICombH>1.388e-05, NA)), size = 1, colour = "steelblue1")+
  geom_hline(yintercept = 1.388e-05, linetype = 2, size = 1) + 
  annotate("text", 0.153, 1.388e-05, vjust = -0.6, label = "ICombH at Current Antibiotic Usage", size = 4.7) +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000035), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Overall Infection ","(I"["CombH"],")"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

ggplot() + 
  geom_line(data=outputtheta12, aes(tau,IResRat), linetype= 1, size = 1, colour = "red") +
  geom_line(data=outputtheta11, aes(tau,IResRat), linetype= 1, size = 1, colour = "firebrick") +
  geom_line(data=outputtheta1, aes(tau,IResRat), linetype= 1, size = 1, colour = "black") + 
  geom_line(data=outputtheta09, aes(tau,IResRat), linetype= 1, size = 1, colour = "slateblue1") +
  geom_line(data=outputtheta08, aes(tau,IResRat), linetype= 1, size = 1, colour = "steelblue1") +
  scale_x_continuous(limits = c(0,0.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x =expression(paste("Antibiotic Usage (",tau,")")), y = expression(paste("Resistance Ratio (ResRat)"))) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position=c(0.8, 0.85), legend.title = element_blank(), text = element_text(size=14),
        legend.spacing.x = unit(0.7, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))


#### Comb IRH and IH Measure Testing - Effect of Treatment - More Elegant Combinatitory Plot ####

#Testing for the Effect of Changing Treatment on the Combined Measure
#parmtau <- seq(0,0.5,by=0.03)
parmtau <- seq(0,0.5,by=0.01)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
output2 <- data.frame()
output3 <- data.frame()
output4 <- data.frame()
output5 <- data.frame()

times <- seq(0, 200000, by = 100)

parmvals <- c(1,0.9,0.8,0.7,0.6)

for (j in 1:length(parmvals)) {
  tempout <- data.frame()
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
    parms2 = c(ra = 52^-1, rh =  (6^-1), ua = 240^-1, uh = 28835^-1, betaAA = parmvals[j], betaAH = 0.000001, betaHH = 0.000001, 
               betaHA = 0.00001, phi = 0.1, tau = parmtau[i], theta = 0.5)
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[i]
    temp[1,2] <- parmvals[j]
    temp[1,3] <- rounding(out[nrow(out),5]) 
    temp[1,4] <- rounding(out[nrow(out),6]) 
    temp[1,5] <- rounding(out[nrow(out),7])
    temp[1,6] <- temp[1,3] + temp[1,4]
    temp[1,7] <- temp[1,4]/temp[1,5]
    print(temp[1,6])  
    tempout <- rbind.data.frame(tempout, temp)
  }
  output1 <- rbind.data.frame(output1, tempout$X6[tempout[1,2]==1])
  output2 <- rbind.data.frame(output2, tempout[tempout[,2]==0.9])
  output3 <- rbind.data.frame(output3, tempout[tempout[,2]==0.8])
  output4 <- rbind.data.frame(output4, tempout[tempout[,2]==0.7])
  output4 <- rbind.data.frame(output4, tempout[tempout[,2]==0.6])
}

colnames(output1)[1:6] <- c("tau", "betaAA", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
icombrat <- output1$ICombH[output1$tau == 0.5]/output1$ICombH[output1$tau == 0]
output1$IResRat <- signif(output1$IResRat, digits = 3)

plot(output1$tau, output1$ICombH, type = "l", lwd = 3)
lines()

plot(output1$tau, output1$IResRat, type = "l", lwd = 3)
lines()



#p10 <- plot_ly(output1, x= ~tau, y = ~InfHumans, type = "bar", name = "Sens Inf Humans") %>%
#  add_trace(y= ~ResInfHumans, name = "Res Inf Humans") %>% 
#  layout(yaxis = list(title = "Proportion Infected", exponentformat= "E", range = c(0,5E-5), showline = TRUE),
#         xaxis = list(title = "Tau (Antibiotic Usage)"),
#         legend = list(orientation = "v", x = 1.0, y=0.5), showlegend = T,
#         barmode = "stack", 
#         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
#                            xshift =3))
#p10

plot_ly(output1, x= ~tau, y = ~InfHumans, type = "bar", name = "Sens Inf Humans") %>%
  add_trace(y= ~ResInfHumans, name = "Res Inf Humans") %>% 
  layout(yaxis = list(title = "Proportion Infected", exponentformat= "E", range = c(0,5E-5), showline = TRUE),
         xaxis = list(title = "Tau (Antibiotic Usage)", range = c(-0.01,0.5)),
         legend = list(orientation = "v", x = 1.0, y=0.5), showlegend = T,
         barmode = "stack",
         annotations = list(
           list(x = 0.45, y = 1.6e-05, text = "Baseline Level of FB Disease", yanchor = "bottom",
                showarrow = FALSE, xshift= -55, yshift = 3, font = list(size = 15)),
           list(x = 0.45, y = output1$ICombH[11], text = "New Level of FB Disease", yanchor = "bottom",
                showarrow = FALSE, xshift= -55, yshift = 3, font = list(size = 15, color = "red")),
           list(ax = 0.1, ay = 1.69e-05, x = 0.1, y = output1$ICombH[11],
                axref = "x", ayref = "y", xref = "x", yref = "y", showarrow= TRUE, arrowhead=1, 
                arrowsize = 1.2, arrowwidth=2.5, arrowcolor= "red"))) %>%
  add_segments(x=-0.01, xend = 0.500001, y=1.6e-05, yend = 1.6e-05, line = list(color = "black", dash = "dot", width = 2),
               showlegend = FALSE) %>%
  add_segments(x=-0.01, xend = 0.500001, y=output1$ICombH[11], yend = output1$ICombH[11], 
               line = list(color = "red", dash = "dot", width = 3), showlegend = FALSE)

#Have put 6 since that is where Tau is equal to 0.05

#### Curve Fitting for Treatment Analysis ####

#Obtaining data for use in model fitting

parmtau <- seq(0,0.5,by=0.01)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

#These parameters are obtained from the FAST parameter analysis

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(ra = 52^-1, rh =  6^-1, ua = 240^-1, uh = 28835^-1, betaAA = 0.1, betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = 0.1, tau = parmtau[i], theta = 0.5)
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
output1$IResRat <- signif(output1$IResRat, digits = 3)

#Curve Fitting - ICombH
plot_ly(output1, x= ~tau, y = ~ICombH, type = "scatter")

ggplot(output1, aes(x=tau, y=ICombH)) %>%  + # tau is from the output1 dataframe
  geom_point() + 
  stat_smooth(method="nls", formula = y ~ SSasymp(x,Asym, R0, lrc), se = FALSE)

fit1 <- nls(output1$ICombH ~ SSasymp(output1$tau, Asym, R0, lrc), data = output1)
summary(fit1)

#Curve Fitting - ResRat
plot_ly(output1, x= ~tau, y = ~IResRat, type = "scatter")

ggplot(output1, aes(x=tau, y=IResRat)) %>%  + # tau is from the output1 dataframe
  geom_point() + 
  stat_smooth(method="nls", formula = y ~ SSasymp(x,Asym, R0, lrc), se = FALSE)

fit2 <- nls(tempoutput1$IResRat ~ SSasymp(tempoutput1$tau, Asym, R0, lrc), data = output1)
summary(fit2)
formula(fit2)

nls(y~a*exp(b*x),start=list(a=13,b=0.1)) 
tempoutput1$ICombH



#Test
fittest <- nls(tempoutput1$ICombH ~ yf + (y0-yf) * exp(-alpha*tempoutput1$tau), data = tempoutput1, start = list(y0 = tempoutput1$ICombH[[1]],
                                                                                                                 yf = tempoutput1$ICombH[[51]]))  


fittest <- nls()
summary(fittest)

ggplot(tempoutput1, aes(x=tau, y=ICombH)) %>%  + # tau is from the output1 dataframe
  geom_point() + 
  stat_smooth(method="nls", formula = y ~ SSasymp(x,Asym, R0, lrc), se = FALSE)

#### YDiff and Alpha Curve Fitting ####

#Test for Bruteforce For Loop - tryCatch() Function
for (i in 1:10) {
  errordat <- data.frame(matrix(NA, nrow = 1, ncol=1))
  temp2 <- data.frame(matrix(NA, nrow = 1, ncol=2))
  errorparms <- data.frame(matrix(NA, nrow = 1, ncol = 11))
  colnames(errorparms)[1:11] <- c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                  "phi", "tau", "theta")
  tryCatch(expr = {
    if (i==7) stop("This is an Error")
    temp2[1,1] <- i
    temp2[1,2] <- i
  },
  error = function(e){
    message("**ERR at ", Sys.time(), "**")
    test <- print(e)
    temp2[1,1] <<- "Error"
    temp2[1,2] <<- "Error"
    errordat <<- e[[1]]
    errorparms <<- parms[1,]
    return(list(temp2, errorparms))
  })
  errordat1 <- rbind.data.frame(errordat1, errordat)
  output1 <- rbind.data.frame(output1, temp2)
  errorparms1 <- rbind.data.frame(errorparms1, errorparms)
}

#YDiff For Loop

start_time <- Sys.time()

#Parameters for the analysis - Tau removed
parms = fast_parameters(minimum = c(5.2^-1,0.6^-1,24^-1,2883.5^-1,0,0,0,0,0,0), maximum = c(520^-1,60^-1,2400^-1,288350^-1,1,0.00001,0.00001,0.0001,1,5), 
                        factor=10, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"))
parmtau <- seq(0,0.5,by=0.01)

#Model Initial Conditions and Specifications
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

for (i in 1:nrow(parms)) {
  tempoutput1 <- data.frame()
  newtemp <- data.frame()
  for (j in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
    parms2 = c(ra = parms$ra[i], rh =  parms$rh[i], ua = parms$ua[i], uh = parms$uh[i], betaAA = parms$betaAA[i], 
               betaAH = parms$betaAH[i], betaHH = parms$betaHH[i], betaHA = parms$betaHA[i], phi = parms$phi[i], 
               tau = parmtau[j], theta = parms$theta[i])
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[j]
    temp[1,2] <- rounding(out[nrow(out),5]) 
    temp[1,3] <- rounding(out[nrow(out),6]) 
    temp[1,4] <- rounding(out[nrow(out),7])
    temp[1,5] <- temp[1,3] + temp[1,4]
    temp[1,6] <- temp[1,4]/temp[1,5]
    tempoutput1 <- rbind.data.frame(tempoutput1, temp)  
  }
  colnames(tempoutput1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
  newtemp[1,1] <- tempoutput1$ICombH[tempoutput1$tau == 0]
  newtemp[1,2] <- tempoutput1$ICombH[tempoutput1$tau == 0.5]
  newtemp[1,3] <- tempoutput1$IResRat[tempoutput1$tau == 0]
  newtemp[1,4] <- tempoutput1$IResRat[tempoutput1$tau == 0.5]
  newtemp[1,5] <- newtemp[1,1] - newtemp[1,2]
  output1 <- rbind.data.frame(output1, newtemp)
  print(i)
  print(newtemp[1,3])
  print()
}

colnames(output1)[1:5] <- c("ICombH0", "ICombH05","IResRatH0", "IResRat05", "YDiff")

run0 <- output1  

end_time <- Sys.time()
time1 <- end_time - start_time

#Alpha For Loop - ICombH

start_time <- Sys.time()

parmtau <- seq(0,0.5,by=0.01)
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 200000, by = 100)

output1 <- data.frame()
errordat1 <- data.frame(matrix(NA, nrow = 0, ncol=1))
colnames(errordat1)[1] <- "data"
errorparms1 <- data.frame(matrix(NA, nrow = 0, ncol=10))
colnames(errorparms1)[1:10] <- c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                 "phi", "theta")
for (i in 1:nrow(parms)) {
  tempoutput1 <- data.frame() #Data.frame for IHComb and Tau values 
  newtemp <- data.frame(matrix(NA, nrow = 1, ncol=2)) #Alpha and Ydiff 
  errordat <- data.frame(matrix(NA, nrow = 1, ncol=1)) # Error Messages
  colnames(errordat)[1] <- "data"
  errorparms <- data.frame(matrix(NA, nrow = 1, ncol = 10)) #Faulty Parameters
  colnames(errorparms)[1:10] <- c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                  "phi","theta")
  for (j in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
    parms2 = c(ra = parms$ra[i], rh =  parms$rh[i], ua = parms$ua[i], uh = parms$uh[i], betaAA = parms$betaAA[i], 
               betaAH = parms$betaAH[i], betaHH = parms$betaHH[i], betaHA = parms$betaHA[i], phi = parms$phi[i], 
               tau = parmtau[j], theta = parms$theta[i])
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[j]
    temp[1,2] <- rounding(out[nrow(out),5]) 
    temp[1,3] <- rounding(out[nrow(out),6]) 
    temp[1,4] <- rounding(out[nrow(out),7])
    temp[1,5] <- temp[1,3] + temp[1,4]
    temp[1,6] <- temp[1,4]/temp[1,5]
    tempoutput1 <- rbind.data.frame(tempoutput1, temp)  
  }
  colnames(tempoutput1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
  tryCatch(expr = {
    fit2 <- nls(tempoutput1$ICombH ~ SSasymp(tempoutput1$tau, Asym, R0, lrc), data = tempoutput1)
    newtemp[1,1] <- (summary(fit2)$parameters[[2]]) - (summary(fit2)$parameters[[1]])
    newtemp[1,2] <- exp(summary(fit2)$parameters[[3]])
    print(newtemp)
  },
  error = function(e){
    message("**ERR at ", Sys.time(), "**")
    newtemp[1,1] <<- "Error"
    newtemp[1,2] <<- "Error"
    errordat <<- e[[1]]
    errorparms[1,] <<- parms[i,]
    return(list(newtemp, errorparms, errordat, errorparms))
  })
  errordat1 <- rbind.data.frame(errordat1, errordat)
  colnames(errordat1)[1] <- "data"
  output1 <- rbind.data.frame(output1, newtemp)
  errorparms1 <- rbind.data.frame(errorparms1, errorparms)
}

run1 <- output1

colnames(run1)[1:2] <- c("Y0-Yf", "AlphaDecay")

run1new <- run1
run1new[run1new == "Error"] <- 0 

squid <- readLines("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/squidward.txt")
print(squid)

end_time <- Sys.time()
time2 <- end_time - start_time

#Alpha For Loop - For Res Ratio

start_time <- Sys.time()

parmtau <- seq(0,0.5,by=0.01)
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 200000, by = 100)

output1 <- data.frame()
errordat1 <- data.frame(matrix(NA, nrow = 0, ncol=1))
colnames(errordat1)[1] <- "data"
errorparms1 <- data.frame(matrix(NA, nrow = 0, ncol=10))
colnames(errorparms1)[1:10] <- c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                 "phi", "theta")

for (i in 1:nrow(parms)) {
  tempoutput1 <- data.frame() #Data.frame for IHComb and Tau values 
  newtemp <- data.frame(matrix(NA, nrow = 1, ncol=2)) #Alpha and Ydiff 
  errordat <- data.frame(matrix(NA, nrow = 1, ncol=1)) # Error Messages
  colnames(errordat)[1] <- "data"
  errorparms <- data.frame(matrix(NA, nrow = 1, ncol = 10)) #Faulty Parameters
  colnames(errorparms)[1:10] <- c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                  "phi","theta")
  for (j in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
    parms2 = c(ra = parms$ra[i], rh =  parms$rh[i], ua = parms$ua[i], uh = parms$uh[i], betaAA = parms$betaAA[i], 
               betaAH = parms$betaAH[i], betaHH = parms$betaHH[i], betaHA = parms$betaHA[i], phi = parms$phi[i], 
               tau = parmtau[j], theta = parms$theta[i])
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[j]
    temp[1,2] <- rounding(out[nrow(out),5]) 
    temp[1,3] <- rounding(out[nrow(out),6]) 
    temp[1,4] <- rounding(out[nrow(out),7])
    temp[1,5] <- temp[1,3] + temp[1,4]
    temp[1,6] <- temp[1,4]/temp[1,5]
    tempoutput1 <- rbind.data.frame(tempoutput1, temp)  
  }
  colnames(tempoutput1)[1:6] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")
  tryCatch(expr = {
    fit2 <- nls(tempoutput1$IResRat ~ SSasymp(tempoutput1$tau, Asym, R0, lrc), data = tempoutput1)
    newtemp[1,1] <- (summary(fit2)$parameters[[1]])
    newtemp[1,2] <- exp(summary(fit2)$parameters[[3]])
    print(newtemp)
  },
  error = function(e){
    message("**ERR at ", Sys.time(), "**")
    newtemp[1,1] <<- "Error"
    newtemp[1,2] <<- "Error"
    errordat <<- e[[1]]
    errorparms[1,] <<- parms[i,]
    return(list(newtemp, errorparms, errordat, errorparms))
  })
  errordat1 <- rbind.data.frame(errordat1, errordat)
  colnames(errordat1)[1] <- "data"
  output1 <- rbind.data.frame(output1, newtemp)
  errorparms1 <- rbind.data.frame(errorparms1, errorparms)
}

run2 <- output1
colnames(run2new)[1:2] <- c("Yf", "AlphaDecay")

run2new <- run2
run2new[run2new == "Error"] <- 0 

end_time <- Sys.time()
time2 <- end_time - start_time

#### Sensitivity Analysis for YDiff and Alpha Parameters ####

sensit7 <- as.numeric(run1$AlphaDecay[run1$AlphaDecay != "Error"])
sensit8 <- as.numeric(run1$YDiff[run1$YDiff != "Error"])

sens<-sensitivity(x=sensit7, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                "phi","theta"))
sens1<-sensitivity(x=sensit8, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi","theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens)
df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi", "theta"), value=sens1)

p1 <- ggplot(df.equilibrium, aes(parameter, value))
p1 + geom_bar(stat="identity", fill="grey23")

p2 <- ggplot(df.equilibrium1, aes(parameter, value))
p2 + geom_bar(stat="identity", fill="grey23")

#Analysis for ResRat

colnames(run2)[1:2] <- c("YDiff", "AlphaDecay")
sensit9 <- as.numeric(run2$AlphaDecay[run2$AlphaDecay != "Error"])
sensit10 <- as.numeric(run2$YDiff[run2$YDiff != "Error"])

sens2<-sensitivity(x=sensit9, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi","theta"))
sens3<-sensitivity(x=sensit10, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                  "phi","theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens2)
df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi", "theta"), value=sens3)

p1 <- ggplot(df.equilibrium, aes(parameter, value))
p1 + geom_bar(stat="identity", fill="grey23")
p2 <- ggplot(df.equilibrium1, aes(parameter, value))
p2 + geom_bar(stat="identity", fill="grey23")

#### Sensitivity Analysis for YDiff and Alpha Parameters - Replace Errors with 0 ####

#Analysis for ICombH
run1new <- run1
run1new[run1new == "Error"] <- 0 

colnames(run1new)[colnames(run1new)=="Y0-Yf"] <- "YDiff"

sensit7 <- as.numeric(run1new$AlphaDecay)
sensit8 <- as.numeric(run1new$YDiff)

sens<-sensitivity(x=sensit7, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                "phi","theta"))
sens1<-sensitivity(x=sensit8, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi","theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens)
df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi", "theta"), value=sens1)

p1 <- ggplot(df.equilibrium, aes(parameter, value))
p1 + geom_bar(stat="identity", fill="grey23")

p2 <- ggplot(df.equilibrium1, aes(parameter, value))
p2 + geom_bar(stat="identity", fill="grey23")

#Analysis for ResRat
run2new[run2new == "Error"] <- 0 

sensit9 <- as.numeric(run2new$AlphaDecay)
sensit10 <- as.numeric(run2new$Yf)

sens2<-sensitivity(x=sensit9, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi","theta"))
sens3<-sensitivity(x=sensit10, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                  "phi","theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens2)
df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                              "phi", "theta"), value=sens3)

p1 <- ggplot(df.equilibrium, aes(parameter, value))
p1 + geom_bar(stat="identity", fill="grey23")
p2 <- ggplot(df.equilibrium1, aes(parameter, value))
p2 + geom_bar(stat="identity", fill="grey23")

#### Testbed Parameter Space Testing ####

#TEST
#Ranges for Parameter Testing
taurange <- seq(0,0.2, by=0.005)
rarange <- seq(0,0.1, by=0.001)

thetarange <- seq(0,1, by=0.01)
betaHArange <- seq(0,0.0001, by=0.00001)
betaAArange <- seq(0,0.15, by=0.005)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(taurange, betaAArange)
colnames(combparm1)[1:2] <- c("tau","betaAA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 100000, by = 500)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 25^-1, rh =  (6^-1), ua = 340^-1, uh = 28835^-1, betaAA = combparm1[i,2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.00001, phi = 0.05, tau = combparm1[i,1], theta = 0.5)
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
plot(surfaceoutput1$InfSensHum, surfaceoutput1$InfResHum, xlim = c(-0.001,1), ylim = c(-0.001,1))

surfaceoutputplot <- surfaceoutput1[,c(1,2,7)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaAA", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

p6 <- plot_ly(x = taurange, y = betaAArange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.0000001, end = 4e-05, size = 0.000002),
              colorbar = list(title = "I<sub>CombH</sub>", 
                              exponentformat= "E", thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Tau (Antibiotic Usage)", autorange = "reversed"),
         yaxis = list(title = "BetaAA", exponentformat= "E"))
p6

surfaceoutputplot <- surfaceoutput1[,c(1,2,6)]
surfaceoutputplot$new <- NULL
mat <- spread(surfaceoutputplot, key = "betaAA", value = "PropRes")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)

p7 <- plot_ly(x = taurange, y = betaAArange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.0000001, end = 1, size = 0.1),
              colorbar = list(title = "PropRes", thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "Tau (Antibiotic Usage)", autorange = "reversed"),
         yaxis = list(title = "BetaAA", exponentformat= "E"))
p7


# put next line after transpose: contours = list(start = -0.0001, end = 1, size = 0.05),

#### Tau and Phi Relationship Exploration - ICOMB ####

#TEST
#Ranges for Parameter Testing
taurange <- seq(0,0.5, by=0.01)
phirange <- seq(0,0.5, by=0.01)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(taurange, phirange)
colnames(combparm1)[1:2] <- c("tau","phi")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- c(0,9999,10000)
times <- seq(0,100000,by= 100)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = 0.1, betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = combparm1[i,2], tau = combparm1[i,1], theta = 0.5)
  out <- ode(y = init, func = amr, times = times, parms = parms1)
  print(out[nrow(out),7])
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- rounding(out[nrow(out),5]) 
  temp[1,4] <- rounding(out[nrow(out),6])
  temp[1,5] <- rounding(out[nrow(out),7])
  temp[1,6] <- temp[1,4] + temp[1,5]
  print(temp[1,6])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

colnames(surfaceoutput1)[1:6] <- c("tau","phi","SuscHum","InfSensHum", "InfResHum", "Comb Inf Res/Sens Humans")
plot(surfaceoutput1$InfSensHum, surfaceoutput1$InfResHum, xlim = c(-0.001,1), ylim = c(-0.001,1))

surfaceoutputplot <- surfaceoutput1[,c(1,2,6)]


#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "phi", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

p9 <- plot_ly(z = mat1, x = phirange, y = taurange) %>% add_surface(
  cmin = 0, cmax = 5e-05,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*", exponentformat= "E")
) %>% layout(
  title = "Equilibrium Prevalence of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    xaxis = list(title = "phi", nticks = 8, range = c(0,1.5)),
    yaxis = list(title = "tau", nticks = 8, range = c(0,0.5)),
    zaxis = list(title = 'IComb*', nticks = 8, exponentformat= "E"),
    aspectratio=list(x=0.8,y=0.8,z=0.8)))
p9 

p8 <- plot_ly(x = taurange, y = phirange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.00000001, end = 5e-05, size = 0.000005),
              colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*", exponentformat= "E")) %>% 
  layout(title = "Comb Equilibrium Prevalence of I<sub>H</sub>*",
         xaxis = list(title = "Tau"),
         yaxis = list(title = "Phi"))
p8

#### Tau and Phi Relationship Exploration - RATIO OF RESISTANCE ####

start_time <- Sys.time()

#TEST
#Ranges for Parameter Testing
taurange <- seq(0.001,0.5, by=0.001)
phirange <- seq(0.001,0.5, by=0.001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(taurange, phirange)
colnames(combparm1)[1:2] <- c("tau","phi")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- c(0,9999,10000)
times1 <- seq(0, 100000, by = 100)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = 0.1, betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = combparm1[i,2], tau = combparm1[i,1], theta = 0.8)
  out <- ode(y = init, func = amr, times = times1, parms = parms1)
  print(out[nrow(out),7])
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- out[nrow(out),5]
  temp[1,4] <- out[nrow(out),6]
  if(temp[1,4] < 1e-10) {temp[1,4] <- 0}
  temp[1,5] <- out[nrow(out),7]
  if(temp[1,5] < 1e-10) {temp[1,5] <- 0}
  temp[1,6] <- temp[1,4] + temp[1,5]
  temp[1,7] <- temp[1,5]/temp[1,6]
  print(temp[1,5])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

#temp[1,4] <- ifelse(temp[1,3] >= 0.999999999 | temp[1,4] <= 0.000000001, 0, temp[1,4])
#temp[1,5] <- ifelse(temp[1,3] >= 0.999999999 | temp[1,5] <= 0.000000001, 0, temp[1,5])

colnames(surfaceoutput1)[1:7] <- c("tau","phi","SuscHum","InfSensHum", "InfResHum", "IHCOMB", "ResRatio")

#surfaceoutput1$ResRatio[is.nan(surfaceoutput1$ResRatio)] <- 0

surfaceoutputplot <- surfaceoutput1[,c(1,2,7)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "phi", value = "ResRatio")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

p9 <- plot_ly(z = mat1, x = phirange, y = taurange) %>% add_surface(
  cmin = 0, cmax = 1,
  colorbar = list(title = "Resistance Ratio")) %>% layout(
    title = "Resistance Ratio",
    scene = list(
      camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
      xaxis = list(title = "phi", nticks = 8, range = c(0,0.15)),
      yaxis = list(title = "tau", nticks = 8, range = c(0,0.15)),
      zaxis = list(title = 'ResRat', nticks = 8),
      aspectratio=list(x=0.8,y=0.8,z=0.8)))
p9 

p8 <- plot_ly(x = taurange, y = phirange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.00000001, end = 1, size = 0.1),
              colorbar = list(title = "Resistance Ratio",exponentformat= "E")) %>% 
  layout(title = "IHCOMB",
         xaxis = list(title = "Tau"),
         yaxis = list(title = "Phi"))
p8

#Phi Tau Ratio

surfaceoutput2 <- surfaceoutput1
surfaceoutput2$phitau <- surfaceoutput2$tau/surfaceoutput2$phi
surfaceoutput2$logphitau <- log10(surfaceoutput2$phitau)

df.unique <- surfaceoutput2[!duplicated(surfaceoutput2$phitau), ]
df.unique

ggplot(surfaceoutput2, aes(logphitau, ResRatio)) + geom_line(size=1.5)
ggplot(surfaceoutput2, aes(logphitau, IHCOMB)) + geom_point(size=1)
plot(surfaceoutput2$logphitau, surfaceoutput2$ResRatio)

end_time <- Sys.time()

end_time - start_time

plot_ly(surfaceoutput2, x= ~logphitau, y = ~InfSensHum, type = "bar", name = "Sens Inf Humans") %>%
  add_trace(y= ~InfResHum, name = "Res Inf Humans") %>% 
  layout(yaxis = list(title = "Proportion Infected", exponentformat= "E", range = c(0,1E-4), showline = TRUE),
         xaxis = list(title = "Log Phi/Tau"),
         legend = list(orientation = "v", x = 1.0, y=0.5), showlegend = T,
         barmode = "stack") 



#         annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
#                            xshift =3))
test <- surfaceoutput2

test2 <- data.frame(tapply(test$InfSensHum, cut(test$logphitau, seq(-2.3, 2.3, by=0.01)), mean))

#test$bin <- bin(surfaceoutput2$logphitau, nbins = 200, labels = NULL, method = "content", na.omit= TRUE)

#### Evaluating the Ratio of Foodborne Infection - Before and After the Intervention ####

betaHArange <- seq(0, 0.00005, by = 0.000001)
betaAArange <- seq(0, 1, by = 0.1)

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times1 <- seq(0, 10000, by = 100)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:length(betaAArange)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=12))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = betaAArange[i], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = 0.1, tau = 0, theta = 0.5)
  parms2 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = betaAArange[i], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = 0.1, tau = 0.05, theta = 0.5)
  out <- ode(y = init, func = amr, times = times1, parms = parms1)
  out1 <- ode(y = init, func = amr, times = times1, parms = parms2)
  temp[1,1] <- betaAArange[i]
  temp[1,2] <- 0
  temp[1,3] <- out[nrow(out),6]
  if(temp[1,3] < 1e-10) {temp[1,3] <- 0}
  temp[1,4] <- out[nrow(out),7]
  if(temp[1,4] < 1e-10) {temp[1,4] <- 0}
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  
  temp[1,7] <- 0.05
  temp[1,8] <- out1[nrow(out1),6]
  if(temp[1,8] < 1e-10) {temp[1,8] <- 0}
  temp[1,9] <- out1[nrow(out1),7]
  if(temp[1,9] < 1e-10) {temp[1,9] <- 0}
  temp[1,10] <- temp[1,8] + temp[1,9]
  temp[1,11] <- temp[1,9]/temp[1,10]
  
  temp[1,12] <- temp[1,10]/temp[1,5]
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

colnames(surfaceoutput1)[1:12] <- c("betaAA","tau","InfSensHum0", "InfResHum0", "IHCOMB0", "ResRatio0",
                                    "tau","InfSensHum005", "InfResHum005", "IHCOMB005", "ResRatio005",
                                    "ICOMBRat0005")

plot_ly(surfaceoutput1, x= ~betaAA, y = ~ICOMBRat0005, type = "scatter")

###-----------------------------------------------

start_time <- Sys.time()

#Ranges for Parameter Testing
betaAArange<- seq(0,0.5, by=0.01)
betaHArange <- seq(0,0.00005, by=0.000001)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(betaAArange, betaHArange)
colnames(combparm1)[1:2] <- c("betaAA","betaHA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times1 <- seq(0, 10000, by = 10)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=13))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = combparm1[i,2], phi = 0.1, tau = 0, theta = 0.5)
  parms2 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = combparm1[i,2], phi = 0.1, tau = 0.1, theta = 0.5)
  out <- ode(y = init, func = amr, times = times1, parms = parms1)
  out1 <- ode(y = init, func = amr, times = times1, parms = parms2)
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- 0
  temp[1,4] <- out[nrow(out),6]
  if(temp[1,4] < 1e-10) {temp[1,4] <- 0}
  temp[1,5] <- out[nrow(out),7]
  if(temp[1,5] < 1e-10) {temp[1,5] <- 0}
  temp[1,6] <- temp[1,4] + temp[1,5]
  temp[1,7] <- temp[1,5]/temp[1,6]
  
  temp[1,8] <- 0.1
  temp[1,9] <- out1[nrow(out1),6]
  if(temp[1,9] < 1e-10) {temp[1,9] <- 0}
  temp[1,10] <- out1[nrow(out1),7]
  if(temp[1,10] < 1e-10) {temp[1,10] <- 0}
  temp[1,11] <- temp[1,9] + temp[1,10]
  temp[1,12] <- temp[1,10]/temp[1,11]
  
  temp[1,13] <- 1-(temp[1,11]/temp[1,6])
  temp[1,14] <- temp[1,12]
  print(temp[1,13])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

colnames(surfaceoutput1)[1:13] <- c("betaAA","betaHA","tau1","InfSensHum0", "InfResHum0", "IHCOMB0", "ResRatio0",
                                    "tau2","InfSensHum005", "InfResHum005", "IHCOMB005", "ResRatio005",
                                    "ICOMBRat0005")
surfaceanal <- surfaceoutput1[,c(1,2,6,11,13)]
surfaceoutputplot <- surfaceoutput1[,c(1,2,13)]

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaHA", value = "ICOMBRat0005")
row.names(mat) <- mat$betaAA
mat$betaAA<- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

plot_ly(z = mat1, x = betaHArange, y = betaAArange) %>% add_surface(
  cmin = 0, cmax = 1,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*")
) %>% layout(
  title = "Equilibrium Prevalence Ratio of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    xaxis = list(title = "BetaHA", nticks = 8, range = c(0,0.00005), exponentformat= "E"),
    yaxis = list(title = "BetaAA", nticks = 8, range = c(0,0.5)),
    zaxis = list(title = '1-(W/W-out)', nticks = 8),
    aspectratio=list(x=0.8,y=0.8,z=0.8)))

end_time <- Sys.time()

end_time - start_time
###----------------------BETAAH-------------------------

#Ranges for Parameter Testing
betaAArange<- seq(0,0.5, by=0.01)
betaAHrange <- seq(0,0.000005, by=0.0000001)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(betaAArange, betaAHrange)
colnames(combparm1)[1:2] <- c("betaAA","betaAH")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times1 <- seq(0, 10000, by = 100)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations

for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=13))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = combparm1[i,2], betaHH = 0.000001, 
             betaHA = 0.00001, phi = 0.1, tau = 0, theta = 0.5)
  parms2 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = combparm1[i,2], betaHH = 0.000001, 
             betaHA = 0.00001, phi = 0.1, tau = 0.1, theta = 0.5)
  out <- ode(y = init, func = amr, times = times1, parms = parms1)
  out1 <- ode(y = init, func = amr, times = times1, parms = parms2)
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- 0
  temp[1,4] <- out[nrow(out),6]
  if(temp[1,4] < 1e-10) {temp[1,4] <- 0}
  temp[1,5] <- out[nrow(out),7]
  if(temp[1,5] < 1e-10) {temp[1,5] <- 0}
  temp[1,6] <- temp[1,4] + temp[1,5]
  temp[1,7] <- temp[1,5]/temp[1,6]
  
  temp[1,8] <- 0.1
  temp[1,9] <- out1[nrow(out1),6]
  if(temp[1,9] < 1e-10) {temp[1,9] <- 0}
  temp[1,10] <- out1[nrow(out1),7]
  if(temp[1,10] < 1e-10) {temp[1,10] <- 0}
  temp[1,11] <- temp[1,9] + temp[1,10]
  temp[1,12] <- temp[1,10]/temp[1,11]
  
  temp[1,13] <- 1-(temp[1,11]/temp[1,6])
  temp[1,14] <- temp[1,12]
  print(temp[1,13])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

colnames(surfaceoutput1)[1:13] <- c("betaAA","betaAH","tau1","InfSensHum0", "InfResHum0", "IHCOMB0", "ResRatio0",
                                    "tau2","InfSensHum005", "InfResHum005", "IHCOMB005", "ResRatio005",
                                    "ICOMBRat0005")
surfaceanal <- surfaceoutput1[,c(1,2,6,11,13)]
surfaceoutputplot <- surfaceoutput1[,c(1,2,13)]

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaAH", value = "ICOMBRat0005")
row.names(mat) <- mat$betaAA
mat$betaAA<- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

plot_ly(z = mat1, x = betaAHrange, y = betaAArange) %>% add_surface(
  cmin = 0, cmax = 1,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*")
) %>% layout(
  title = "Equilibrium Prevalence Ratio of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    xaxis = list(title = "BetaAH", nticks = 8, range = c(0,0.000005), exponentformat= "E"),
    yaxis = list(title = "BetaAA", nticks = 8, range = c(0,0.5)),
    zaxis = list(title = '1-(W/W-out)', nticks = 8),
    aspectratio=list(x=0.8,y=0.8,z=0.8)))

###----------------------Phi-------------------------

#Ranges for Parameter Testing
betaAArange<- seq(0,0.5, by=0.01)
phirange <- seq(0,1, by=0.1)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(betaAArange, phirange)
colnames(combparm1)[1:2] <- c("betaAA","phi")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times1 <- seq(0, 10000, by = 100)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations

for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=13))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = combparm1[i,2], tau = 0, theta = 0.5)
  parms2 = c(ra = 52^-1, rh =  6^-1, ua = 28835^-1, uh = 240^-1, betaAA = combparm1[i,1], betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = 0.00001, phi = combparm1[i,2], tau = 0.1, theta = 0.5)
  out <- ode(y = init, func = amr, times = times1, parms = parms1)
  out1 <- ode(y = init, func = amr, times = times1, parms = parms2)
  temp[1,1] <- combparm1[i,1]
  temp[1,2] <- combparm1[i,2]
  temp[1,3] <- 0
  temp[1,4] <- out[nrow(out),6]
  if(temp[1,4] < 1e-10) {temp[1,4] <- 0}
  temp[1,5] <- out[nrow(out),7]
  if(temp[1,5] < 1e-10) {temp[1,5] <- 0}
  temp[1,6] <- temp[1,4] + temp[1,5]
  temp[1,7] <- temp[1,5]/temp[1,6]
  
  temp[1,8] <- 0.1
  temp[1,9] <- out1[nrow(out1),6]
  if(temp[1,9] < 1e-10) {temp[1,9] <- 0}
  temp[1,10] <- out1[nrow(out1),7]
  if(temp[1,10] < 1e-10) {temp[1,10] <- 0}
  temp[1,11] <- temp[1,9] + temp[1,10]
  temp[1,12] <- temp[1,10]/temp[1,11]
  
  temp[1,13] <- 1-(temp[1,11]/temp[1,6])
  temp[1,14] <- temp[1,12]
  print(temp[1,13])
  surfaceoutput1 <- rbind.data.frame(surfaceoutput1, temp)
}

colnames(surfaceoutput1)[1:13] <- c("betaAA","phi","tau1","InfSensHum0", "InfResHum0", "IHCOMB0", "ResRatio0",
                                    "tau2","InfSensHum005", "InfResHum005", "IHCOMB005", "ResRatio005",
                                    "ICOMBRat0005")
surfaceanal <- surfaceoutput1[,c(1,2,6,11,13)]
surfaceoutputplot <- surfaceoutput1[,c(1,2,13)]

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "phi", value = "ICOMBRat0005")
row.names(mat) <- mat$betaAA
mat$betaAA<- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

plot_ly(z = mat1, x = phirange, y = betaAArange) %>% add_surface(
  cmin = 0, cmax = 1,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*")
) %>% layout(
  title = "Equilibrium Prevalence Ratio of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    xaxis = list(title = "phi", nticks = 8, range = c(0,1)),
    yaxis = list(title = "BetaAA", nticks = 8, range = c(0,0.5)),
    zaxis = list(title = '1-(W/W-out)', nticks = 8),
    aspectratio=list(x=0.8,y=0.8,z=0.8)))

#### Testbed Parameter Space Testing ####

#TEST
#Ranges for Parameter Testing

betaHArange <- seq(0.000005,0.00001, by=0.0000001)
betaAArange <- seq(0.045,0.09, by=0.001)
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
  parms1 = c(ra = 25^-1, rh =  (6^-1), ua = 240^-1, uh = 28835^-1, betaAA = combparm1[i,2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = combparm1[i,1], phi = 0.05, tau = 0, theta = 0.5)
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

p6 <- plot_ly(x = betaHArange, y = betaAArange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.0000001, end = 4e-05, size = 0.000002),
              colorbar = list(title = "I<sub>CombH</sub>", 
                              exponentformat= "E", thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "BetaHA", exponentformat= "E", autorange = "reversed"),
         yaxis = list(title = "BetaAA"))
p6

mat2 <- mat1
mat2[mat2 < 100] <- 1.388e-05


p9 <- plot_ly(z = mat1, x = betaAArange, y = betaHArange) %>% add_surface()

p9 <- plot_ly(z = mat1, x = betaAArange, y = betaHArange) %>% add_surface(
  cmin = 0, cmax = 3e-05,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*", exponentformat= "E")
) %>% layout(
  title = "Equilibrium Prevalence of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    yaxis = list(title = "betaHA", nticks = 8, range = c(0.000005,0.00001), exponentformat= "E"),
    xaxis = list(title = "betaAA", nticks = 8, range = c(0.045,0.09)),
    zaxis = list(title = 'IComb*', nticks = 8, exponentformat= "E"),
    aspectratio=list(x=0.8,y=0.8,z=0.8))) %>% 
  add_surface(z = mat2, opacity = 0.8, color = "black")

p9 

mat3 <- mat1*100000
mat4 <- mat2*100000


plot_ly(z = mat3, x = betaaaperc, y = betahaperc) %>% add_surface(
  cmin = 0, cmax = 3,
  colorbar = list(title = "ICombH")
) %>% layout(
  title = "Equilibrium Prevalence of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    yaxis = list(nticks = 8, range = c(50,100)),
    xaxis = list( nticks = 8, range = c(50,100)),
    zaxis = list(nticks = 8),
    aspectratio=list(x=0.8,y=0.8,z=0.8))) %>% 
  add_surface(z = mat4, opacity = 0.8, color = "black")

plot_ly(z = mat3, x = betaaaperc, y = betahaperc) %>% add_surface(
  cmin = 0, cmax = 3,
  colorbar = list(title = "ICombH")
) %>% layout(
  title = "Equilibrium Prevalence of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    yaxis = list(nticks = 8, range = c(50,100)),
    xaxis = list(nticks = 8, range = c(50,100)),
    zaxis = list(nticks = 8, range = c(1.388,3.1)),
    aspectratio=list(x=0.8,y=0.8,z=0.8))) 

p15 <- plot_ly(x = betaaaperc, y = betahaperc, z = mat3, type = "contour",
               contours = list(start = 1.388, end = 3, size = 0.2, showlabels = TRUE,
                               labelfont = list(size = 12, color = 'white')),
               colorbar = list(title = "I<sub>CombH</sub>", 
                               thickness = 35, len = 1)) %>% 
  layout(xaxis = list(title = "BetaAA (%)"),
         yaxis = list(title = "BetaHA (%)"))
p15
