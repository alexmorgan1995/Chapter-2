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

init <- c(Sa=0.99, Ia=0.01, Ira=0.01, Sh=1, Ih=0, Irh=0)
times1 <- seq(0,1000000,by=1000)

#Need to Specify Model Parameters
parms = c(ra = 52^-1, rh =  6^-1, uh = 28835^-1, ua = 240^-1, betaAA = 0.10, betaAH = 0.000001, betaHH = 0.000001, 
          betaHA = 0.00001, phi = 0.1, tau = 0.1, theta = 0.5)

out <- ode(y = init, func = amr, times = times1, parms = parms)

head(out,100) #Provides the first 10 timesteps for the model output

par(mfrow=c(1,2))
plot(out[,"time"],out[,"Sa"], xlab = "Time (Days)",ylab="Proportion (Animals)", type="l", lwd=2, col="green", ylim = c(-0.01, 1))
lines(out[,"time"],out[,"Ia"], type="l", lwd=5, col="red")
lines(out[,"time"],out[,"Ira"],type="l", lwd=2, col="blue")
legend(x=450, y=0.97, legend= c("Susceptible", "Infected (S)", "Infected (R)"),col=c("green","red","blue"),lty=1,cex=0.9)

plot(out[,"time"],out[,"Sh"], xlab = "Time (Days)",ylab="Proportion(Humans)", type="l", lwd=2, col="green", ylim = c(0, 0.00002))
lines(out[,"time"],out[,"Ih"], type="l", lwd=5, col="red")
lines(out[,"time"],out[,"Irh"], type="l", lwd=2, col="blue")
legend(x=450, y=1.95e-05, legend= c("Susceptible", "Infected (S)", "Infected (R)"),col=c("green","red","blue"),lty=1,cex=0.9)

Icomb <- as.numeric(out[nrow(out),6]) + as.numeric(out[nrow(out),7])  
print(Icomb)
#plot(out[,"time"], out[,"Ia"]+ out[,"Ira"], type="l", lwd=2, col="black", ylab="IComb*", xlab="Time")

#### Function for Parameter Combinations - From FAST ####
start_time <- Sys.time()

parms = fast_parameters(minimum = c(5.2^-1,0.6^-1,24^-1,2883.5^-1,0,0,0,0,0,0,0), maximum = c(520^-1,60^-1,240^-1,288350^-1,1,0.00001,0.00001,0.0001,1,1,5), 
                        factor=11, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "tau", "theta"))

#### Creating Model Output for Parameter Combinations from FAST Parameters #### 

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0,100000, by = 1000) 

output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
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
  print(temp[1,3])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:5] <- c("SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat")

test <- output[complete.cases(output),] # to get rid of NA found in the dataframe - make it analysable int he sensitivity analysis

#### Sensitivity Analysis ####

sensit <- test$ICombH #Creating Variable for the output variable of interest
sens<-sensitivity(x=sensit, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                               "phi", "tau", "theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi","tau", "theta"), value=sens)

p <- ggplot(df.equilibrium, aes(parameter, value))
p + geom_bar(stat="identity", fill="grey23")

end_time <- Sys.time()

end_time - start_time

#### Comb IRH and IH Measure Testing - Effect of Treatment ####

#Testing for the Effect of Changing Treatment on the Combined Measure

#parmtau <- seq(0,0.5,by=0.03)
parmtau <- seq(0,0.5,by=0.01)

init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

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
           list(x = 0.45, y = 1.82e-05, text = "Baseline Level of FB Disease", yanchor = "bottom",
                showarrow = FALSE, xshift= -55, yshift = 3, font = list(size = 15)),
           list(x = 0.45, y = output1$ICombH[11], text = "New Level of FB Disease", yanchor = "bottom",
                showarrow = FALSE, xshift= -55, yshift = 3, font = list(size = 15, color = "red")),
           list(ax = 0.1, ay = 1.9e-05, x = 0.1, y = output1$ICombH[11],
                axref = "x", ayref = "y", xref = "x", yref = "y", showarrow= TRUE, arrowhead=1, 
                arrowsize = 1.2, arrowwidth=2.5, arrowcolor= "red"))) %>%
  add_segments(x=-0.01, xend = 0.500001, y=1.82e-05, yend = 1.82e-05, line = list(color = "black", dash = "dot", width = 2),
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

fit2 <- nls(output1$IResRat ~ SSasymp(output1$tau, Asym, R0, lrc), data = output1)
summary(fit2)
formula(fit1)

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
parms = fast_parameters(minimum = c(5.2^-1,0.6^-1,24^-1,2883.5^-1,0,0,0,0,0,0), maximum = c(520^-1,60^-1,240^-1,288350^-1,1,0.00001,0.00001,0.0001,1,5), 
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
  newtemp[1,3] <- newtemp[1,1] - newtemp[1,2]
  output1 <- rbind.data.frame(output1, newtemp)
  print(i)
  print(newtemp[1,3])
}

colnames(output1)[1:3] <- c("ICombH0", "ICombH05", "YDiff")

run0 <- output1  

end_time <- Sys.time()
time1 <- end_time - start_time

#Alpha For Loop

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


end_time <- Sys.time()
time2 <- end_time - start_time

#### Sensitivity Analysis for YDiff and Alpha Parameters ####

sensit1 <- run0$YDiff
sensit3 <- as.numeric(run1$AlphaDecay[run1$AlphaDecay != "Error"])
                              
sens<-sensitivity(x=sensit1, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                               "phi","theta"))
sens1<-sensitivity(x=sensit3, numberf=10, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                "phi","theta"))

df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens)
df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta"), value=sens1)

p1 <- ggplot(df.equilibrium, aes(parameter, value))
p1 + geom_bar(stat="identity", fill="grey23")

p2 <- ggplot(df.equilibrium1, aes(parameter, value))
p2 + geom_bar(stat="identity", fill="grey23")

#### Testbed Parameter Space Testing ####

#TEST
#Ranges for Parameter Testing
taurange <- seq(0,0.5, by=0.01)
thetarange <- seq(0,1, by=0.01)
betaHArange <- seq(0,0.00002, by=0.000001)

#betaHArange <- seq(0,0.0001, by=0.00001)

#Creating Possible Combinations of Parameters
combparm1 <- NULL
combparm1 <- expand.grid(taurange, betaHArange)
colnames(combparm1)[1:2] <- c("tau","betaHA")

#Setting up the initial Conditions for the Model
init <- c(Sa=0.99, Ia=0.01, Ira=0, Sh=1, Ih=0, Irh=0)
times <- seq(0, 100000, by = 1000)

#Creating Dummy Data Frame for For Loop
surfaceoutput1 <- data.frame()

#For Loop to Create Output for the Parameter Combinations
for (i in 1:nrow(combparm1)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms1 = c(ra = 52^-1, rh =  6^-1, ua = 240^-1, uh = 28835^-1, betaAA = 0.05, betaAH = 0.000001, betaHH = 0.000001, 
             betaHA = combparm1[i,2], phi = 0.1, tau = combparm1[i,1], theta = 0.5)
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

#temp[1,6] <- out[nrow(out),6] + out[nrow(out),7]
colnames(surfaceoutput1)[1:6] <- c("tau","betaHA","SuscHum","InfSensHum", "InfResHum", "Comb Inf Res/Sens Humans")
plot(surfaceoutput1$InfSensHum, surfaceoutput1$InfResHum, xlim = c(-0.001,1), ylim = c(-0.001,1))

surfaceoutputplot <- surfaceoutput1[,c(1,2,6)]

#surfaceoutputplot$new <- round(surfaceoutputplot$`Comb Inf Res/Sens Humans`, digits = 6)

#Spread the Parameter Combinations out so it can be plotted
surfaceoutputplot$new <- NULL

mat <- spread(surfaceoutputplot, key = "betaHA", value = "Comb Inf Res/Sens Humans")
row.names(mat) <- mat$tau
mat$tau <- NULL
mat1 <- data.matrix(mat)

#Plotting a vector plot and a surface plot
#cmin = 0, cmax = 0.00045,
#contours = list(start = -0.0000001, end = 0.00045, size = 0.00005)

p5 <- plot_ly(z = mat1, x = betaHArange, y = taurange) %>% add_surface(
  cmin = 0, cmax = 0.8e-04,
  colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*", exponentformat= "E")
) %>% layout(
  title = "Equilibrium Prevalence of I<sub>H</sub>*",
  scene = list(
    camera = list(eye = list(x = -1.25, y = 1.25, z = 0.5)),
    xaxis = list(title = "betaHA", nticks = 8, range = c(0,0.00002), exponentformat= "E"),
    yaxis = list(title = "tau", nticks = 8, range = c(0,0.5)),
    zaxis = list(title = 'IComb*', nticks = 8, exponentformat= "E"),
    aspectratio=list(x=0.8,y=0.8,z=0.8)))
p5 

p6 <- plot_ly(x = taurange, y = betaHArange, z = mat1, type = "contour", transpose = TRUE,
              contours = list(start = -0.0000001, end = 0.8e-04, size = 0.00001),
              colorbar = list(title = "I<sub>RH</sub>* + I<sub>H</sub>*", exponentformat= "E")) %>% 
  layout(title = "Comb Equilibrium Prevalence of I<sub>H</sub>*",
         xaxis = list(title = "Tau", autorange = "reversed"),
         yaxis = list(title = "BetaHA", exponentformat= "E"))
p6

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
taurange <- seq(0.001,0.5, by=0.005)
phirange <- seq(0.001,0.5, by=0.005)

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
             betaHA = 0.00001, phi = combparm1[i,2], tau = combparm1[i,1], theta = 0.5)
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

plot_ly(surfaceoutput2, x= ~logphitau, y = ~InfSensHum, type = "bar", name = "Sens Inf Humans") %>%
  add_trace(y= ~InfResHum, name = "Res Inf Humans") %>% 
  layout(yaxis = list(title = "Proportion Infected", exponentformat= "E", range = c(0,1E-4), showline = TRUE),
         xaxis = list(title = "Log Phi/Tau"),
         legend = list(orientation = "v", x = 1.0, y=0.5), showlegend = T,
         barmode = "stack") 

end_time <- Sys.time()

end_time - start_time

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
