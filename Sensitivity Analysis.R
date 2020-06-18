library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast")

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

# Sensitivity Analysis ----------------------------------------------------
# Joint Parameters

#We do this for Tetracycline for the Parameter Bounds
#All other parameters will be included in the ranges

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
tau_range <- c(0, 0.0106) # Comparing Baseline Average with Curtailment

#These PArameters Are Based on MAP from Model Fitting
parms = fast_parameters(minimum = c(600^-1, 55^-1, 2400^-1, 288350^-1, 
                                    0.0074716, 0.000001, 0.000001, 0.000001, 
                                    0, 0, 0), 
                        maximum = c(6^-1, 0.55^-1, 24^-1, 2883.5^-1, 
                                    0.74716, 0.0001, 0.0001, 0.0001, 
                                    0.10948457, 0.08345866, 1), 
                        factor=11, names = c("ra", "rh" ,"ua", "uh", 
                                             "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "theta", "alpha"))

# What Parameters Cause the Largest Relative Increase? --------------------

tauoutput <- data.frame(matrix(nrow = 0, ncol = 3))

for (j in 1:nrow(parms)) {
  temp <- numeric(2)
  for (i in 1:length(tau_range)) {
    parms2 = c(ra = parms$ra[j], rh = parms$rh[j], ua = parms$ua[j], uh = parms$uh[j], betaAA = parms$betaAA[j],
               betaAH = parms$betaAH[j], betaHH = parms$betaHH[j], betaHA = parms$betaHA[j], phi=parms$phi[j],
               theta=parms$theta[j], alpha = parms$alpha[j], tau = tau_range[i])
    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[i] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
  }
  tauoutput <- rbind(tauoutput, c(temp[1], temp[2], abs(temp[1] - temp[2]), parms2[parms2 != "tau"] ))
  print(j/nrow(parms))
}

colnames(tauoutput) <- c("curt", "usage", "diff") 

#Running the FAST Sensitivity Analysis
tauoutput1 <- tauoutput 
tauoutput1$inc <- ((tauoutput1$curt / tauoutput1$usage) - 1)* 100 # % Change from the current usage scenario

tauoutput1$inc[is.nan(tauoutput1$inc)] <- 0; neg <- tauoutput1[tauoutput1$inc < 0,] 
tauanalysis <- tauoutput1$inc[!is.infinite(tauoutput1$inc)]
tauanalysis <- tauanalysis[tauanalysis < quantile(tauanalysis, 0.99)]
#This step changes all NA input to 0 and removes infinities
#removes all negative changes - might need to review
#The tail of the distribution has also been trimmed to prevent massive artificial increases from showing up

#We then view the Distribution of Increases Above Baseline
hist(tauanalysis, xlab = "% Increase above Baseline (Tau = 0.0106)", breaks = 50)

# What Parameters Can Compensate? ------------------------------------------

tauoutput <- data.frame(matrix(nrow = 0, ncol = 3))

for (j in 1:nrow(parms)) {
  parms2 = c(ra = parms$ra[j], rh = parms$rh[j], ua = parms$ua[j], uh = parms$uh[j], betaAA = parms$betaAA[j],
             betaAH = parms$betaAH[j], betaHH = parms$betaHH[j], betaHA = parms$betaHA[j], phi=parms$phi[j],
             theta=parms$theta[j], alpha = parms$alpha[j], tau = 0)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
  tauoutput <- rbind(tauoutput, c(temp ,abs(temp - 3.262)))
  print(j/nrow(parms))
}

colnames(tauoutput) <- c("IComb0","diff") 

tauoutput1 <- tauoutput 

tauoutput1$inc <- ((tauoutput1$IComb0 / 3.262) - 1)* 100 # % Increase from the current usage scenario
tauoutput1$inc[is.nan(tauoutput1$inc)] <- 0; neg <- tauoutput1[tauoutput1$inc < 0,] 
tauanalysis2 <- tauoutput1$inc[!is.infinite(tauoutput1$inc)]
tauanalysis2 <- tauanalysis2[tauanalysis2 < quantile(tauanalysis2, 0.99)]
#This step changes all NA input to 0 and removes infinities
#removes all negative changes - might need to review
#The tail of the distribution has also been trimmed to prevent massive artificial increases from showing up

hist(tauanalysis2, xlab = "% Increase above Baseline (Tau = 0.032)", breaks = 50)

# Plotting Sensitivity Analysis -------------------------------------------

#Increase
sensit <- tauanalysis 
sens <- sensitivity(x=sensit, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi", "theta", "alpha"))
df.equilibrium <- NULL; df.equilibrium <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                     "phi", "theta", "alpha"), value=sens)

#Compensation
sensit1 <- tauanalysis2 #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=11, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                 "phi", "theta", "alpha"))
df.equilibrium1 <- NULL; df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                     "phi", "theta", "alpha"), value=sens1)

#Plotting
p1 <- ggplot(df.equilibrium, aes(x = reorder(parameter, -value), y = value)) + geom_bar(stat="identity", fill="lightgrey", col = "black", width  = 0.8) + theme_bw() + 
  scale_y_continuous(limits = c(0,  max(df.equilibrium$value)*1.1), expand = c(0, 0), name = "Variance") + 
  scale_x_discrete(expand = c(0, 0.7), name = "Parameter", 
                   labels = c(expression(alpha), expression(phi), expression(r[A]), expression(beta[AA]), expression(mu[A]), 
                              expression(theta), expression(beta[HA]), expression(beta[AH]), expression(r[H]), expression(beta[HH]), expression(mu[H]))) +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5"))+
  theme(legend.text=element_text(size=14), axis.text=element_text(size=14),
        axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), plot.margin = unit(c(0.15,0.4,0.15,0.55), "cm"))

p2 <- ggplot(df.equilibrium1, aes(x = reorder(parameter, -value), y = value)) + geom_bar(stat="identity", fill="lightgrey", col = "black", width  = 0.8) + theme_bw() + 
  scale_y_continuous(limits = c(0,  max(df.equilibrium1$value)*1.1), expand = c(0, 0), name = "Variance") + 
  scale_x_discrete(expand = c(0, 0.7), name = "Parameter", 
                   labels = c(expression(r[H]), expression(beta[HA]), expression(beta[AA]), expression(r[A]), expression(beta[AH]),
                              expression(alpha), expression(mu[A]), expression(phi), expression(mu[H]), expression(beta[HH]),  expression(theta))) +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5"))+
  theme(legend.text=element_text(size=14), axis.text=element_text(size=14),
        axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), plot.margin = unit(c(0.15,0.4,0.15,0.55), "cm"))

sensplot <- ggarrange(p1,p2, nrow = 2, ncol = 1,
                      align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(sensplot, filename = "Sensitivity.png", dpi = 300, type = "cairo", width = 5, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft Figures")


# Effect on Parameters ----------------------------------------------------

parmdetails <- rbind(data.frame("Parameter" = "betaAA", "Value" = seq(0, 0.74716, by = 0.74716/50)),
                     data.frame("Parameter" = "betaHA", "Value" = seq(0, 0.0001, by = 0.0001/50)),
                     data.frame("Parameter" = "betaHH", "Value" = seq(0, 0.0001, by = 0.0001/50)),
                     data.frame("Parameter" = "betaAH", "Value" = seq(0, 0.0001, by = 0.0001/50)),
                     data.frame("Parameter" = "phi", "Value" = seq(0, 0.10948457, by = 0.10948457/50)),
                     data.frame("Parameter" = "theta", "Value" = seq(0, 0.08345866, by = 0.08345866/50)),
                     data.frame("Parameter" = "alpha", "Value" = seq(0, 1, by = 1/50)),
                     data.frame("Parameter" = "rh", "Value" = seq(0, 0.55^-1, by = 0.55^-1/50)),
                     data.frame("Parameter" = "ra", "Value" = seq(0, 6^-1, by = 6^-1/50)),
                     data.frame("Parameter" = "uh", "Value" = seq(0, 2883.5^-1, by = 2883.5^-1/50)),
                     data.frame("Parameter" = "ua", "Value" = seq(0, 24^-1, by = 24^-1/50)))

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
tau_range <- c(0, 0.0106)

parms = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.074716), betaAH = 0.00001, betaHH = 0.00001, 
          betaHA = (0.00001), phi = 0.010948457, theta = 0.008345866, alpha = 0.28247322)

for (j in 1:length(unique(parmdetails[,1]))) { 
  output <- data.frame()
  for (x in 1:length(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]),2])) { #for the individual parameter values in the sequence
    temp1 <- data.frame()
    for (i in 1:length(tau_range)) {
      temp <- data.frame(matrix(nrow = 0, ncol=3))
      parmstemp <- c(parms, tau = tau_range[i])
      parmstemp[as.character(unique(parmdetails[,1])[j])] <- parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]
      out <- ode(y = init, func = amr, times = times, parms = parmstemp)
      temp[1,1] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
      temp[1,2] <- as.character(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]) #what is the parameter value used
      temp[1,3] <- as.character(unique(parmdetails[,1])[j]) # what is the parameter explored 
      temp1 <- rbind.data.frame(temp1, temp)
    }
    print(temp1)
    output <- rbind(output, stringsAsFactors = FALSE,
                    c(as.numeric(temp1[1,1]), 
                      as.numeric(temp1[2,1]), 
                      as.numeric(abs(temp1[1,1] - temp1[2,1])),
                      as.numeric(abs(((temp1[1,1] / temp1[2,1]) - 1)* 100)),
                      as.numeric(abs(((temp1[1,1] / 3.382) - 1)* 100)),
                      as.numeric(temp1[i,2]),
                      temp1[i,3]))
  }
  colnames(output)[1:7] <- c("ICombHCurt", "ICombHUsage","IDiff", "PercInc", "RelInc", "ParmValue","Parm")
  output$Parm <- as.factor(output$Parm)
  assign(paste("output", as.character(unique(parmdetails[,1])[j]), sep=""), output) 
}

#Absolute Diff
pbetaAA <- ggplot(outputbetaAA, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="BetaAA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaHH <- ggplot(outputbetaHH, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="BetaHH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaHA <- ggplot(outputbetaHA, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="BetaHA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaAH <- ggplot(outputbetaAH, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="BetaAH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

palpha <- ggplot(outputalpha, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Alpha Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pphi <- ggplot(outputphi, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Phi Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
ptheta <- ggplot(outputtheta, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Theta Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

pra <- ggplot(outputra, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,500)) +
  labs(x ="rA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
prh <- ggplot(outputrh, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,500)) +
  labs(x ="rH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pua <- ggplot(outputua, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,500)) +
  labs(x ="uA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
puh <- ggplot(outputuh, aes(x = as.numeric(ParmValue), y = as.numeric(PercInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,500)) +
  labs(x ="uH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

plot_grid(plot_grid(pbetaAA, pbetaHH, pbetaHA,pbetaAH, palpha, pphi, ptheta, pra, prh, pua, puh, nrow = 4, ncol =3), scale=0.95) + 
  draw_label("% Increase in ICombH Relative to Baseline Usage", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

#Relative Increase from 3.382
pbetaAA <- ggplot(outputbetaAA, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="BetaAA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaHH <- ggplot(outputbetaHH, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="BetaHH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaHA <- ggplot(outputbetaHA, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1000)) +
  labs(x ="BetaHA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pbetaAH <- ggplot(outputbetaAH, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="BetaAH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

palpha <- ggplot(outputalpha, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Alpha Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pphi <- ggplot(outputphi, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Phi Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
ptheta <- ggplot(outputtheta, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="Theta Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

pra <- ggplot(outputra, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,200)) +
  labs(x ="rA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
prh <- ggplot(outputrh, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,500)) +
  labs(x ="rH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
pua <- ggplot(outputua, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="uA Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())
puh <- ggplot(outputuh, aes(x = as.numeric(ParmValue), y = as.numeric(RelInc))) + geom_line(lwd = 1.02, col ="darkblue") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) +
  labs(x ="uH Parameter Value") + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title.y=element_blank())

plot_grid(plot_grid(pbetaAA, pbetaHH, pbetaHA,pbetaAH, palpha, pphi, ptheta, pra, prh, pua, puh, nrow = 4, ncol =3), scale=0.95) + 
  draw_label("% Increase in ICombH Relative to Case Study Baseline (3.382 per 100,000)", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)
