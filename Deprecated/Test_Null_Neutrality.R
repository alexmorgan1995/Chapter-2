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
    
    dSa = ra*(Isa + Ira) - (betaAA*Isa*Sa) - (betaAA*Ira*Sa) - zeta*Sa - zeta*Sa 
    dIsa = betaAA*Isa*Sa - ra*Isa + zeta*Sa
    dIra = betaAA*Ira*Sa - ra*Ira  + zeta*Sa
    
    return(list(c(dSa,dIsa,dIra)))
  })
}



# Running the Model -------------------------------------------------------
#The hypotehtis behind null neutrality is that if there were no differences between the strains, then both strains would be treated equally
#So for example with the strain with 1% resistant and 99% sensitive and allowing dual colonisation - there wouldn't be a random increase to 50% carriage
#We can run this experiment with our model. 

#BASELINE
init <- c(Sa=0.5, Isa=0.3, Ira=0.2, Sh=0, Ish=0, Irh=0)
times <- seq(0, 1000, by = 1)

parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = 0.02944741, betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.01310852, theta = 1.130831, alpha = 0.4285830, tau = 0.0123,
           zeta = 0.04968898)
out <- data.frame(ode(y = init, func = amr, times = times, parms = parms2))

test <- melt(out, id.vars = "time", measure.vars = c("Ira", "Isa"))

ggplot(test, aes(x = time, y = value, color = variable)) + geom_line()

#SIMILAR STRAINS

init <- c(Sa=0.5, Isa=0.3, Ira=0.2)
times <- seq(0, 10000, by = 1)

parms2 = c(ra = 0.1, rh =  (0), ua = 0, uh = 0, betaAA = 0.5, betaAH = 0, betaHH = 0, 
           betaHA = 0, phi = 0, theta = 0, alpha = 1, tau = 0,
           zeta = .1)
out <- data.frame(ode(y = init, func = amr, times = times, parms = parms2))

test <- melt(out, id.vars = "time", measure.vars = c("Ira", "Isa"))

ggplot(test, aes(x = time, y = value, color = variable)) + geom_line() + scale_y_continuous(limits = c(0, 0.001))

       