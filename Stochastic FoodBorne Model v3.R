rm(list=ls())

# set parameters
parms = c(ra = 52^-1, rh =  6^-1, betaAA = 0.1, betaAH = 0.000001, betaHH = 0.000001, 
           betaHA = 0.00001, phi = 0.1, tau1 = 0.05, theta = 0.5)
initial=c(Sa=99, Ia=40, Ira=0, Sh = 100, Ih = 0, Irh = 0)
time.window=c(0, 20)

# initialize state and time variables and write them into output
state <- initial
time <-  time.window[1]

# define output dataframe
output <- data.frame(t=time,
                     Sa=state["Sa"], Ia=state["Ia"], Ira=state["Ira"], Sh=state["Sh"], Ih=state["Ih"], Irh=state["Irh"],
                     row.names=1)

# define how state variables S, I and R change for each process
processes <- matrix(0, nrow=15, ncol=6,
                    dimnames=list(c("HHInfSusc", "HHInfRes","HAInfSusc", "HAInfRes","AAInfSusc","AAInfRes",
                                    "AHInfSusc", "AHInfRes", "RecoveryIh", "RecoveryIrh", "RecoveryIa", "RecoveryIra",
                                    "Treatment", "Conversion", "Reversion"),
                                  c("dSh","dIh","dIrh", "dSa", "dIa", "dIra")))

#Human Transmission Events                                    
processes[1,] <- c(-1,1,0,0,0,0)
processes[2,] <- c(-1,0,1,0,0,0)
processes[3,] <- c(-1,1,0,0,0,0)
processes[4,] <- c(-1,0,1,0,0,0)
#Animal Transmission Events
processes[5,] <- c(0,0,0,-1,1,0)
processes[6,] <- c(0,0,0,-1,0,1)
processes[7,] <- c(0,0,0,-1,1,0)
processes[8,] <- c(0,0,0,-1,0,1)
#Recovery Rates
processes[9,] <- c(1,-1,0,0,0,0)
processes[10,] <- c(1,0,-1,0,0,0)
processes[11,] <- c(0,0,0,1,-1,0)
processes[12,] <- c(0,0,0,1,0,-1)
#Treatment and Transitions
processes[13,] <- c(0,0,0,1,-1,0)
processes[14,] <- c(0,0,0,0,-1,1)
processes[15,] <- c(0,0,0,0,1,-1)

# process probabilities
probabilities <- function(state,parms){
  with(as.list(c(state,parms)), {
    c(HHInfSusc = betaHH*Ih*Sh,
      HHInfRes = betaHH*Irh*Sh,
      HAInfSusc = betaHA*Ia*Sh,
      HAInfRes = betaHA*Ira*Sh,
      AAInfSusc = betaAA*Ia*Sa,
      AAInfRes = betaAA*Ira*Sa,
      AHInfSusc = betaAH*Ih*Sa, 
      AHInfRes = betaAH*Irh*Sa,
      
      RecoveryIh = rh*Ih,
      RecoveryIrh = rh*Irh,
      RecoveryIa = ra*Ia,
      RecoveryIra = ra*Ira,
      
      Treatment = tau1*Ia,
      Conversion = tau1*theta*Ia,
      Reversion = phi*Ira
      )
  })
}

while(time < time.window[2]){
  prob <- probabilities(state,parms) # calculate process probabilities for current state
  tau <- rexp(1, rate=sum(prob))   # WHEN does the next process happen? - exponentially distributed random number scaled by the sum of all process rates
  w <- sample(length(prob),1,prob=(prob)) # WHICH process happens after tau? - needs to have length 2 to choose either 1 or 2
  state <- state+processes[w,] # update states
  time <- time + tau # update time
  output <- rbind(output,c(time,state)) # write into output
}


plot(output$t,output$Sa, type="l", col="blue", ylab= "Number of Individuals", xlab="Time",lwd=2, ylim = c(0, 100))
lines(output$t,output$Ia, col="red",lwd=2)
lines(output$t,output$Ira, col="green", lwd=2)
legend(x=65,y=120, legend= c("Susceptible - Animals", "Infected - Sens - Animals", "Infected - Res - Animals"),col=c("blue","red","green"),lty=1,cex=0.9)

plot(output$t,output$Sh, type="l", col="blue", ylab= "Number of Individuals", xlab="Time",lwd=2, ylim=c(0,101))
lines(output$t,output$Ih, col="red",lwd=2)
lines(output$t,output$Irh, col="green", lwd=2)
legend(x=65,y=120, legend= c("Susceptible - Humans", "Infected - Sens - Humans", "Infected - Res - Humans"),col=c("blue","red","green"),lty=1,cex=0.9)