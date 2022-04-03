library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel")

rm(list=ls())

setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Models/Github/Chapter-2/NewFits_041021/data")
# Model FUnctions ---------------------------------------------------------
#Model

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      (0.5*zeta)*Sa*(1-alpha) - (0.5*zeta)*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + (0.5*zeta)*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + (0.5*zeta)*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    
    CumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    CumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), CumS, CumR))
  })
}


# Run the Model -----------------------------------------------------------



parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0.01)
init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)

#Make the cluster 

cl <- makeCluster(8,type="SOCK")
clusterEvalQ(cl, {library("deSolve")})

system.time(test <- parLapply(cl, 1:500, function(x, amr) {
  init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
  tau_t <- runif(1, 0, 0.05)
  parms["tau"] <- tau_t 
  out <- ode(y = init.state, func = amr, times = seq(0, 1000), parms = parms)
  return(out)
}, amr))

stopCluster(cl)

system.time(test <- lapply(1:10, function(x, amr) {
  init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
  tau_t <- runif(1, 0, 0.05)
  parms["tau"] <- tau_t 
  out <- ode(y = init.state, func = amr, times = seq(0, 1000), parms = parms)
  return(out)
}, amr)
)



# Test --------------------------------------------------------------------

cl <- makeCluster(14,type="SOCK")
clusterEvalQ(cl, {library("rootSolve")})

system.time(test <- parLapply(cl, 1:100, function(x, amr) {
  i <- 0
  while(i <= 1) {
    init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
    parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
               betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
    tau_t <- runif(1, 0, 100)
    parms["tau"] <- tau_t 
    out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
    if(out$y[["Ira"]] <  0.9)
    {
      i <- i +1
    }
  }
  return(parms["tau"])
}, amr))

stopCluster(cl)



system.time(lapply(1:100, function(x, amr) {
  i <- 0
  while(i <= 1) {
    init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
    parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
               betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
    tau_t <- runif(1, 0, 100)
    parms["tau"] <- tau_t 
    out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
    if(out$y[["Ira"]] <  0.9)
    {
      i <- i +1
    }
  }
  return(parms["tau"])
}, amr))


detectCores(all.tests = FALSE, logical = TRUE)

lapply(1:2, function(x, amr) {
  i <- 0
  while(i <= 1) {
    init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
    parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
               betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
    tau_t <- runif(1, 0, 100)
    parms["tau"] <- tau_t 
    out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
    if(out$y[["Ira"]] <  0.9)
    {
      i <- i +1
    }
  }
  return(parms["tau"])
}, amr)


# Timings -----------------------------------------------------------------

cl <- makeCluster(14,type="SOCK")
clusterEvalQ(cl, {library("rootSolve")})

timings <- matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100) {
  time_cluster <- system.time(test <- parLapply(cl, 1:i, function(x, amr) {
    i <- 0
    while(i <= 1) {
      init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
      parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
      tau_t <- runif(1, 0, 100)
      parms["tau"] <- tau_t 
      out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
      if(out$y[["Ira"]] <  0.9)
      {
        i <- i +1
      }
    }
    return(parms["tau"])
  }, amr))
  time_normal <- system.time(lapply(1:i, function(x, amr) {
    i <- 0
    while(i <= 1) {
      init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
      parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
      tau_t <- runif(1, 0, 100)
      parms["tau"] <- tau_t 
      out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
      if(out$y[["Ira"]] <  0.9)
      {
        i <- i +1
      }
    }
    return(parms["tau"])
  }, amr))
  
  timings[i,] <-c(i, time_cluster[[3]], time_normal[[3]])
}

stopCluster(cl)

timings <- data.frame(timings[!is.na(timings[,2]),])

timings_m <- melt(timings, id.vars = c("X1"), measure.vars = c("X2","X3"))


ggplot(timings_m, aes(x = X1, y = value, col = variable)) +geom_line() + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))

time_normal <- system.time(lapply(1:100, function(x, amr) {
  i <- 0
  while(i <= 1) {
    init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
    parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.3, betaAH = 0.00001, betaHH = 0.00001, 
               betaHA = 0.001, phi = 0.05, kappa = 0.8, alpha = 0.3, zeta = 1, tau = 0)
    tau_t <- runif(1, 0, 100)
    parms["tau"] <- tau_t 
    out <- runsteady(y = init.state, func = amr, times = c(0, Inf), parms = parms)
    if(out$y[["Ira"]] <  0.9)
    {
      i <- i +1
    }
  }
  return(parms["tau"])
}, amr))
