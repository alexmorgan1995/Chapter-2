library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2"); library("ggrepel")
library("tidyr"); library("bayestestR"); library("profvis"); library("tmvtnorm")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Basic Model  ####
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

#Function to round large D.P numbers
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 1000, by = 1)

parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = (0.0709), betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0141, theta = 0.0303, alpha = 0.232, tau = 0.032)
out <- data.frame(ode(y = init, func = amr, times = times, parms = parms2))
meltedout <- melt(out, id.vars = c("time"), measure.vars = c("Isa", "Ira"))

#For the Animal Population
ggplot(data = meltedout, aes(x = time, y = value, col = variable)) + 
  geom_line(size = 1.02) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Data Import ####
datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetrahum <- read.csv("resistanceprofHum_v1.csv")
datatetra$mgpcuuseage <- datatetra$mgpcuuseage / 1000; datatetrahum$mgpcuuseage <- datatetrahum$mgpcuuseage / 1000
datatetra$pig_tetra_sales <- datatetra$pig_tetra_sales / 1000; datatetrahum$pig_tetra_sales <- datatetrahum$pig_tetra_sales / 1000
datatetrahum$ResPropHum <- datatetrahum$ResPropHum/ 100 
datatetra <- datatetra[!datatetra$N < 5,]

ggplot()  + geom_point(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim)) +
  geom_text(data = datatetra, aes(x = pig_tetra_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

ggplot()  + geom_point(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum)) +
  geom_text(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Human Carriage")

#### Approximate Bayesian Computation - Rejection Algorithm ####

summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol=4)
  tau_range <- append(tau_range, 0.0106)
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol=4)
    parms2 = c(ra = thetaparm[["ra"]], rh =  thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], 
               betaAA = thetaparm[["betaAA"]], betaAH = thetaparm[["betaAH"]], betaHH = thetaparm[["betaHH"]], 
               betaHA = thetaparm[["betaHA"]], phi = thetaparm[["phi"]], tau = tau_range[i], theta = thetaparm[["theta"]], 
               alpha = thetaparm[["alpha"]])
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == 0.0106,]),
           abs(tauoutput$ICombH[tauoutput$tau == 0.0106] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == 0.0106] - 0.32)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0)
lm.upp <- c(0.2, 0.03, 0.4, 1)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=4,nrow=N)
res.new<-matrix(ncol=4,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c()
epsilon_food <- c()
epsilon_AMR <- c()

#Test Environment

#Where G is the number of generations
par <- data_ABC1[1,]
prior.non.zero<-function(par){
  prod(sapply(1:4, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}
prior.non.zero(par)
rtmvnorm(1,mean=data_ABC1[1,], sigma=cov(data_ABC1), lower=lm.low, upper=lm.upp)


ABC_algorithm <- function(N, epsilon, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  for(g in 1:G) {
    i <- 0
    
    while(i < N) {
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.2)
        d_phi <- runif(1, min = 0, max = 0.03)
        d_theta <- runif(1, min = 0, max = 0.4)
        d_alpha <- rbeta(1, 1.5, 8.5)
      } else{ 
        p<-sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        d_phi<-par[1]
        d_theta<-par[2]
        d_betaAA<-par[3]
        d_alpha<-par[4]
      }
      if(prior.non.zero(c(d_phi,d_theta,d_betaAA,d_alpha))) {
        m <- 0
        thetaparm <- c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
                       betaHA = 0.00001, phi = d_phi, theta = d_theta, alpha = d_alpha)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_food[g]) && (dist[3] <= epsilon_AMR[g]) && (!is.na(dist))) {
          m <- m+1
        }
        
        ####
        if (m>0){
          # Store results
          res.new[i,]<-c(d_phi, d_theta, d_betaAA, d_alpha)  
          # Calculate weights
          w1<-prod(sapply(1:5, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b]))
                   #how will this work for a beta distribution?
                   #Have an if function here if alpha then do dbeta
                   #else (theta, betaAA or phi then do a dunif)
                   ) 
          
          if(g==1){
            
            w2<-1
            
          } else {
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          }
          w.new[i] <- (m/n)*w1/w2
          # Update counter
          i <- i+1
          print(paste0('Generation: ', g, ", particle: ", i))
        }
      }
    }#
    sigma <- cov(res.new) 
    res.old<-res.new
    w.old<-w.new/sum(w.new)
    
    write.csv(res.new, file = paste("results_ABC_SMC_gen_",g,".csv",sep=""), row.names=FALSE)
    
    ####
  }
}


data_ABC <- ABC_algorithm(N = 3, 
                          epsilon = c(1.5, 3.26*0.15, 0.32*0.15),
                          sum.stats = summarystatprev, 
                          distanceABC = sum_square_diff_dist, 
                          fitmodel = amr, 
                          tau_range = datatetra$pig_tetra_sales, 
                          init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                          times = seq(0, 2000, by = 50), 
                          data = datatetra)

end_time <- Sys.time(); end_time - start_time

colnames(data_ABC) <- c("phi", "theta", "d_alpha", "d_betaAA")
