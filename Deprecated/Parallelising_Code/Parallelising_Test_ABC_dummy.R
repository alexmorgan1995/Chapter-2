library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel")

rm(list=ls())

#setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Models/Github/Chapter-2/NewFits_041021/data")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data")

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

#### Data Import ####

#Import Data
dataamp_pigs <- read.csv("Amp_FatPigs_Comb.csv")
dataamp_hum <- read.csv("Hum_FatPigs.csv")

#Cleaning Data - Animals
dataamp_pigs[,(2+5):(6+5)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
dataamp_pigs[,(2+10):(6+10)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
dataamp_pigs[,2:6][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
dataamp_pigs <- dataamp_pigs[!(is.na(dataamp_pigs$N_2015) & is.na(dataamp_pigs$N_2016) & is.na(dataamp_pigs$N_2017) & 
                                 is.na(dataamp_pigs$N_2018) & is.na(dataamp_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(dataamp_pigs), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset

# NON-AGGREGATED - AMP PIGS  -------------------------------------------------------------
colnames(dataamp_pigs)[12:16] <- pig_yrs

#Create dataset where each row is a different observation. 
melt_amp_pigs <- melt(dataamp_pigs, id.vars = "Country", measure.vars = pig_yrs)
melt_amp_pigs$usage <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("scale_ampusage_2015", "scale_ampusage_2016", 
                                                                                "scale_ampusage_2017", "scale_ampusage_2018", "scale_ampusage_2019"))[,3]
melt_amp_pigs$N <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                            "N_2017", "N_2018", "N_2019"))[,3]
melt_amp_pigs$IsolPos <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                                  "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_amp_pigs)[c(2,3)] <- c("Year", "Resistance")

#Cleaning Data - Humans
#only include countries/years which are present in the resistance dataset
dataamp_hum <- dataamp_hum[dataamp_hum$Country %in% intersect(dataamp_hum$Country, dataamp_pigs$Country),]
colnames(dataamp_hum)[26:31] <- as.character(2014:2019)
dataamp_hum_melt <- melt(dataamp_hum, id.vars = "Country", measure.vars = pig_yrs)
colnames(dataamp_hum_melt)[c(2,3)] <- c("Year", "Resistance")

# Combine Human and Livestock Dataset -----------------------------------------------------------------
melt_amp_pigs$ResPropHum <- dataamp_hum_melt[,3] #Obtain the melted human resistances

melt_amp_pigs <- melt_amp_pigs[!is.na(melt_amp_pigs$Resistance),]
melt_amp_pigs <- melt_amp_pigs[!is.na(melt_amp_pigs$usage),] # Remove all rows with NAs for usage and resistance

#Add 95% CIs for each datapoint
melt_amp_pigs$lower_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[1]]))
melt_amp_pigs$upper_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[2]]))

#Rename the columns
colnames(melt_amp_pigs) <- c("Country", "Year", "ResPropAnim", "Usage", "N", "IsolPos", "ResPropHum", "Lower_Amp", "Upper_Amp")
melt_amp_pigs$Usage <- melt_amp_pigs$Usage/1000 #Change from mg/PCU to g/PCU

ggplot(melt_amp_pigs, aes(x = Usage, y= ResPropAnim, color = Country)) + geom_point() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.055)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

#Find the average EU antibiotic usage and average human resistance for model fitting
avg_EU_usage <- mean(melt_amp_pigs$Usage)
avg_hum_res <- mean(melt_amp_pigs$ResPropHum, na.rm = TRUE)

#### Approximate Bayesian Computation - Rejection Algorithm ####

prior.non.zero<-function(par, lm.low, lm.upp){
  prod(sapply(1:6, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

#Return the sum of squares between resistance and the model output
sum_square_diff_dist <- function(data.obs, model.obs) {
  sumsquare <- (data.obs$ResPropAnim - model.obs$ResPropAnim)^2
  return(sum(sumsquare))
}


#Compute the distances for all 3 summary statistics - this section involves running the model
computeDistanceABC_ALEX <- function(distanceABC, fitmodel, tau_range, thetaparm, init.state, data) {
  tauoutput <- data.frame(matrix(nrow = length(tau_range), ncol=4))
  tau_range <- append(tau_range, avg_EU_usage)
  parms2 = thetaparm
  
  for (i in 1:length(tau_range)) {
    
    parms2["tau"] = tau_range[i]
    out <- runsteady(y = init.state, func = fitmodel, times = c(0, Inf), parms = parms2)
    
    tauoutput[i,] <- c(tau_range[i],
                       ((out[[2]] + out[[3]])*(446000000))/100000,
                       out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                       out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }
  
  colnames(tauoutput) <- c("tau", "IncH", "ResPropAnim", "ResPropHum") 
  return(c(distanceABC(data, tauoutput[!tauoutput$tau == avg_EU_usage,]),
           abs(tauoutput$IncH[tauoutput$tau == avg_EU_usage] - 0.593),
           abs(tauoutput$ResPropHum[tauoutput$tau == avg_EU_usage] - avg_hum_res)))
}


####

singlerun <- function(x, G, init.state, distanceABC, fitmodel, thetaparm, epsilon, 
                      tau_range, data, lm.low, lm.upp, w.old, sigma, res.old, N) {
  i <- 0
  m <- 0
  w.new <- 0
  
  while(i <= 1) {
    m <- m + 1
    if(G == 1) {
      d_betaAA <- runif(1, min = 0, max = 0.25)
      d_phi <- runif(1, min = 0, max = 0.1)
      d_kappa <- runif(1, min = 0, max = 2)
      d_alpha <- rbeta(1, 1.5, 8.5)
      d_zeta <- runif(1, 0, 1)
      d_betaHA <- runif(1, 0, 0.0005)
    } else { 
      p <- sample(seq(1,N),1,prob = w.old) # check w.old here
      par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
      d_betaAA<-par[1]
      d_phi<-par[2]
      d_kappa<-par[3]
      d_alpha<-par[4]
      d_zeta <- par[5]
      d_betaHA <-par[6]
    }
    
    new.parms = c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHA)
    
    if(prior.non.zero(new.parms, lm.low, lm.upp)) {
      
      thetaparm[c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")] <- new.parms
      
      dist <- computeDistanceABC_ALEX(distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
      
      if((dist[1] <= epsilon[["dist"]][G]) && (dist[2] <= epsilon[["food"]][G]) && (dist[3] <= epsilon[["AMR"]][G]) && (!is.na(dist))) {
        if(G==1){
          w.new <- 1
        } else {
          w1 <- prod(c(sapply(c(1:3,5:6), function(b) dunif(new.parms[b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(new.parms[4], 1.5, 8.5))) 
          w2 <- sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(new.parms, mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          w.new <- w1/w2
        }
        i <- i + 1
        return(list(dist, m, new.parms, w.new))
      }
    }
  }
}

# Large Function ----------------------------------------------------------

ABC_algorithm <- function(N, G, distanceABC, fitmodel, tau_range, init.state, data, epsilon, lm.low, lm.upp, thetaparm)  {
  out <- list()
  
  for(g in 1:G) {
    
    print(paste0("Generation ", g, " | Time: ", Sys.time()))
    
    if(g == 1) {
      sigma <- 0
      res.old <- 0
      w.old <- 0
    }
    clusterExport(cl, varlist = c("amr", "computeDistanceABC_ALEX", "prior.non.zero", "sum_square_diff_dist",
                                  "melt_amp_pigs", "avg_EU_usage", "avg_hum_res"))
    
    particles <- parLapply(cl, 1:N, 
                      singlerun, 
                      G = g, 
                      init.state = init.state,
                      distanceABC = sum_square_diff_dist,
                      fitmodel = amr, 
                      thetaparm = thetaparm, 
                      epsilon = epsilon,
                      tau_range = melt_amp_pigs$Usage,
                      data = melt_amp_pigs,
                      lm.low = lm.low,
                      lm.upp = lm.upp,
                      w.old = w.old, 
                      sigma = sigma, 
                      res.old = res.old,
                      N = N)
    
    dat_dist <- as.matrix(do.call(rbind, lapply(particles, "[[", 1)))
    dat_nruns <- do.call(sum, lapply(particles, "[[", 2))
    res.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 3)))
    w.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 4)))
    
    sigma <- cov(res.new) 
    res.old <- res.new
    w.old <- w.new/sum(w.new)
    
    out[[g]] <- list(dat_nruns, dat_dist, res.old, w.old)
    colnames(res.old) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
    write.csv(res.old, file = paste("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/Parallelising_Code/out/ABC_post_amppigs_",g,".csv",sep=""), row.names=FALSE)
  }
  return(out)
}

# Run the model -----------------------------------------------------------

cl <- makeCluster(7,type="SOCK")

clusterEvalQ(cl, {c(library("rootSolve"), library("tmvtnorm"))})

test <- ABC_algorithm(N = 1000,
              G = 10,
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = melt_amp_pigs$Usage, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              data = melt_amp_pigs, 
              epsilon = list("dist" = c(3.5, 3, 2.5, 2, 1.7, 1.5, 1.3, 1.2, 1.1, 1),
                             "food" = c(0.593*1, 0.593*0.8, 0.593*0.6, 0.593*0.4, 0.593*0.2, 0.593*0.1, 0.593*0.08, 0.593*0.07, 0.593*0.06, 0.593*0.05),
                             "AMR" = c(avg_hum_res*1, avg_hum_res*0.8, avg_hum_res*0.6, avg_hum_res*0.4, avg_hum_res*0.2, avg_hum_res*0.1, avg_hum_res*0.08, 
                                       avg_hum_res*0.07, avg_hum_res*0.06, avg_hum_res*0.05)), 
              lm.low = c(0, 0, 0, 0, 0, 0), 
              lm.upp = c(0.25, 0.1, 2, 1, 1, 0.0005), 
              thetaparm = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAH = 0.00001, betaHH = 0.00001))

stopCluster(cl)


# Test Distributions ------------------------------------------------------

post_dist_names <- grep("ABC_post_amppigs_",
                        list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/Parallelising_Code/out"), value = TRUE)

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/Parallelising_Code/out")

post_dist <- lapply(post_dist_names, read.csv)

post_dist <- mapply(cbind, post_dist, "gen" = sapply(1:length(post_dist), function(x) paste0("gen_", x)), 
                    SIMPLIFY=F)
post_dist <- do.call("rbind", post_dist)

maps_est <- map_estimate(post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][,1:6])

p_list <- list()

for(i in 1:(length(post_dist)-1)) {
  p_list[[i]] <- local ({
    name_exp <- post_dist[,c(i,7)]
    p <- ggplot(name_exp, aes(x= name_exp[,1], fill=gen)) + geom_density(alpha=.5) + 
      geom_vline(xintercept = maps_est[i,2], size = 1.2, col = "red") +
      scale_x_continuous(expand = c(0, 0), name = colnames(post_dist)[-7][i]) + 
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.text=element_text(size=14),axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    return(p)
  })
}


