library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data")

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
datatet_pigs <- read.csv("Tet_FatPigs_Comb.csv")
datatet_hum <- read.csv("Hum_FatPigs.csv")

#Cleaning Data - Animals
datatet_pigs[,(2+5):(6+5)][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
datatet_pigs[,(2+10):(6+10)][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
datatet_pigs[,2:6][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
datatet_pigs <- datatet_pigs[!(is.na(datatet_pigs$N_2015) & is.na(datatet_pigs$N_2016) & is.na(datatet_pigs$N_2017) & 
                         is.na(datatet_pigs$N_2018) & is.na(datatet_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(datatet_pigs), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset

# NON-AGGREGATED - AMP PIGS  -------------------------------------------------------------
colnames(datatet_pigs)[12:16] <- pig_yrs

#Create dataset where each row is a different observation. 
melt_tet_pigs <- melt(datatet_pigs, id.vars = "Country", measure.vars = pig_yrs)
melt_tet_pigs$usage <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("scale_tetusage_2015", "scale_tetusage_2016", 
                                                                              "scale_tetusage_2017", "scale_tetusage_2018", "scale_tetusage_2019"))[,3]
melt_tet_pigs$N <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                        "N_2017", "N_2018", "N_2019"))[,3]
melt_tet_pigs$IsolPos <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                              "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_tet_pigs)[c(2,3)] <- c("Year", "Resistance")

#Cleaning Data - Humans
#only include countries/years which are present in the resistance dataset
datatet_hum <- datatet_hum[datatet_hum$Country %in% intersect(datatet_hum$Country, datatet_pigs$Country),]
colnames(datatet_hum)[32:37] <- as.character(2014:2019)
datatet_hum_melt <- melt(datatet_hum, id.vars = "Country", measure.vars = pig_yrs)
colnames(datatet_hum_melt)[c(2,3)] <- c("Year", "Resistance")

# Combine Human and Livestock Dataset -----------------------------------------------------------------
melt_tet_pigs$ResPropHum <- datatet_hum_melt[,3] #Obtain the melted human resistances

melt_tet_pigs <- melt_tet_pigs[!is.na(melt_tet_pigs$Resistance),]
melt_tet_pigs <- melt_tet_pigs[!is.na(melt_tet_pigs$usage),] # Remove all rows with NAs for usage and resistance

#Add 95% CIs for each datapoint
melt_tet_pigs$lower_tet <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[1]]))
melt_tet_pigs$upper_tet <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[2]]))

#Rename the columns
colnames(melt_tet_pigs) <- c("Country", "Year", "ResPropAnim", "Usage", "N", "IsolPos", "ResPropHum", "Lower_Tet", "Upper_Tet")
melt_tet_pigs$Usage <- melt_tet_pigs$Usage/1000 #Change from mg/PCU to g/PCU

ggplot(melt_tet_pigs, aes(x = Usage, y= ResPropAnim, color = Country)) + geom_point() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

#Find the average EU antibiotic usage and average human resistance for model fitting
avg_EU_usage <- mean(melt_tet_pigs$Usage)
avg_hum_res <- mean(melt_tet_pigs$ResPropHum, na.rm = TRUE)

#### Approximate Bayesian Computation - Rejection Algorithm ####

#Obtain the Resistance 
summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

#Return the sum of squares between resistance and the model output
sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- (x(data.obs) - x(model.obs))^2
  })
  return(sum(sumsquare))
}

#Compute the distances for all 3 summary statistics - this section involves running the model
computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data) {
  tauoutput <-data.frame(matrix(nrow = length(tau_range), ncol=4))
  tau_range <- append(tau_range, avg_EU_usage)
  
  for (i in 1:length(tau_range)) {
    parms2 = thetaparm
    parms2["tau"] = tau_range[i]
    
    out <- runsteady(y = init.state, func = fitmodel, times = c(0, Inf), parms = parms2)
    
    tauoutput[i,] <- c(tau_range[i],
                       ((out[[2]] + out[[3]])*(446000000))/100000,
                       out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                       out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "IncH", "ResPropAnim", "ResPropHum") 
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == avg_EU_usage,]),
           abs(tauoutput$IncH[tauoutput$tau == avg_EU_usage] - 0.593),
           abs(tauoutput$ResPropHum[tauoutput$tau == avg_EU_usage] - avg_hum_res)))
}

#Where G is the number of generations
#Function to 100% make sure the sampled particles for all parameters are non zero
prior.non.zero <- function(par){
  prod(sapply(1:6, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

#Wrapper function for all of the functions to output the distance measures and the diagnostics
#Saving of the accepted particles in each generation done within the function 
ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, data)  {
  N_ITER_list <- list()
  fit_parms <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
  thetaparm <- c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAH = 0.00001, betaHH = 0.00001)
  
  for(g in 1:G) {
    i <- 1
    dist_data <- data.frame(matrix(nrow = 1000, ncol = 3))
    N_ITER <- 1
    
    while(i <= N) {
      
      N_ITER <- N_ITER + 1
      
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 1)
        d_phi <- runif(1, min = 0, max = 0.5)
        d_kappa <- runif(1, min = 0, max = 100)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 5)
        d_betaHA <- runif(1, 0, 0.00075)
        
      } else{ 
        p <- sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        d_betaAA<-par[1]
        d_phi<-par[2]
        d_kappa<-par[3]
        d_alpha<-par[4]
        d_zeta <- par[5]
        d_betaHA <-par[6]
      }
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHA))) {
        m <- 0
        thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
                       betaHA = d_betaHA, phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
        print(dist)
        
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_food[g]) && (dist[3] <= epsilon_AMR[g]) && (!is.na(dist))) {
          # Store results
          
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHA)  
          dist_data[i,] <- dist
          # Calculate weights
          if(g==1){
            
            w.new[i] <- 1
            
          } else {
            w1<-prod(c(sapply(c(1:3,5:6), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(res.new[i,4], 1.5, 8.5))) 
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
            w.new[i] <- w1/w2
          }
          # Update counter
          print(paste0('Generation: ', g, ", particle: ", i,", weights: ", w.new[i]))
          print(dist)
          i <- i+1
        }
      }
    }
    N_ITER_list[[g]] <- list(N_ITER, dist_data)
    
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
    write.csv(res.new, file = paste("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/new/ABC_post_tetpigs_",g,".csv",sep=""), row.names=FALSE)
  }
  return(N_ITER_list)
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0)
lm.upp <- c(1, 0.5, 100, 1, 5, 0.00075) #Upper and lower bounds for the priors - for the multivariate normal dist pert kernel

# Empty matrices to store results (6 model parameters)
res.old<-matrix(ncol=6,nrow=N)
res.new<-matrix(ncol=6,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

#Thresholds 
epsilon_dist <- c(4, 3.5, 3, 2.5, 2.25, 2, 1.8, 1.75, 1.725, 1.7)
epsilon_food <- c(0.593*1, 0.593*0.8, 0.593*0.6, 0.593*0.4, 0.593*0.2, 0.593*0.1, 0.593*0.09, 0.593*0.08, 0.593*0.07, 0.593*0.06)
epsilon_AMR <- c(avg_hum_res*1, avg_hum_res*0.8, avg_hum_res*0.6, avg_hum_res*0.4, avg_hum_res*0.2, avg_hum_res*0.1, avg_hum_res*0.09, avg_hum_res*0.08, avg_hum_res*0.07, avg_hum_res*0.06)
#Run the model 
start_time <- Sys.time()

dist_save <- ABC_algorithm(N = 1000, 
              G = 10,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = melt_tet_pigs$Usage, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              data = melt_tet_pigs)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "dist_tetpigs_list.rds")

# Examining Posteriors ----------------------------------------------------

post_dist_names <- grep("ABC_post_tetpigs_",
                      list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/new"), value = TRUE)

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/new")

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
