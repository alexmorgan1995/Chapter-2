library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data")

#Model
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
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
datatet_pigs[,(2+5):(6+5)][datatet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
datatet_pigs[,(2+10):(6+10)][datatet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
datatet_pigs[,2:6][datatet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
datatet_pigs <- datatet_pigs[!(is.na(datatet_pigs$N_2015) & is.na(datatet_pigs$N_2016) & is.na(datatet_pigs$N_2017) & 
                         is.na(datatet_pigs$N_2018) & is.na(datatet_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(datatet_pigs), value = TRUE)) #Find the number of years included 

#Cleaning Data - Humans
datatet_hum <- datatet_hum[datatet_hum$Country %in% intersect(datatet_hum$Country, datatet_pigs$Country),]
datatet_hum <- datatet_hum[,c(1,grep(paste(pig_yrs,collapse="|"), colnames(datatet_hum)))]
datatet_hum <- datatet_hum[,-grep("AMP", colnames(datatet_hum))] # Remove AMP from dataframe

#Combine the Data 

data_pig_fit <- data.frame("Country" = datatet_hum$Country,
                           "N" = rowSums(datatet_pigs[,2:6], na.rm = TRUE),
                           "N_PosIsol" = rowSums(datatet_pigs[,7:11], na.rm = TRUE),
                           "Usage" = rowMeans(datatet_pigs[,grep("usage_20", colnames(datatet_pigs))], na.rm = TRUE)/1000,
                           "ResPropHum" = rowSums(datatet_hum[,7:11], na.rm = TRUE)/
                            rowSums(datatet_hum[,2:6], na.rm = TRUE),
                           "ResPropAnim" = rowSums(datatet_pigs[,7:11], na.rm = TRUE)/
                             rowSums(datatet_pigs[,2:6], na.rm = TRUE)) 

data_pig_fit$lower <- unlist(lapply(1:nrow(data_pig_fit), function(i) prop.test(data_pig_fit$N_PosIsol[i], data_pig_fit$N[i])[[6]][[1]]))
data_pig_fit$upper <- unlist(lapply(1:nrow(data_pig_fit), function(i) prop.test(data_pig_fit$N_PosIsol[i], data_pig_fit$N[i])[[6]][[2]]))

ggplot(data_pig_fit, aes(x = Usage, y= ResPropAnim)) + geom_point() +
  geom_text(aes(x = Usage, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05, inherit.aes = TRUE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.055)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

avg_EU_usage <- mean(data_pig_fit$Usage)
avg_hum_res <- mean(data_pig_fit$ResPropHum, na.rm = TRUE)

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

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data) {
  tauoutput <-data.frame(matrix(nrow = length(tau_range), ncol=4))
  tau_range <- append(tau_range, avg_EU_usage)
  
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol = 4)
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

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:6, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

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
        d_betaAA <- runif(1, min = 0, max = 1.5)
        d_phi <- runif(1, min = 0, max = 1)
        d_kappa <- runif(1, min = 0, max = 250)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 1)
        d_betaHA <- runif(1, 0, 0.0015)
        
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
    write.csv(res.new, file = paste("results_ABC_SMC_gen_tetpigs_",g,".csv",sep=""), row.names=FALSE)
  }
  return(N_ITER_list)
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0)
lm.upp <- c(1.5, 1, 250, 1, 1, 0.0015)

# Empty matrices to store results (6 model parameters)
res.old<-matrix(ncol=6,nrow=N)
res.new<-matrix(ncol=6,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c(2.5, 2, 1.5, 1.25, 1, 0.9, 0.8, 0.75, 0.7, 0.65)
epsilon_food <- c(0.593*1, 0.593*0.75, 0.593*0.5, 0.593*0.25, 0.593*0.2, 0.593*0.15, 0.593*0.1, 0.593*0.075, 0.593*0.065, 0.593*0.05)
epsilon_AMR <- c(avg_hum_res*1, avg_hum_res*0.75, avg_hum_res*0.5, avg_hum_res*0.25, avg_hum_res*0.2, avg_hum_res*0.15, avg_hum_res*0.1, avg_hum_res*0.075, avg_hum_res*0.065, avg_hum_res*0.05)

dist_save <- ABC_algorithm(N = 1000, 
              G = 10,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = data_pig_fit$Usage, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              data = data_pig_fit)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "dist_tetpigs_list.rds")

#### Test Data ####
data1 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_1.csv", header = TRUE), "group" = "data1")
data2 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_2.csv", header = TRUE), "group" = "data2")
data3 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_3.csv", header = TRUE), "group" = "data3")
data4 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_4.csv", header = TRUE), "group" = "data4") 
data5 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_5.csv", header = TRUE), "group" = "data5") 
data6 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_6.csv", header = TRUE), "group" = "data6")
data7 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_7.csv", header = TRUE), "group" = "data7")
data8 <- cbind(read.csv("results_ABC_SMC_gen_tetpigs_8.csv", header = TRUE), "group" = "data8")


#Plotting the Distributions

testphi <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "phi"); testphi$group <- factor(testphi$group, levels = unique(testphi$group))
testkappa <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "kappa"); testkappa$group <- factor(testkappa$group, levels = unique(testkappa$group))
testbetaAA <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "betaAA"); testbetaAA$group <- factor(testbetaAA$group, levels = unique(testbetaAA$group))
testalpha <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "alpha"); testalpha$group <- factor(testalpha$group, levels = unique(testalpha$group))
testzeta <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "zeta"); testzeta$group <- factor(testzeta$group, levels = unique(testzeta$group))
testbetaHA <- melt(rbind(data4, data5, data6, data7, data8), id.vars = "group", measure.vars = "betaHA"); testbetaHA$group <- factor(testbetaHA$group, levels = unique(testbetaHA$group))

p1 <- ggplot(testphi, aes(x=value, fill=group)) + geom_density(alpha=.5) + 
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p2 <- ggplot(testkappa, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p3<- ggplot(testbetaAA, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p4 <- ggplot(testalpha, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels =  c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p5 <- ggplot(testzeta, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels =c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p6<- ggplot(testbetaHA, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Human Transmission (", beta[HA], ")"))) + 
  scale_y_continuous(expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 4", "Generation 5", "Generation 6", "Generation 7", "Generation 8")) +
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


#### Plotting ####

plot <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol =2, 
                  labels = c("A","B","C","D","E", "F"), font.label = c(size = 20), common.legend = TRUE, legend = "bottom")

#### Testing the Model #### 

map_parms <- map_estimate(data8[,1:6])

parmtau <- seq(0, 0.05,by=0.001)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
temp <- data.frame(matrix(NA, nrow = length(parmtau), ncol=6))

for (i in 1:length(parmtau)) {
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_parms[map_parms$Parameter == "betaAA", 2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = map_parms[map_parms$Parameter == "betaHA", 2], phi = map_parms[map_parms$Parameter == "phi", 2], 
             kappa = map_parms[map_parms$Parameter == "kappa", 2], alpha = map_parms[map_parms$Parameter == "alpha", 2], tau = parmtau[i], 
             zeta = map_parms[map_parms$Parameter == "zeta", 2])
  
  out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms2)
  temp[i,] <- c(parmtau[i],
                ((out[[2]] + out[[3]])*(446000000))/100000,
                (out[[2]]*(446000000))/100000,
                (out[[3]]*(446000000))/100000,
                out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
}

colnames(temp)[1:6] <- c("tau", "Incidence","SuscInc", "ResInc", "IResRatA","IResRatH")

plotdata <- melt(temp,
                 id.vars = c("tau"), measure.vars = c("ResInc","SuscInc")) 

ggplot(plotdata[plotdata$tau <= 0.02,], aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = avg_EU_usage, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,1.6), expand = c(0, 0))  + 
  geom_text(label= c(round(temp$IResRatH[temp$tau <= 0.02],digits = 2), rep("",length(temp$IResRatH[temp$tau <= 0.02]))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Tetracycline Antibiotic Usage (g/PCU)", y = "Tetracycline-Resistant Fattening Pig Carriage") 


ggplot(temp, aes(x = tau, y = IResRatA)) + geom_line() +  theme_bw() +
  geom_text(data = data_pig_fit, aes(x = Usage, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05, inherit.aes = FALSE) +
  geom_point(data = data_pig_fit, aes(x = Usage, y = ResPropAnim), inherit.aes = FALSE)  +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))+
  labs(x ="Tetracycline Antibiotic Usage (g/PCU)", y = "Ampicllin-Resistant Fattening Pig Carriage") 
