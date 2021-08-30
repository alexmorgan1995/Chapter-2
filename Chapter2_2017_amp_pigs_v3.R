library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
#setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

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
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

#### Data Import ####
dataamp <- read.csv("resistanceprofAnim_amp.csv")

dataamp$mgpcuuseage <- dataamp$mgpcuuseage / 1000
dataamp$pig_amp_sales <- dataamp$pig_amp_sales / 1000
dataamp <- dataamp[!dataamp$N < 10,]

dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

mean(dataamp$hum_res_amp ,na.rm = TRUE)
mean(dataamp$pig_amp_sales ,na.rm = TRUE)

ggplot(dataamp, aes(x = pig_amp_sales, y= ResPropAnim))  + geom_point() +
  geom_text(aes(x = pig_amp_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE)

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
  tau_range <- append(tau_range, 0.01156391)
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol=4)
    parms2 = c(ra = thetaparm[["ra"]], rh =  thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], 
               betaAA = thetaparm[["betaAA"]], betaAH = thetaparm[["betaAH"]], betaHH = thetaparm[["betaHH"]], 
               betaHA = thetaparm[["betaHA"]], phi = thetaparm[["phi"]], tau = tau_range[i], kappa = thetaparm[["kappa"]], 
               alpha = thetaparm[["alpha"]], zeta = thetaparm[["zeta"]])
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == 0.01156391,]),
           abs(tauoutput$ICombH[tauoutput$tau == 0.01156391] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == 0.01156391] - 0.3108)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:5, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  N_ITER_LIST <- list()
  for(g in 1:G) {
    i <- 1
    dist_data <- data.frame(matrix(nrow = 1000, ncol = 3))
    N_ITER <- 1
    
    while(i <= N) {
      
      N_ITER <- N_ITER + 1
      
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.2)
        d_phi <- runif(1, min = 0, max = 0.5)
        d_kappa <- runif(1, min = 0, max = 50)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.25)
      } else{ 
        p <- sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        d_betaAA<-par[1]
        d_phi<-par[2]
        d_kappa<-par[3]
        d_alpha<-par[4]
        d_zeta <-par[5]
      }
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta))) {
        m <- 0
        thetaparm <- c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
                       betaHA = 0.00001, phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        #print(dist)
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_food[g]) && (dist[3] <= epsilon_AMR[g]) && (!is.na(dist))) {
          # Store results
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta)  
          dist_data[i,] <- dist
          print(dist)
          # Calculate weights

          if(g==1){
            
            w.new[i] <- 1
            
          } else {
            w1<-prod(c(sapply(c(1:3,5), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(res.new[i,4], 1.5, 8.5))) 
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
            w.new[i] <- w1/w2
          }
          
          # Update counter
          print(paste0('Generation: ', g, ", particle: ", i,", weights: ", w.new[i]))
          i <- i+1
        }
      }
    }#
    N_ITER_list[[g]] <- list(N_ITER, dist_data)
    
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta")
    write.csv(res.new, file = paste("results_ABC_SMC_gen_amp_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
  return(N_ITER_list)
}


N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0)
lm.upp <- c(0.2, 0.5, 50, 1, 0.25)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=5,nrow=N)
res.new<-matrix(ncol=5,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c(2.5, 2, 1.5, 1.25, 1, 0.9, 0.8, 0.7, 0.6, 0.55)
epsilon_food <- c(3.26*0.25, 3.26*0.2, 3.26*0.15, 3.26*0.125, 3.26*0.1, 3.26*0.075, 3.26*0.05, 3.26*0.03, 3.26*0.01, 3.26*0.005)
epsilon_AMR <- c(0.3108*0.25, 0.3108*0.2, 0.3108*0.15, 0.3108*0.125, 0.3108*0.1, 0.3108*0.075, 0.3108*0.05, 0.3108*0.03, 0.3108*0.01, 0.3108*0.005)

dist_save <- ABC_algorithm(N = 1000, 
              G = 10,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = dataamp$pig_amp_sales, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              times = seq(0, 2000, by = 50), 
              data = dataamp)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "dist_amppigs_list.rds")
 
#### Test Data ####
data1 <- cbind(read.csv("results_ABC_SMC_gen_amp_1.csv", header = TRUE), "group" = "data1")
data2 <- cbind(read.csv("results_ABC_SMC_gen_amp_2.csv", header = TRUE), "group" = "data2")
data3 <- cbind(read.csv("results_ABC_SMC_gen_amp_3.csv", header = TRUE), "group" = "data3")
data4 <- cbind(read.csv("results_ABC_SMC_gen_amp_4.csv", header = TRUE), "group" = "data4") 
data5 <- cbind(read.csv("results_ABC_SMC_gen_amp_5.csv", header = TRUE), "group" = "data5") 
data6 <- cbind(read.csv("results_ABC_SMC_gen_amp_6.csv", header = TRUE), "group" = "data6")
data7 <- cbind(read.csv("results_ABC_SMC_gen_amp_7.csv", header = TRUE), "group" = "data7")
data8 <- cbind(read.csv("results_ABC_SMC_gen_amp_8.csv", header = TRUE), "group" = "data8")
data9 <- cbind(read.csv("results_ABC_SMC_gen_amp_9.csv", header = TRUE), "group" = "data9") 
data10 <- cbind(read.csv("results_ABC_SMC_gen_amp_10.csv", header = TRUE), "group" = "data10") 

map_phi <- map_estimate(data10[,"phi"], precision = 20) 
map_kappa <- map_estimate(data10[,"kappa"], precision = 20) 
map_betaAA <- map_estimate(data10[,"betaAA"], precision = 20) 
map_alpha <- map_estimate(data10[,"alpha"], precision = 20) 
map_zeta <- map_estimate(data10[,"zeta"], precision = 20) 

map_phi <- mean(data10[,"phi"]) 
map_kappa <- mean(data10[,"kappa"]) 
map_betaAA <- mean(data10[,"betaAA"]) 
map_alpha <- mean(data10[,"alpha"]) 
map_zeta <- mean(data10[,"zeta"]) 

#Plotting the Distributions

testphi <- melt(rbind(data6, data7, data8, data9, data10), id.vars = "group",measure.vars = "phi"); testphi$group <- factor(testphi$group, levels = unique(testphi$group))
testkappa <- melt(rbind(data6, data7, data8, data9,data10), id.vars = "group",measure.vars = "kappa"); testkappa$group <- factor(testkappa$group, levels = unique(testkappa$group))
testbetaAA <- melt(rbind(data6, data7, data8, data9, data10), id.vars = "group",measure.vars = "betaAA"); testbetaAA$group <- factor(testbetaAA$group, levels = unique(testbetaAA$group))
testalpha <- melt(rbind(data6, data7, data8, data9, data10), id.vars = "group",measure.vars = "alpha"); testalpha$group <- factor(testalpha$group, levels = unique(testalpha$group))
testzeta <- melt(rbind(data6, data7, data8, data9, data10), id.vars = "group",measure.vars = "zeta"); testzeta$group <- factor(testzeta$group, levels = unique(testzeta$group))

p1 <- ggplot(testphi, aes(x=value, fill=group)) + geom_density(alpha=.5) + 
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")"))) + 
  scale_y_continuous(limits = c(0,100), expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 6", "Generation 7", "Generation8", "Generation 9", "Generation 10"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p2 <- ggplot(testkappa, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")"))) + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 6", "Generation 7", "Generation8", "Generation 9", "Generation 10"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p3<- ggplot(testbetaAA, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")"))) + 
  scale_y_continuous(limits = c(0,60),expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 6", "Generation 7", "Generation8", "Generation 9", "Generation 10"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p4 <- ggplot(testalpha, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")"))) + 
  scale_y_continuous(limits = c(0,9),expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels =c("Generation 6", "Generation 7", "Generation8", "Generation 9", "Generation 10"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p5 <- ggplot(testzeta, aes(x=value, fill=group)) + geom_density(alpha=.5)+
  scale_x_continuous(expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) + 
  scale_y_continuous(limits = c(0,100),expand = c(0, 0), name = "") +
  labs(fill = NULL) + scale_fill_discrete(labels = c("Generation 6", "Generation 7", "Generation8", "Generation 9", "Generation 10"))+
  theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


#### Plotting ####

plot <- ggarrange(p1, p2, p3, p4,p5, nrow = 3, ncol =2, 
                  labels = c("A","B","C","D","E"), font.label = c(size = 20), common.legend = TRUE, legend = "bottom")

ggsave(plot, filename = "ABCSMC_salm_amp_pigs.png", dpi = 300, type = "cairo", width = 13, height = 11, units = "in")

# Pairs Plot --------------------------------------------------------------

pairs_data <- data10[data10$group == "data10",1:4]

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) + 
    stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours=viridis::viridis(100))
  p
}

pairs_amp <- GGally::ggpairs(pairs_data, lower=list(continuous=my_fn)) + theme_bw()

ggsave(pairs_amp, filename = "pairs_plot_amp.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

#### Testing the Model #### 

parmtau <- c(seq(0,0.035,by=0.001), 0.01156391)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = map_phi, kappa = map_kappa, alpha = map_alpha, tau = parmtau[i], zeta = map_zeta)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:7] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

output2 <- output1
output2[,2:5] <- output2[,2:5]*100000 #Scaling the prevalence (per 100,000)

plot_ly(output2, x= ~tau, y = ~InfHumans, type = "bar", name = "Antibiotic-Sensitive Infection (Human)") %>%
  add_trace(y= ~ResInfHumans, name = "Antibiotic-Resistant Infection (Human)", textfont = list(size = 30)) %>% 
  layout(yaxis = list(title = "Proportion of Infected Humans (ICombH) per 100,000", range = c(0,6), showline = TRUE),
         xaxis = list(title = "Livestock Antibiotic Usage (g/PCU)"),
         legend = list(orientation = "v", x = 0.6, y=1), showlegend = T,
         barmode = "stack", annotations = list(x = ~tau, y = ~ICombH, text = ~IResRat, yanchor = "bottom", showarrow = FALSE, textangle = 310,
                                               xshift =3))

#### Showing the Fit with Data #### 

parmtau <- seq(0,0.035, by = 0.001)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame()
times <- seq(0, 200000, by = 100)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = map_phi, kappa = map_kappa, alpha = map_alpha, tau = parmtau[i], zeta = map_zeta)
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:7] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

output2 <- output1
output2[,2:5] <- output2[,2:5]*100000 #Scaling the prevalence (per 100,000)

dataamp <- read.csv("resistanceprofAnim_amp.csv")

dataamp$mgpcuuseage <- dataamp$mgpcuuseage / 1000
dataamp$pig_amp_sales <- dataamp$pig_amp_sales / 1000
dataamp <- dataamp[!dataamp$N < 10,]

dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

ggplot(dataamp, aes(x = pig_amp_sales, y= ResPropAnim))  + geom_point() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.03)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = output2, aes(x = tau, y= IResRatA), col = "red", size = 1.02)
