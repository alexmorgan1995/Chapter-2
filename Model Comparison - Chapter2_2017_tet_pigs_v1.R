library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}


amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(m == 1) {#Model 1 - with Zeta 
      dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
        zeta*Sa*(1-alpha) - zeta*Sa 
      dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
      dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
      
      dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
      dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
      dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    }
    
    if(m == 2) {#Model 2 - no Zeta 
      dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa 
      dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa 
      dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira
      
      dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
      dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
      dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    }

    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}


#### Data Import ####
datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetrahum <- read.csv("resistanceprofHum_v1.csv")
datatetra$mgpcuuseage <- datatetra$mgpcuuseage / 1000; datatetrahum$mgpcuuseage <- datatetrahum$mgpcuuseage / 1000
datatetra$pig_tetra_sales <- datatetra$pig_tetra_sales / 1000; datatetrahum$pig_tetra_sales <- datatetrahum$pig_tetra_sales / 1000
datatetrahum$ResPropHum <- datatetrahum$ResPropHum/ 100 
datatetra <- datatetra[!datatetra$N < 10,]

mean(datatetra$ResPropHum ,na.rm = TRUE)
mean(datatetra$pig_tetra_sales ,na.rm = TRUE)

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
  tauoutput <- matrix(nrow = 0, ncol = 4)
  tau_range <- append(tau_range, 0.0122887)
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol = 4)
    parms2 = c(ra = thetaparm[["ra"]], rh =  thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], 
               betaAA = thetaparm[["betaAA"]], betaAH = thetaparm[["betaAH"]], betaHH = thetaparm[["betaHH"]], 
               betaHA = thetaparm[["betaHA"]], phi = thetaparm[["phi"]], tau = tau_range[i], kappa = thetaparm[["kappa"]], 
               alpha = thetaparm[["alpha"]], zeta = thetaparm[["zeta"]], m = thetaparm[["m"]])
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == 0.0122887,]),
           abs(tauoutput$ICombH[tauoutput$tau == 0.0122887] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == 0.0122887] - 0.35)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:5, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}



ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  for(g in 1:G) {

    i <- 1
    while(i <= N) {
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.2)
        d_phi <- runif(1, min = 0, max = 0.04)
        d_kappa <- runif(1, min = 0, max = 2)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.15)
        d_m <- ceiling(runif(1, min=0, max=2))
      } else{ 
        p <- sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp) #Check what this is - the perturbation kernel? 
        d_betaAA <- par[1]
        d_phi <- par[2]
        d_kappa <- par[3]
        d_alpha <- par[4]
        d_zeta <- par[5]
        d_m <- round(par[6],0)
      }
      
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_m))) {
        thetaparm <- c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
                       betaHA = 0.00001, phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta, m = d_m)

        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_food[g]) && (dist[3] <= epsilon_AMR[g]) && (!is.na(dist))) {
          # Store results
          
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_m)  
          # Calculate weights
          
          if(g == 1) {
            w.new[i] <- 1
            
          } else {
            
            w1 <- prod(c(sapply(c(1:3,5:6), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
                         dbeta(res.new[i,4], 1.5, 8.5))) 
            w2<-sum(sapply(1:N, function(a) w.old[a] * dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
            
            w.new[i] <- w1/w2
          }
          
          # Update counter
          print(paste0('Generation: ', g, ", particle: ", i,", weights: ", w.new[i]))
          print(dist)
          i <- i+1
        }
      }
    }#
    sigma <- cov(res.new)#create a covariance matrxi for the multivariate normal distribution used for the weights 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "m")

    write.csv(res.new, file = paste("COMP_results_ABC_SMC_gen_tet_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 1)
lm.upp <- c(0.2, 0.04, 2, 1, 0.15, 2)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=6,nrow=N)
res.new<-matrix(ncol=6,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c(2.5, 2, 1.75, 1.5, 1.25, 1, 0.9, 0.8, 0.75)
epsilon_food <- c(3.26*0.25, 3.26*0.2, 3.26*0.175, 3.26*0.15, 3.26*0.125, 3.26*0.1, 3.26*0.08, 3.26*0.06, 3.26*0.04)
epsilon_AMR <- c(0.35*0.25, 0.35*0.2, 3.26*0.175, 0.35*0.15, 0.35*0.125, 0.35*0.1, 0.35*0.08, 0.35*0.06, 3.26*0.04)

ABC_algorithm(N = 1000, 
              G = 9,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = datatetra$pig_tetra_sales, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              times = seq(0, 2000, by = 100), 
              data = datatetra)

end_time <- Sys.time(); end_time - start_time

#### Test Data ####
post_dist_m <- list()
post_dist_m[[1]] <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_1.csv", header = TRUE), "group" = "data1")
post_dist_m[[2]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_2.csv", header = TRUE), "group" = "data2")
post_dist_m[[3]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_3.csv", header = TRUE), "group" = "data3")
post_dist_m[[4]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_4.csv", header = TRUE), "group" = "data4") 
post_dist_m[[5]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_5.csv", header = TRUE), "group" = "data5") 
post_dist_m[[6]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_6.csv", header = TRUE), "group" = "data6") 
post_dist_m[[7]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_7.csv", header = TRUE), "group" = "data7") 
post_dist_m[[8]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_8.csv", header = TRUE), "group" = "data8") 
post_dist_m[[9]]  <- cbind(read.csv("COMP_results_ABC_SMC_gen_tet_9.csv", header = TRUE), "group" = "data9") 

post_dist_m_data <- data.frame("gen" = seq(1,9, by =1), "Model_1_Prop" = matrix(unlist(lapply(post_dist_m, function(x) prop.table(table(x$m))))[seq(1,17,by = 2)]))

post_dist_m_data$Model_2_Prop <- 1 - post_dist_m_data$Model_1_Prop

post_dist_m_data_melt <- as.data.frame(melt(post_dist_m_data, id.vars = c("gen"), measure.vars = c("Model_1_Prop", "Model_2_Prop")))
post_dist_m_data_melt$scen <- "tet_pigs"

write.csv(post_dist_m_data_melt, file = "comp_tetpigs.csv")

comp_plot <- ggplot(post_dist_m_data_melt, aes(x = gen, y = value, fill = as.factor(variable))) + geom_bar(position = "fill", stat="identity") + 
  scale_x_continuous(breaks = seq(1,10), expand = c(0, 0), name = expression(paste("Generation"))) + 
  scale_y_continuous(expand = c(0, 0), name = bquote("Proportion of " ~ italic(m) ~ "particles accepted (n = 1000)")) +
  theme(legend.text=element_text(size=14),  axis.text = element_text(size=14),
        axis.title.y=element_text(size=14), axis.title.x = element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.title=element_text(size=14), title= element_text(size= 15)) +
  scale_fill_manual(values = c("black","grey"),name = "Model", labels = c("1 (Zeta)", "2 (No Zeta)")) + 
  labs(title = "Model Comparison (Tetracycline usage in Fattening Pigs)")

ggsave(comp_plot, filename = "comp_plot_tetpigs.png", dpi = 300, type = "cairo", width = 10, height = 9, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

