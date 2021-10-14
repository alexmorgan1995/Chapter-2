library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("rootSolve")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data")

# Model Functions ----------------------------------------------------------

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model ODEs
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

#Importing in the Datasets

#for(j in 1:4) {
#  fit_names <- unique(gsub(".{6}$","",grep("ABC_post_",list.files(), value = TRUE)))
#  nam <- grep(fit_names[1], list.files(), value = TRUE)
#}

import <- function(id) {
  data <- data.frame(matrix(ncol = 6, nrow = 0))
  
  for(i in 1:length(grep(id, list.files(), value = TRUE))) {
    test  <- cbind(read.csv(paste0("ABC_post_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

# Data Import -------------------------------------------------------------
#Import of Posterior Distributions
data <- list(import("tetbroil"), import("ampbroil"), import("tetpigs"), import("amppigs"))
lapply(1:4, function(x) data[[x]]$group = factor(data[[x]]$group, levels = unique(data[[x]]$group)))

#Import of Fitting Data
#Ampicillin in Broiler
dataamp_broil <- read.csv("Amp_Broil_Comb.csv")
dataamp_broil[,(2+4):(5+4)][dataamp_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
dataamp_broil[,(2+8):(5+8)][dataamp_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
dataamp_broil[,2:5][dataamp_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
dataamp_broil <- dataamp_broil[!(is.na(dataamp_broil$N_2014) & is.na(dataamp_broil$N_2016) & is.na(dataamp_broil$N_2017) & 
                                   is.na(dataamp_broil$N_2018)),]
broil_yrs <- sub("N_", "", grep("N_20",colnames(dataamp_broil), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset
colnames(dataamp_broil)[10:13] <- broil_yrs

melt_amp_broil <- melt(dataamp_broil, id.vars = "Country", measure.vars = broil_yrs)
melt_amp_broil$usage <- melt(dataamp_broil, id.vars = "Country", measure.vars = c("scale_ampusage_2014", "scale_ampusage_2016", 
                                                                                  "scale_ampusage_2017", "scale_ampusage_2018"))[,3]
melt_amp_broil$N <- melt(dataamp_broil, id.vars = "Country", measure.vars = c("N_2014", "N_2016", 
                                                                              "N_2017", "N_2018"))[,3]
melt_amp_broil$IsolPos <- melt(dataamp_broil, id.vars = "Country", measure.vars = c("PosIsol_2014", "PosIsol_2016", 
                                                                                    "PosIsol_2017", "PosIsol_2018"))[,3]
colnames(melt_amp_broil)[c(2,3)] <- c("Year", "Resistance"); rm(dataamp_broil)
melt_amp_broil <- melt_amp_broil[!(is.na(melt_amp_broil$Resistance) | is.na(melt_amp_broil$usage)),]

#Tetracycline in Broiler
datatet_broil <- read.csv("Tet_Broil_Comb.csv")
datatet_broil[,(2+4):(5+4)][datatet_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
datatet_broil[,(2+8):(5+8)][datatet_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
datatet_broil[,2:5][datatet_broil[,2:5] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
datatet_broil <- datatet_broil[!(is.na(datatet_broil$N_2014) & is.na(datatet_broil$N_2016) & is.na(datatet_broil$N_2017) & 
                                   is.na(datatet_broil$N_2018)),]
broil_yrs <- sub("N_", "", grep("N_20",colnames(datatet_broil), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset
colnames(datatet_broil)[10:13] <- broil_yrs

melt_tet_broil <- melt(datatet_broil, id.vars = "Country", measure.vars = broil_yrs)
melt_tet_broil$usage <- melt(datatet_broil, id.vars = "Country", measure.vars = c("scale_tetusage_2014", "scale_tetusage_2016", 
                                                                                  "scale_tetusage_2017", "scale_tetusage_2018"))[,3]
melt_tet_broil$N <- melt(datatet_broil, id.vars = "Country", measure.vars = c("N_2014", "N_2016", 
                                                                              "N_2017", "N_2018"))[,3]
melt_tet_broil$IsolPos <- melt(datatet_broil, id.vars = "Country", measure.vars = c("PosIsol_2014", "PosIsol_2016", 
                                                                                    "PosIsol_2017", "PosIsol_2018"))[,3]
colnames(melt_tet_broil)[c(2,3)] <- c("Year", "Resistance"); rm(datatet_broil)
melt_tet_broil <- melt_tet_broil[!(is.na(melt_tet_broil$Resistance) | is.na(melt_tet_broil$usage)),]

#Ampicillin in Pigs
dataamp_pigs <- read.csv("Amp_FatPigs_Comb.csv")
dataamp_pigs[,(2+5):(6+5)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
dataamp_pigs[,(2+10):(6+10)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
dataamp_pigs[,2:6][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
dataamp_pigs <- dataamp_pigs[!(is.na(dataamp_pigs$N_2015) & is.na(dataamp_pigs$N_2016) & is.na(dataamp_pigs$N_2017) & 
                                 is.na(dataamp_pigs$N_2018) & is.na(dataamp_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(dataamp_pigs), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset
colnames(dataamp_pigs)[12:16] <- pig_yrs

melt_amp_pigs <- melt(dataamp_pigs, id.vars = "Country", measure.vars = pig_yrs)
melt_amp_pigs$usage <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("scale_ampusage_2015", "scale_ampusage_2016", 
                                                                                "scale_ampusage_2017", "scale_ampusage_2018", "scale_ampusage_2019"))[,3]
melt_amp_pigs$N <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                            "N_2017", "N_2018", "N_2019"))[,3]
melt_amp_pigs$IsolPos <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                                  "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_amp_pigs)[c(2,3)] <- c("Year", "Resistance"); rm(dataamp_pigs)
melt_amp_pigs <- melt_amp_pigs[!(is.na(melt_amp_pigs$Resistance) | is.na(melt_amp_pigs$usage)),]

#Tetracycline in Pigs
datatet_pigs <- read.csv("Tet_FatPigs_Comb.csv")
datatet_pigs[,(2+5):(6+5)][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
datatet_pigs[,(2+10):(6+10)][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
datatet_pigs[,2:6][datatet_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
datatet_pigs <- datatet_pigs[!(is.na(datatet_pigs$N_2015) & is.na(datatet_pigs$N_2016) & is.na(datatet_pigs$N_2017) & 
                                 is.na(datatet_pigs$N_2018) & is.na(datatet_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(datatet_pigs), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset
colnames(datatet_pigs)[12:16] <- pig_yrs

melt_tet_pigs <- melt(datatet_pigs, id.vars = "Country", measure.vars = pig_yrs)
melt_tet_pigs$usage <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("scale_tetusage_2015", "scale_tetusage_2016", 
                                                                                "scale_tetusage_2017", "scale_tetusage_2018", "scale_tetusage_2019"))[,3]
melt_tet_pigs$N <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                            "N_2017", "N_2018", "N_2019"))[,3]
melt_tet_pigs$IsolPos <- melt(datatet_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                                  "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_tet_pigs)[c(2,3)] <- c("Year", "Resistance"); rm(datatet_pigs)

melt_tet_pigs <- melt_tet_pigs[!(is.na(melt_tet_pigs$Resistance) | is.na(melt_tet_pigs$usage)),]

# Create 95% CIs ----------------------------------------------------------
melt_tet_broil$lower <- unlist(lapply(1:nrow(melt_tet_broil), function(i) prop.test(melt_tet_broil$IsolPos[i],melt_tet_broil$N[i])[[6]][[1]]))
melt_tet_broil$upper <- unlist(lapply(1:nrow(melt_tet_broil), function(i) prop.test(melt_tet_broil$IsolPos[i],melt_tet_broil$N[i])[[6]][[2]]))

melt_amp_broil$lower <- unlist(lapply(1:nrow(melt_amp_broil), function(i) prop.test(melt_amp_broil$IsolPos[i],melt_amp_broil$N[i])[[6]][[1]]))
melt_amp_broil$upper <- unlist(lapply(1:nrow(melt_amp_broil), function(i) prop.test(melt_amp_broil$IsolPos[i],melt_amp_broil$N[i])[[6]][[2]]))

melt_amp_pigs$lower <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[1]]))
melt_amp_pigs$upper <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[2]]))

melt_tet_pigs$lower <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[1]]))
melt_tet_pigs$upper <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[2]]))

#Obtain the MAPs for each dataset

MAP <- rbind(c(map_estimate(data[[1]][which(data[[1]]$group == tail(unique(data[[1]]$group),1)),][,1:6])[,2]),
             c(map_estimate(data[[2]][which(data[[2]]$group == tail(unique(data[[2]]$group),1)),][,1:6])[,2]),
             c(map_estimate(data[[3]][which(data[[3]]$group == tail(unique(data[[3]]$group),1)),][,1:6])[,2]),
             c(map_estimate(data[[4]][which(data[[4]]$group == tail(unique(data[[4]]$group),1)),][,1:6])[,2]))

colnames(MAP) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
rownames(MAP) <- c("ampbroil", "tetbroil","amppigs", "tetpigs")

# ABC-SMC Posterior -------------------------------------------------------

data_bind <- do.call(rbind,data)

plotlist1 <- list()

data_test <- data_bind[data_bind$fit == unique(data_bind$fit)[1],]
data_test[data_test$group == tail(unique(data_test$group),1),]

for(j in 1:length(unique(data_bind$fit))) {
  
  plotlist2 <- list()
  
  data_test <- data_bind[data_bind$fit == unique(data_bind$fit)[j],]
  
  for (i in 1:length(colnames(MAP))) { # Loop over loop.vector

    dens <- density(data_test[data_test$group == tail(unique(data_test$group),1),][,i])
    
    plotlist2[[i]] <- local({
      i = i
      p1 <- ggplot(data_test, aes(x=get(colnames(MAP)[i]), fill=group)) + geom_density(alpha=.5) + theme_bw()  +
        scale_fill_discrete(labels = unique(data_test$group)) +
        theme(legend.text=element_text(size=10), axis.text.x=element_text(size=10),axis.ticks.y=element_blank(), axis.text.y=element_blank(),
              axis.title.y=element_text(size=10), axis.title.x= element_text(size=10), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
              plot.title = element_text(size = 12, vjust = 3, hjust = 0.5, face = "bold"))
      
      if(colnames(MAP)[i] == "phi") {
        p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Resistance Reversion (", phi, ")"))) + 
          scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "kappa") {
        p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Recovery (", kappa, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "betaAA") {
        p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")))+ scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ") +
          labs(fill = NULL, title = c("Ampicillin Sales in Broilers", "Tetracycline Sales in Broilers", 
                                      "Ampicillin Sales in Pigs", "Tetracycline Sales in Pigs")[j]) 
      }
      if(colnames(MAP)[i] == "alpha") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.8),expand = c(0, 0), name = expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      if(colnames(MAP)[i] == "zeta") {
        p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      if(colnames(MAP)[i] == "betaHA") {
        p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Human Transmission (", beta[HA], ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      
      p1 <- p1 + geom_vline(xintercept = MAP[j,i], size  = 1.2, color = "red", alpha = 0.5)
      
      return(p1)
      
    })
    
    print(paste0("Plot Parameter: ",colnames(MAP)[i], " | Data: ", unique(data_bind$fit)[j] ))
  }
  plotlist1[[unique(data_bind$fit)[j]]] <- plotlist2
}

abc <- ggarrange(plotlist1[[1]][[1]], plotlist1[[2]][[1]], plotlist1[[3]][[1]], plotlist1[[4]][[1]],
                 plotlist1[[1]][[2]], plotlist1[[2]][[2]], plotlist1[[3]][[2]], plotlist1[[4]][[2]],
                 plotlist1[[1]][[3]], plotlist1[[2]][[3]], plotlist1[[3]][[3]], plotlist1[[4]][[3]],
                 plotlist1[[1]][[4]], plotlist1[[2]][[4]], plotlist1[[3]][[4]], plotlist1[[4]][[4]],
                 plotlist1[[1]][[5]], plotlist1[[2]][[5]], plotlist1[[3]][[5]], plotlist1[[4]][[5]],
                 plotlist1[[1]][[6]], plotlist1[[2]][[6]], plotlist1[[3]][[6]], plotlist1[[4]][[6]],
                 nrow = 6, ncol =4, 
                 labels =  c("A","B","C","D",
                             "","","", "",
                             "","","", "",
                             "","","", "",
                             "","","", "",
                             "","","", ""),
                 font.label = c(size = 20), common.legend = TRUE, legend = "bottom",
                 align = "hv", vjust = 1.05)

ggsave(abc, filename = "ABC_SMC_Post_v2.png", dpi = 300, type = "cairo", width = 13, height = 12, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

# Model Fit with Data - Ribbon -----------------------------------------------------

start_time <- Sys.time()

data_bind <- do.call(rbind,data)

parmtau <- seq(0,0.08, by = 0.002)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
ribbon_final <- data.frame()

for(j in 1:nrow(MAP)) {
  output1 <- data.frame(matrix(NA, nrow = length(parmtau), ncol = 8))
  output_ribbon <- data.frame()
  
  ribbondata <- data_bind[data_bind$fit == unique(data_bind$fit)[j],] 
  ribbondata <- ribbondata[ribbondata$group == tail(unique(ribbondata$group),1),] 
  
  for (i in 1:length(parmtau)) {
    
    temp_ribbon <- data.frame(matrix(NA, nrow = 1000, ncol=4))
    
    if(rownames(MAP)[j] == "ampbroil") {
      parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = MAP[j,"betaHA"], phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])} 
    if(rownames(MAP)[j] == "tetbroil") {
      parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = MAP[j,"betaHA"], phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])} 
    if(rownames(MAP)[j] == "amppigs") {
      parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = MAP[j,"betaHA"], phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])} 
    if(rownames(MAP)[j] == "tetpigs") {
      parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = MAP[j,"betaHA"], phi = MAP[j,"phi"], kappa = MAP[j,"kappa"], alpha = MAP[j,"alpha"], tau = parmtau[i],
                 zeta = MAP[j,"zeta"])} 
    
    out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms2)
    output1[i,1] <- parmtau[i]
    output1[i,2] <- out$y["Sh"]
    output1[i,3] <- (out[[2]]*(446000000))/100000
    output1[i,4] <- (out[[3]]*(446000000))/100000
    output1[i,5] <- ((out[[2]] + out[[3]])*(446000000))/100000
    output1[i,6] <- out[[3]] / (out[[2]] + out[[3]])
    output1[i,7] <- out$y["Ira"] / (out$y["Isa"] + out$y["Ira"])
    output1[i,8] <- rownames(MAP)[j]
    
    for(z in 1:1000) {
      if(rownames(MAP)[j] == "ampbroil") {
        parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])} 
      if(rownames(MAP)[j] == "tetbroil") {
        parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])} 
      if(rownames(MAP)[j] == "amppigs") {
        parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])} 
      if(rownames(MAP)[j] == "tetpigs") {
        parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = ribbondata[z,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                   betaHA = (0.00001), phi =  ribbondata[z,"phi"], kappa = ribbondata[z,"kappa"], alpha = ribbondata[z,"alpha"], tau = parmtau[i],
                   zeta = ribbondata[z,"zeta"])}
      
      out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms2)
      temp_ribbon[z,1] <- parmtau[i]
      temp_ribbon[z,2] <- z
      temp_ribbon[z,3] <- out$y["Ira"] / (out$y["Isa"] + out$y["Ira"])
      temp_ribbon[z,4] <- rownames(MAP)[j]
      
      print(paste0(temp_ribbon[z,4], ", tau: ", temp_ribbon[z,1], ", ", (z/1000)*100, "%"))
      
    }
    output_ribbon <- rbind.data.frame(output_ribbon, temp_ribbon)
  }
  icombhdata <- rbind(icombhdata, output1)
  ribbon_final <- rbind(ribbon_final, output_ribbon)
}

colnames(icombhdata)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
colnames(ribbon_final)[1:4] <- c("tau","particle","IResRatA", "group")

HDI_ribbon <- data.frame()
for(j in 1:length(unique(ribbon_final$group))) {
  for(i in 1:length(unique(ribbon_final$tau))) {
    HDI_ribbon <- rbind(HDI_ribbon, 
                        data.frame("tau" = unique(ribbon_final$tau)[i],
                                   "lowHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[2]],
                                   "highHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[3]],
                                   "scen" = unique(ribbon_final$group)[j]))
  }
}

end_time <- Sys.time(); end_time - start_time

amp_broil <- ggplot(melt_amp_broil, aes(x = usage/1000, y= Resistance))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,(max(melt_amp_broil$usage)/1000)*1.05), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Broiler Poultry Ampicillin Sales (g/PCU)", y = "Ampicillin-Resistant Broiler Poultry Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "ampbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.1)+
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "ampbroil",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

tet_broil <- ggplot(melt_tet_broil, aes(x = usage/1000, y= Resistance))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,(max(melt_tet_broil$usage)/1000)*1.05), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Broiler Poultry Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Broiler Poultry Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "tetbroil",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  geom_line(data = icombhdata[icombhdata$group == "tetbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.1)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

amp_pig <- ggplot(melt_amp_pigs, aes(x = usage/1000, y= Resistance))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,(max(melt_amp_pigs$usage)/1000)*1.05),expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Ampicillin Sales (g/PCU)", y = "Ampicillin-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "amppigs",], aes(x = tau, y= IResRatA), col = "red", size = 1.1) +
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "amppigs",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

tet_pig <- ggplot(melt_tet_pigs, aes(x = usage/1000, y= Resistance))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,(max(melt_tet_pigs$usage)/1000)*1.05), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "tetpigs",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  geom_line(data = icombhdata[icombhdata$group == "tetpigs",], aes(x = tau, y= IResRatA), col = "red", size = 1.1) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

fits <- ggarrange(amp_broil, tet_broil,amp_pig, tet_pig, nrow = 2, ncol = 2, align = "h", labels = c("A","B", "C", "D"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(fits, filename = "Model_Fits_v1.png", dpi = 300, type = "cairo", width = 14, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

# Tau and ICombH Base Plot - Requires you to run previous section ----------------------------------------------------------

icombhlist <- list()

#amp_broil, tet_broil, amp_pigs, tet_pigs

averagesales <- sapply(list(melt_amp_broil$usage, melt_tet_broil$usage, melt_amp_pigs$usage,melt_tet_pigs$usage), mean)/1000

for(i in 1:4){
  icombhlist[[i]] <- local({
    
    plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[i],],
                     id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 
    
    p1 <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
      geom_vline(xintercept = averagesales[i], alpha = 0.3, size = 2) + 
      geom_col(color = "black",position= "stack", width  = 0.004) + scale_x_continuous(expand = c(0, 0.0005)) + 
      scale_y_continuous(limits = c(0,2.5), expand = c(0, 0))  + 
      geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[i]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
                position = "stack", angle = 45) +
      theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm')) + 
      scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) 
    
    if(unique(icombhdata$group)[i] == "ampbroil") {
      p1 <- p1 + labs(x ="Ampicillin Usage in Broiler Poultry (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    
    if(unique(icombhdata$group)[i] == "tetbroil") {
      p1 <- p1 + labs(x ="Tetracycline Usage in Broiler Poultry (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    
    if(unique(icombhdata$group)[i] == "amppigs") {
      p1 <- p1 + labs(x ="Ampicillin Usage in Fattening Pigs (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    
    if(unique(icombhdata$group)[i] == "tetpigs") {
      p1 <- p1 + labs(x ="Tetracycline Usage in Fattening Pigs (g/PCU)", y = "Inf Humans (per 100,000)")  
    }
    
    return(p1)
  })
}

icombh <- ggarrange(icombhlist[[1]], icombhlist[[2]], 
                    icombhlist[[3]], icombhlist[[4]], 
                    nrow = 2, ncol = 2, align = "v", labels = c("A","B", "C", "D"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(icombh, filename = "Icombh.png", dpi = 300, type = "cairo", width = 8, height = 11, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

ggsave(icombh, filename = "Icombh_poster.png", dpi = 300, type = "cairo", width = 10, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# HDI for parameter estimates ---------------------------------------------

#Data is defined at the beginning 

HDI_parms <- list()

for(i in 1:length(unique(data$fit))) {
  temp_data <- data[data$group == "data10" & data$fit == c("tet", "amp", "broil")[i],]
  HDI_parms[[i]] <- data.frame("parm" = c("phi", "kappa", "betaAA", "alpha", "zeta"),
                               "mean_parm" = unlist(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                                temp_data[["alpha"]], temp_data[["zeta"]]), mean)),
                               "lowHDI" = sapply(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                             temp_data[["alpha"]], temp_data[["zeta"]]), hdi), "[[", 2),
                               "highHDI" = sapply(lapply(list(temp_data[["phi"]], temp_data[["kappa"]], temp_data[["betaAA"]], 
                                                              temp_data[["alpha"]], temp_data[["zeta"]]), hdi), "[[", 3))
}; HDI_parms

# Distances for Fits ------------------------------------------------------

summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, parms, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol=4)
  for (i in 1:length(tau_range)) {
    parms1 <- parms
    parms1[["tau"]] <- tau_range[[i]]
    temp <- matrix(NA, nrow = 1, ncol=4)
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms1)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput)))
}


datatetra <- read.csv("resistanceprofAnim_v1.csv"); datatetra <- datatetra[!datatetra$N < 10,]
colnames(datatetra)[12] <- "usage"; datatetra$usage <- datatetra$usage/1000 
databroil <- read.csv("salm_broilers_2018.csv"); databroil <- databroil[!databroil$N < 10,]
colnames(databroil)[9] <- "usage"; databroil$usage <- databroil$usage/1000 
dataamp <- read.csv("resistanceprofAnim_amp.csv"); dataamp <- dataamp[!dataamp$N < 10,]
colnames(dataamp)[9] <- "usage"; dataamp$usage <- dataamp$usage/1000 

parmstet_pigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPtet","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPtet","phi"] , 
                  kappa = MAP["MAPtet","kappa"], alpha = MAP["MAPtet","alpha"] , tau = 0, zeta = MAP["MAPtet","zeta"])

parmstet_broil = c(ra = 0, rh = (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP["MAPbroil","betaAA"], 
                   betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPbroil","phi"], 
                   kappa = MAP["MAPbroil","kappa"], alpha = MAP["MAPbroil","alpha"] , tau = 0, 
                   zeta = MAP["MAPbroil","zeta"])

parms_amppigs = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP["MAPamp","betaAA"], 
                  betaAH = 0.00001, betaHH = 0.00001, betaHA = (0.00001), phi = MAP["MAPamp","phi"], 
                  kappa = MAP["MAPamp","kappa"], alpha = MAP["MAPamp","alpha"], tau = 0, zeta = MAP["MAPamp","zeta"])

for(i in 1:3) {
  casestudy <- list(datatetra, databroil, dataamp)[[i]]
  parms <- list(parmstet_pigs, parmstet_broil, parms_amppigs)[[i]]
  #print(parms)
  test<- computeDistanceABC_ALEX(sum.stats = summarystatprev, 
                          distanceABC = sum_square_diff_dist, 
                          fitmodel = amr, 
                          tau_range = casestudy$usage, 
                          parms = parms,
                          init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                          times = seq(0, 2000, by = 100), 
                          data = casestudy)
  print(test)
}

#Identify the Exact Values of the Fit
tauoutput <- matrix(nrow = 0, ncol=4)

averagesales <- c(0.0122887, 0.01156391, 0.006666697)

exploredtau <- c(datatetra$usage, averagesales[1], 0)

init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times = seq(0, 2000, by = 100)

for (i in 1:length(exploredtau)) {
  parms1 <- parmstet_pigs
  parms1[["tau"]] <- exploredtau[[i]]
  temp <- matrix(NA, nrow = 1, ncol=4)
  out <- ode(y = init.state, func = amr, times = times, parms = parms1)
  temp[1,1] <- exploredtau[i]
  temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
  temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
  temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
  tauoutput <- rbind(tauoutput, temp)
}

tauoutput <- data.frame(tauoutput)
colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  

tauoutput$ICombH[tauoutput$tau == 0]/tauoutput$ICombH[tauoutput$tau == averagesales[1]]

# Plotting Prior Distribution ---------------------------------------------

prior.data <- setNames(data.frame(runif(1000, min = 0, max = 0.2), runif(1000, min = 0, max = 0.04), runif(1000, min = 0, max = 2),
                                  "p_alpha" <- rbeta(1000, 1.5, 8.5), runif(1000, min = 0, max = 0.15)), 
                       c("p_betaAA", "p_phi", "p_kappa", "p_alpha", "p_zeta"))

prior.plot.list <- list()

for (i in 1:length(colnames(prior.data))) { # Loop over loop.vector
  
  data <- data.frame("data" = prior.data[,colnames(prior.data)[i]])
  p1 <- ggplot(data, aes(x=data)) + geom_histogram(fill="darkgrey", bins = 10, size = 0.7) + theme_bw() +
    theme(legend.text=element_text(size=14), axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
          plot.title = element_text(size = 14, vjust = 3, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  
  max.y <- max(ggplot_build(p1)$data[[1]][1])*1.2
  
    
  if(colnames(prior.data)[i] == "p_phi") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Resistance Reversion (", phi, ")"))) +
        labs(fill = NULL, title = "") + 
        scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency") 
    }
    if(colnames(prior.data)[i] == "p_kappa") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Recovery (", kappa, ")"))) +
        labs(fill = NULL, title = "") + 
        scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency") 
    }
    if(colnames(prior.data)[i] == "p_betaAA") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")))+
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y), expand = c(0, 0), name = "Frequency")
    }
    if(colnames(prior.data)[i] == "p_alpha") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")"))) +
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y),expand = c(0, 0), name = "Frequency")
    }
    if(colnames(prior.data)[i] == "p_zeta") {
      p1 <- p1 + scale_x_continuous(expand = c(0, 0), name = expression(paste("Background Infection Rate (", zeta, ")"))) +
        labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0,max.y),expand = c(0, 0), name = "Frequency")
    }
  prior.plot.list[[i]] <- p1
}

comb.prior.plot <- ggarrange(prior.plot.list[[1]],prior.plot.list[[2]],prior.plot.list[[3]],prior.plot.list[[4]],prior.plot.list[[5]],
                             ncol = 1,nrow = 5)

ggsave(comb.prior.plot, filename = "prior_parms.png", dpi = 300, type = "cairo", width = 8, height = 13, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")

# Supplementary Plots (kappa Responsible for Co-Existence) -----------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
output1 <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)


for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
  
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[1,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = MAP[1,"phi"], kappa = 0, alpha = 0, tau = parmtau[i],
             zeta = MAP[1,"zeta"])
  
  out <- ode(y = init, func = amr, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),5]) 
  temp[1,3] <- rounding(out[nrow(out),6]) 
  temp[1,4] <- rounding(out[nrow(out),7])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
  temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  temp[1,8] <- rownames(MAP)[1]
  print(temp[1,3])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
output1[,2:5] <- output1[,2:5]*100000 #Scaling the prevalence (per 100,000)

plotdata <- melt(output1,
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

averagesales <- 0.0122887

p1 <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = averagesales, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,7), expand = c(0, 0))  + 
  geom_text(label= c(round(output1$IResRat,digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Tetracycline Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)")  


ggsave(p1, filename = "Icombh_nokappa.png", dpi = 300, type = "cairo", width = 7, height = 4, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")
