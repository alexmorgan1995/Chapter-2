library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("rootSolve")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Models/Github/Chapter-2/NewFits_041021/data")

# Model Functions ----------------------------------------------------------

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
  
  for(i in 1:length(grep(paste0("post_",id), list.files(), value = TRUE))) {
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

MAP <- rbind(c(colMeans(data[[1]][which(data[[1]]$group == tail(unique(data[[1]]$group),1)),][,1:6])),
             c(colMeans(data[[2]][which(data[[2]]$group == tail(unique(data[[2]]$group),1)),][,1:6])),
             c(colMeans(data[[3]][which(data[[3]]$group == tail(unique(data[[3]]$group),1)),][,1:6])),
             c(colMeans(data[[4]][which(data[[4]]$group == tail(unique(data[[4]]$group),1)),][,1:6])))

colnames(MAP) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHA")
rownames(MAP) <- c("ampbroil", "tetbroil","amppigs", "tetpigs")
