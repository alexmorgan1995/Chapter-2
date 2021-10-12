library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data")


# Functions ---------------------------------------------------------------

sum_square_diff_dist <- function(data.obs, model.obs) {
  sumsquare <- abs((data.obs - model.obs)^2)
  return(sum(sumsquare))
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
    
    CumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    CumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), CumS, CumR))
  })
}


# Import the Data ---------------------------------------------------------

#NON-AGGREGATED DATA

#Import the data and make sure that countries and years with less than 10 samples 
amp_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_FatPigs_Years.csv")
amp_pigs[,(2+5):(6+5)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,(2+10):(6+10)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,2:6][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs <- amp_pigs[!(is.na(amp_pigs$N_2015) & is.na(amp_pigs$N_2016) & is.na(amp_pigs$N_2017) & 
                         is.na(amp_pigs$N_2018) & is.na(amp_pigs$N_2019)),]

#Rename certain columns
colnames(amp_pigs)[12:16] <- as.character(seq(2015, 2019))

#Import in the usage data
usage_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/pig_usage_years.csv")
usage_pigs <- usage_pigs[usage_pigs$Country %in% intersect(usage_pigs$Country, amp_pigs$Country),] # only include years and data points which are in the usage and resistance years

non_aggrepig <- melt(amp_pigs, id.vars = "Country", measure.vars = c("2015", "2016", "2017", "2018", "2019"))
non_aggrepig$usage <- melt(usage_pigs, id.vars = "Country", measure.vars = c("scale_ampusage_2015", "scale_ampusage_2016", 
                                                                              "scale_ampusage_2017", "scale_ampusage_2018", "scale_ampusage_2019"))[,3]
non_aggrepig$N <- melt(amp_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                        "N_2017", "N_2018", "N_2019"))[,3]
non_aggrepig$IsolPos <- melt(amp_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                              "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(non_aggrepig)[c(2,3)] <- c("Year", "Resistance")
non_aggrepig <- non_aggrepig[!is.na(non_aggrepig$N),]

non_aggrepig$lower_amp <- unlist(lapply(1:nrow(non_aggrepig), function(i) prop.test(non_aggrepig$IsolPos[i],non_aggrepig$N[i])[[6]][[1]]))
non_aggrepig$upper_amp <- unlist(lapply(1:nrow(non_aggrepig), function(i) prop.test(non_aggrepig$IsolPos[i],non_aggrepig$N[i])[[6]][[2]]))

non_aggrepig <- non_aggrepig[!is.na(non_aggrepig$usage),]

#Agregated Data

aggre_pig <- data.frame("Country" = usage_pigs$Country,
                      "usage_amp" = rowMeans(usage_pigs[,27:31], na.rm = TRUE),
                      "N" = rowSums(amp_pigs[,2:6], na.rm = TRUE),
                      "isolpos_amp" = rowSums(amp_pigs[,7:11], na.rm = TRUE))

aggre_pig$propres_amp <- aggre_pig$isolpos_amp /  aggre_pig$N

aggre_pig$lower_amp <- unlist(lapply(1:nrow(aggre_pig), function(i) prop.test(aggre_pig$isolpos_amp[i],aggre_pig$N[i])[[6]][[1]]))
aggre_pig$upper_amp <- unlist(lapply(1:nrow(aggre_pig), function(i) prop.test(aggre_pig$isolpos_amp[i],aggre_pig$N[i])[[6]][[2]]))

# Import Parameter Fits ---------------------------------------------------

nonagg_zeta <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/results_ABC_SMC_gen_tetpigs_zeta_7.csv")
nonagg_nozeta <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/results_ABC_SMC_gen_tetpigs_nozeta_5.csv")

agg_zeta <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/results_ABC_SMC_gen_tetpigs_zetaagg_10.csv")
agg_nozeta <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/results_ABC_SMC_gen_tetpigs_nozetaagg_8.csv")


parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = map_betaAA, betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = map_betaHA, phi = map_phi, kappa = map_kappa, alpha = map_alpha, zeta = 0)


# Run for each scenario ----------------------------------------------------


MAP_list <- list(
  "nonagg_zeta_MAP" = data.frame("Parameters" = colnames(nonagg_zeta),
                                 "MAP_Estimate" = colMeans(nonagg_zeta)),
  "nonagg_nozeta_MAP" = data.frame("Parameters" = colnames(nonagg_nozeta),
                                   "MAP_Estimate" = colMeans(nonagg_nozeta)),
  "agg_zeta_MAP" = data.frame("Parameters" = colnames(agg_zeta),
                              "MAP_Estimate" = colMeans(agg_zeta)),
  "agg_nozeta_MAP" = data.frame("Parameters" = colnames(agg_nozeta),
                                "MAP_Estimate" = colMeans(agg_nozeta))
)


MAP_list <- list(
  "nonagg_zeta_MAP" = map_estimate(nonagg_zeta),
  "nonagg_nozeta_MAP" = map_estimate(nonagg_nozeta),
  "agg_zeta_MAP" = map_estimate(agg_zeta),
  "agg_nozeta_MAP" = map_estimate(agg_nozeta)
)

colnames(aggre_pig)[5] <- "Resistance"
colnames(aggre_pig)[2] <- "usage"

data_list <- list(non_aggrepig,
                  aggre_pig)

scen_list <- list()

for(i in 1:4) {
  fit_parms <- MAP_list[[i]]
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAH = 0.00001, betaHH = 0.00001, zeta = 0)
  parms2[fit_parms[,1]] <-  fit_parms[,2]
  
  data <- ifelse(grepl("nonagg",names(MAP_list)[i]), data_list[1], data_list[2])
  data <- data[[1]]
  print(data$usage)
  
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  parmtau <- data$usage/1000
  
  temp <- data.frame(matrix(NA, nrow = length(parmtau), ncol=6))
  
  for(j in 1:length(parmtau)) {
    parms2["tau"] <- parmtau[j] 
    out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms2)
    temp[j,] <- c(parmtau[j],
                  ((out[[2]] + out[[3]])*(446000000))/100000,
                  (out[[2]]*(446000000))/100000,
                  (out[[3]]*(446000000))/100000,
                  out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                  out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }
  
  colnames(temp) <- c("tau", "Incidence","SuscInc", "ResInc", "IResRatA","IResRatH")

  ss <- round(sum_square_diff_dist(temp$IResRatA,data$Resistance), digits = 5)
  
  p <- ggplot(data, aes(x = usage/1000, y = Resistance, color = Country)) + geom_point() + theme_bw() + 
    geom_line(data = temp, aes(tau, IResRatA), colour = "red", size = 1.2, inherit.aes = FALSE) + 
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme(legend.position =  "bottom") +
    scale_x_continuous(expand = c(0,0)) +
    geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
    labs(title = paste0("Ampicillin Usage in Fat Pigs: 2015-2019 - ", c("non-Aggregated Data", "non-Aggregated Data", "Aggregated Data", "Aggregated Data")[i]), 
         x = "Ampicillin Usage", y = "Proportion Pigs Resistant")  +
    annotate("text",label = paste0("Sum of squares: ", ss), y = 0.95, 
             x = max(data$usage/1000, na.rm = TRUE)*0.75, size = 5)

  scen_list[[i]] <- list(temp, p)
}


for(i in 1:4) {
  ggsave(scen_list[[i]][[2]], filename = paste0(c("nonagg_zeta","nonagg_nozeta","agg_zeta","agg_nozeta")[i],".png"), 
         dpi = 300, type = "cairo", width = 8, height = 6, units = "in",
         path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")
}


# Plot the data -----------------------------------------------------------


