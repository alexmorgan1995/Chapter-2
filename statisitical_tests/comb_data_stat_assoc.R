library(reshape2); library(ggplot2); library(ggpubr); library(betareg)
rm(list=ls())

# Import and Clean Data - Broilers (Ampicillin and Tetracycline) -----------
# Ampicillin Resistance
amp_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_Broil_Years.csv")
amp_broil[,(2+4):(5+4)][amp_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
amp_broil[,2:5][amp_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
amp_broil <- amp_broil[!(is.na(amp_broil$N_2014) & is.na(amp_broil$N_2016) & is.na(amp_broil$N_2017) & 
                           is.na(amp_broil$N_2018)),]

colnames(amp_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

# Tetracycline Resistance
tet_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Tet_Broil_Years.csv")
tet_broil[,(2+4):(5+4)][tet_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
tet_broil[,2:5][tet_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
tet_broil <- tet_broil[!(is.na(tet_broil$N_2014) & is.na(tet_broil$N_2016) & is.na(tet_broil$N_2017) & 
                           is.na(tet_broil$N_2018)),]

colnames(tet_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

#Check that we have the same countries for tetracycline and ampicillin usage
tet_broil$Country %in% amp_broil$Country

#Broiler Usage
usage_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/broil_usage_years.csv")
usage_broil <- usage_broil[usage_broil$Country %in% intersect(usage_broil$Country, amp_broil$Country),]

# Import and Clean Data - Fattening Pigs (Ampicillin and Tetracycline) -----------
# Ampicillin Resistance
amp_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_FatPigs_Years.csv")
amp_pigs[,(2+5):(6+5)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,(2+10):(6+10)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,2:6][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs <- amp_pigs[!(is.na(amp_pigs$N_2015) & is.na(amp_pigs$N_2016) & is.na(amp_pigs$N_2017) & 
                         is.na(amp_pigs$N_2018) & is.na(amp_pigs$N_2019)),]

colnames(amp_pigs)[12:16] <- as.character(seq(2015, 2019))

# Tetracycline Resistance
tet_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Tet_FatPigs_Years.csv")
tet_pigs[,(2+5):(6+5)][tet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
tet_pigs[,(2+10):(6+10)][tet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
tet_pigs[,2:6][tet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with

tet_pigs <- tet_pigs[!(is.na(tet_pigs$N_2015) & is.na(tet_pigs$N_2016) & is.na(tet_pigs$N_2017) & 
                         is.na(tet_pigs$N_2018) & is.na(tet_pigs$N_2019)),]

colnames(tet_pigs)[12:16] <- as.character(seq(2015, 2019))

#Check that we have the same countries for tetracycline and ampicillin usage
tet_pigs$Country %in% amp_pigs$Country

#Broiler Usage
usage_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/pig_usage_years.csv")
usage_pigs <- usage_pigs[usage_pigs$Country %in% intersect(usage_pigs$Country, tet_pigs$Country),]

# Stat Testing BROILERS -------------------------------------------------------------

#Create the Combined Dataset 

broil_rel <- data.frame("country" = usage_broil$Country,
                        "usage_amp" = rowMeans(usage_broil[,22:25], na.rm = TRUE),
                        "usage_tet" = rowMeans(usage_broil[,26:29], na.rm = TRUE),
                        "N" = rowSums(amp_broil[,2:5], na.rm = TRUE),
                        "isolpos_amp" = rowSums(amp_broil[,6:9], na.rm = TRUE),
                        "isolpos_tet" = rowSums(tet_broil[,6:9], na.rm = TRUE))


broil_rel$propres_amp <- broil_rel$isolpos_amp /  broil_rel$N
broil_rel$propres_tet <- broil_rel$isolpos_tet /  broil_rel$N

broil_rel$propres_amp_trans <- ((broil_rel$propres_amp*(nrow(broil_rel)-1)) + 0.5)/nrow(broil_rel)
broil_rel$propres_tet_trans <- ((broil_rel$propres_tet*(nrow(broil_rel)-1)) + 0.5)/nrow(broil_rel)

broil_rel$lower_amp <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_amp[i],broil_rel$N[i])[[6]][[1]]))
broil_rel$upper_amp <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_amp[i],broil_rel$N[i])[[6]][[2]]))
broil_rel$lower_tet <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_tet[i],broil_rel$N[i])[[6]][[1]]))
broil_rel$upper_tet <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_tet[i],broil_rel$N[i])[[6]][[2]]))

#Stat Test
broil_stat_amp <- betareg(propres_amp_trans ~ usage_amp, broil_rel)
broil_stat_tet <- betareg(propres_tet_trans ~ usage_tet, data = broil_rel)
summary(broil_stat_amp)
summary(broil_stat_tet)

broil_amp <- ggplot(broil_rel, aes(x = usage_amp, y = propres_amp)) + geom_point() + theme_bw() +  
  geom_line(aes(y = predict(broil_stat_amp, broil_rel), x = usage_amp), colour = "red", size = 1.2) + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
  labs(title = "Ampicillin Usage in Broilers: 2014-2018", x = "Ampicillin Usage", y = "Proportion Broilers Resistant")

broil_tet <- ggplot(broil_rel, aes(x = usage_tet, y = propres_tet)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(broil_stat_tet, broil_rel), x = usage_tet), colour = "red", size = 1.2) + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet),  size=0.5, width = 0)+ 
  labs(title = "Tetracycline Usage in Broilers: 2014-2018", x = "Tetracycline Usage", y = "Proportion Broilers Resistant")

# Stat Testing - FATTENING PIGS -------------------------------------------------------------

pig_rel <- data.frame("Country" = usage_pigs$Country,
                      "usage_amp" = rowMeans(usage_pigs[,27:31], na.rm = TRUE),
                      "usage_tet" = rowMeans(usage_pigs[,32:36], na.rm = TRUE),
                      "N" = rowSums(amp_pigs[,2:6], na.rm = TRUE),
                      "isolpos_amp" = rowSums(amp_pigs[,7:11], na.rm = TRUE),
                      "isolpos_tet" = rowSums(tet_pigs[,7:11], na.rm = TRUE))

pig_rel$propres_amp <- pig_rel$isolpos_amp /  pig_rel$N
pig_rel$propres_tet <- pig_rel$isolpos_tet /  pig_rel$N

pig_rel$propres_amp_trans <- ((pig_rel$propres_amp*(nrow(pig_rel)-1)) + 0.5)/nrow(pig_rel)
pig_rel$propres_tet_trans <- ((pig_rel$propres_tet*(nrow(pig_rel)-1)) + 0.5)/nrow(pig_rel)

pig_rel$lower_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[1]]))
pig_rel$upper_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[2]]))
pig_rel$lower_tet <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_tet[i],pig_rel$N[i])[[6]][[1]]))
pig_rel$upper_tet <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_tet[i],pig_rel$N[i])[[6]][[2]]))

#Stat Test
pig_stat_amp <- betareg(propres_amp_trans ~ usage_amp, pig_rel)
pig_stat_tet <- betareg(propres_tet_trans ~ usage_tet, data = pig_rel)
summary(pig_stat_amp)
summary(pig_stat_tet)

pig_amp <- ggplot(pig_rel, aes(x = usage_amp, y = propres_amp)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_amp, pig_rel), x = usage_amp), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
  labs(title = "Ampicillin Usage in Fat Pigs: 2015-2019", x = "Ampicillin Usage", y = "Proportion Pigs Resistant")

pig_tet <- ggplot(pig_rel, aes(x = usage_tet, y = propres_tet)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_tet, pig_rel), x = usage_tet), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet),  size=0.5, width = 0) + 
  labs(title = "Tetracycline Usage in Fat Pigs: 2015-2019", x = "Tetracycline Usage", y = "Proportion Pigs Resistant")


# Comb Plots --------------------------------------------------------------

broilers <- ggarrange(broil_amp, broil_tet, ncol = 2, nrow = 1)
pigs <- ggarrange(pig_amp, pig_tet, ncol = 2, nrow = 1)

ggsave(broilers, filename = "broil_comb_stat.png", dpi = 300, type = "cairo", width = 14, height = 5, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

ggsave(pigs, filename = "pigs_comb_stat.png", dpi = 300, type = "cairo", width = 14, height = 5, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

