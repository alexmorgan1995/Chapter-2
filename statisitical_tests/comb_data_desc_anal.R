library(reshape2); library(ggplot2); library(ggpubr); library(betareg)


# Ampicillin in Fattening Pigs ------------------------------------------

amp_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Amp_FatPigs_Years.csv")
amp_pigs[,2:6][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with

amp_pigs <- amp_pigs[!(is.na(amp_pigs$N_2015) & is.na(amp_pigs$N_2016) & is.na(amp_pigs$N_2017) & 
                         is.na(amp_pigs$N_2018) & is.na(amp_pigs$N_2019)),]

colnames(amp_pigs)[12:16] <- as.character(seq(2015, 2019))

amp_pigs_melt <- melt(amp_pigs, id.vars = "Country", measure.vars = c("2015","2016","2017","2018","2019"))

amp_pigs_melt$N <- unlist(amp_pigs[2:6])
amp_pigs_melt$PosIsol <- unlist(amp_pigs[7:11])
amp_pigs_melt$value[is.na(amp_pigs_melt$N)] <- NA
amp_pigs_melt$PosIsol[is.na(amp_pigs_melt$N)] <- NA

amp_pigs_melt$lower <- 0; amp_pigs_melt$upper <- 0

for(i in 1:nrow(amp_pigs_melt)) {
  if(!is.na(amp_pigs_melt[i, "N"])) {
    amp_pigs_melt[i,"lower"] <- prop.test(amp_pigs_melt$PosIsol[i],amp_pigs_melt$N[i])[[6]][[1]]
    amp_pigs_melt[i,"upper"] <- prop.test(amp_pigs_melt$PosIsol[i],amp_pigs_melt$N[i])[[6]][[2]]
  } else{
    amp_pigs_melt[i,"lower"] <- NA
    amp_pigs_melt[i,"upper"] <- NA
  }
}

p1 <- ggplot(amp_pigs_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Country),  size=0.5, width = 0) + coord_flip() +
  facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Proportion of Ampicillin Resistant Isolates") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Ampicillin in Broiler Poultry ------------------------------------------

amp_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Amp_Broil_Years.csv")
amp_broil[,2:5][amp_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with

amp_broil <- amp_broil[!(is.na(amp_broil$N_2014) & is.na(amp_broil$N_2016) & is.na(amp_broil$N_2017) & 
                         is.na(amp_broil$N_2018)),]

colnames(amp_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

amp_broil_melt <- melt(amp_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))

amp_broil_melt$N <- unlist(amp_broil[2:5])
amp_broil_melt$PosIsol <- unlist(amp_broil[6:9])
amp_broil_melt$value[is.na(amp_broil_melt$N)] <- NA
amp_broil_melt$PosIsol[is.na(amp_broil_melt$N)] <- NA

amp_broil_melt$lower <- 0; amp_broil_melt$upper <- 0

for(i in 1:nrow(amp_broil_melt)) {
  if(!is.na(amp_broil_melt[i, "N"])) {
    amp_broil_melt[i,"lower"] <- prop.test(amp_broil_melt$PosIsol[i],amp_broil_melt$N[i])[[6]][[1]]
    amp_broil_melt[i,"upper"] <- prop.test(amp_broil_melt$PosIsol[i],amp_broil_melt$N[i])[[6]][[2]]
  } else{
    amp_broil_melt[i,"lower"] <- NA
    amp_broil_melt[i,"upper"] <- NA
  }
}

p2 <- ggplot(amp_broil_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Country),  size=0.5, width = 0) + coord_flip() +
  facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Proportion of Ampicillin Resistant Isolates") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Tetracycline in Fattening Pigs ------------------------------------------

tet_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Tet_FatPigs_Years.csv")
tet_pigs[,2:6][tet_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with

tet_pigs <- tet_pigs[!(is.na(tet_pigs$N_2015) & is.na(tet_pigs$N_2016) & is.na(tet_pigs$N_2017) & 
                         is.na(tet_pigs$N_2018) & is.na(tet_pigs$N_2019)),]

colnames(tet_pigs)[12:16] <- as.character(seq(2015, 2019))

tet_pigs_melt <- melt(tet_pigs, id.vars = "Country", measure.vars = c("2015","2016","2017","2018","2019"))

tet_pigs_melt$N <- unlist(tet_pigs[2:6])
tet_pigs_melt$PosIsol <- unlist(tet_pigs[7:11])
tet_pigs_melt$value[is.na(tet_pigs_melt$N)] <- NA
tet_pigs_melt$PosIsol[is.na(tet_pigs_melt$N)] <- NA

tet_pigs_melt$lower <- 0; tet_pigs_melt$upper <- 0

for(i in 1:nrow(tet_pigs_melt)) {
  if(!is.na(tet_pigs_melt[i, "N"])) {
    tet_pigs_melt[i,"lower"] <- prop.test(tet_pigs_melt$PosIsol[i],tet_pigs_melt$N[i])[[6]][[1]]
    tet_pigs_melt[i,"upper"] <- prop.test(tet_pigs_melt$PosIsol[i],tet_pigs_melt$N[i])[[6]][[2]]
  } else{
    tet_pigs_melt[i,"lower"] <- NA
    tet_pigs_melt[i,"upper"] <- NA
  }
}

p3 <- ggplot(tet_pigs_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Country),  size=0.5, width = 0) + coord_flip() +
  facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Proportion of Tetracycline Resistant Isolates") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Tetracycline in Broiler Poultry ------------------------------------------

tet_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Tet_Broil_Years.csv")
tet_broil[,2:5][tet_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with

tet_broil <- tet_broil[!(is.na(tet_broil$N_2014) & is.na(tet_broil$N_2016) & is.na(tet_broil$N_2017) & 
                           is.na(tet_broil$N_2018)),]

colnames(tet_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

tet_broil_melt <- melt(tet_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))

tet_broil_melt$N <- unlist(tet_broil[2:5])
tet_broil_melt$PosIsol <- unlist(tet_broil[6:9])
tet_broil_melt$value[is.na(tet_broil_melt$N)] <- NA
tet_broil_melt$PosIsol[is.na(tet_broil_melt$N)] <- NA

tet_broil_melt$lower <- 0; tet_broil_melt$upper <- 0

for(i in 1:nrow(tet_broil_melt)) {
  if(!is.na(tet_broil_melt[i, "N"])) {
    tet_broil_melt[i,"lower"] <- prop.test(tet_broil_melt$PosIsol[i],tet_broil_melt$N[i])[[6]][[1]]
    tet_broil_melt[i,"upper"] <- prop.test(tet_broil_melt$PosIsol[i],tet_broil_melt$N[i])[[6]][[2]]
  } else{
    tet_broil_melt[i,"lower"] <- NA
    tet_broil_melt[i,"upper"] <- NA
  }
}

p4 <- ggplot(tet_broil_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Country),  size=0.5, width = 0) + coord_flip() +
  facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Proportion of Tetracycline Resistant Isolates") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Saving the Plots --------------------------------------------------------

p_fatpig <- ggarrange(p1,p3, ncol = 2)

ggsave(p_fatpig, filename = "prev_fatpig.png", dpi = 300, type = "cairo", width = 8, height = 14, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures")

p_broil <- ggarrange(p2,p4, ncol = 2)

ggsave(p_broil, filename = "prev_broil.png", dpi = 300, type = "cairo", width = 8, height = 14, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures")


# Broil Usage - Stat Testing -------------------------------------------------------------

#Amp Data - Broil
amp_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Amp_Broil_Years.csv")
amp_broil[,(2+4):(5+4)][amp_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
amp_broil[,2:5][amp_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
amp_broil <- amp_broil[!(is.na(amp_broil$N_2014) & is.na(amp_broil$N_2016) & is.na(amp_broil$N_2017) & 
                           is.na(amp_broil$N_2018)),]
colnames(amp_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

#Tet Data - Broil
tet_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Tet_Broil_Years.csv")
tet_broil[,(2+4):(5+4)][tet_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
tet_broil[,2:5][tet_broil[,2:5] < 10] <- NA # replace anything under a sample size of 10 with
tet_broil <- tet_broil[!(is.na(tet_broil$N_2014) & is.na(tet_broil$N_2016) & is.na(tet_broil$N_2017) & 
                           is.na(tet_broil$N_2018)),]
colnames(tet_broil)[10:13] <- as.character(c(2014, 2016, 2017, 2018))

#Analysis 
amp_broil_melt <- melt(amp_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))

amp_broil_melt$N <- unlist(amp_broil[2:5])
amp_broil_melt$PosIsol <- unlist(amp_broil[6:9])
amp_broil_melt$value[is.na(amp_broil_melt$N)] <- NA
amp_broil_melt$PosIsol[is.na(amp_broil_melt$N)] <- NA

#Ampicillin
amp_usage_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/broil_usage_years.csv")
amp_usage_broil <- amp_usage_broil[amp_usage_broil$Country %in% intersect(amp_usage_broil$Country, amp_broil$Country),]

amp_usage_broil_test <- data.frame("usage_amp" = rowMeans(amp_usage_broil[,22:25], na.rm = TRUE),
                                   "usage_tet" = rowMeans(amp_usage_broil[,26:29], na.rm = TRUE),
                                   "N" = rowSums(amp_broil[,2:5], na.rm = TRUE),
                                   "isolpos_amp" = rowSums(amp_broil[,6:9], na.rm = TRUE),
                                   "isolpos_tet" = rowSums(tet_broil[,6:9], na.rm = TRUE))
amp_usage_broil_test$propres_amp <- amp_usage_broil_test$isolpos_amp /  amp_usage_broil_test$N
amp_usage_broil_test$propres_tet <- amp_usage_broil_test$isolpos_tet /  amp_usage_broil_test$N

amp_usage_broil_test$propres_amp_trans <- ((amp_usage_broil_test$propres_amp*(nrow(amp_usage_broil_test)-1)) + 0.5)/nrow(amp_usage_broil_test)
amp_usage_broil_test$propres_tet_trans <- ((amp_usage_broil_test$propres_tet*(nrow(amp_usage_broil_test)-1)) + 0.5)/nrow(amp_usage_broil_test)

#Stat Test
t1 <- lm(propres_amp ~ usage_amp, amp_usage_broil_test)
t1 <- betareg(propres_amp_trans ~ usage_amp, data = amp_usage_broil_test)
summary(t1)

t2 <- lm(propres_tet ~ usage_tet, amp_usage_broil_test)
t2 <- betareg(propres_tet_trans ~ usage_tet, data = amp_usage_broil_test)

plot(amp_usage_broil_test$usage_amp, amp_usage_broil_test$propres_amp)
plot(amp_usage_broil_test$usage_tet, amp_usage_broil_test$propres_tet)
summary(t2)

# Pig Usage - Stat Testing -------------------------------------------------------------

#Amp Data - Broil
amp_pig <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Amp_FatPigs_Years.csv")
amp_pig[,(2+5):(6+5)][amp_pig[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pig[,2:6][amp_pig[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pig <- amp_pig[!(is.na(amp_pig$N_2015) & is.na(amp_pig$N_2016) & is.na(amp_pig$N_2017) & 
                           is.na(amp_pig$N_2018) & is.na(amp_pig$N_2019)),]
colnames(amp_pig)[12:16] <- as.character(c(2015, 2016, 2017, 2018, 2019))

#Tet Data - Broil
tet_pig <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/Tet_FatPigs_Years.csv")
tet_pig[,(2+5):(6+5)][tet_pig[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
tet_pig[,2:6][tet_pig[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
tet_pig <- tet_pig[!(is.na(tet_pig$N_2015) & is.na(tet_pig$N_2016) & is.na(tet_pig$N_2017) & 
                                   is.na(tet_pig$N_2018) & is.na(tet_pig$N_2019)),]
colnames(tet_pig)[12:16] <- as.character(c(2015, 2016, 2017, 2018, 2019))

#Analysis 

amp_pig_melt <- melt(amp_pig, id.vars = "Country", measure.vars = c("2015", "2016", "2017", "2018", "2019"))

amp_pig_melt$N <- unlist(amp_pig[2:6])
amp_pig_melt$PosIsol <- unlist(amp_pig[7:11])
amp_pig_melt$value[is.na(amp_pig_melt$N)] <- NA
amp_pig_melt$PosIsol[is.na(amp_pig_melt$N)] <- NA

#Ampicillin
amp_usage_pig <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data/pig_usage_years.csv")
amp_usage_pig <- amp_usage_pig[amp_usage_pig$Country %in% intersect(amp_usage_pig$Country, amp_pig$Country),]

amp_usage_pig_test <- data.frame("usage_amp" = rowMeans(amp_usage_pig[,27:31], na.rm = TRUE),
                                   "usage_tet" = rowMeans(amp_usage_pig[,32:36], na.rm = TRUE),
                                   "N" = rowSums(amp_pigs[,2:6], na.rm = TRUE),
                                   "isolpos_amp" = rowSums(amp_pig[,6:9], na.rm = TRUE),
                                   "isolpos_tet" = rowSums(tet_pig[,6:9], na.rm = TRUE))
amp_usage_pig_test$propres_amp <- amp_usage_pig_test$isolpos_amp /  amp_usage_pig_test$N
amp_usage_pig_test$propres_tet <- amp_usage_pig_test$isolpos_tet /  amp_usage_pig_test$N

amp_usage_pig_test <- amp_usage_pig_test[-c(11,18),] # Need to check out number 18

amp_usage_pig_test$propres_amp_trans <- ((amp_usage_pig_test$propres_amp*(nrow(amp_usage_pig_test)-1)) + 0.5)/nrow(amp_usage_pig_test)
amp_usage_pig_test$propres_tet_trans <- ((amp_usage_pig_test$propres_tet*(nrow(amp_usage_pig_test)-1)) + 0.5)/nrow(amp_usage_pig_test)

#Stat Test
t1_pigs <- lm(propres_amp ~ usage_amp, amp_usage_pig_test)
t1_pigs <- betareg(propres_amp_trans ~ usage_amp, data = amp_usage_pig_test)
summary(t1_pigs)

t2_pigs <- lm(propres_tet ~ usage_tet, amp_usage_broil_test)
t2_pigs <- betareg(propres_tet_trans ~ usage_tet, data = amp_usage_pig_test)

plot(amp_usage_pig_test$usage_amp, amp_usage_pig_test$propres_amp)
plot(amp_usage_pig_test$usage_tet, amp_usage_pig_test$propres_tet)
summary(t2_pigs)

