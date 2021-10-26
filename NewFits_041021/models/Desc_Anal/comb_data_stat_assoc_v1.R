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

# Create the non-Aggregated Plots  -------------------------------------------------------------

#Create the Combined Dataset - Amp Fattening Pigs 

melt_amp_pigs <- melt(amp_pigs, id.vars = "Country", measure.vars = c("2015", "2016", "2017", "2018", "2019"))
melt_amp_pigs$usage <- melt(usage_pigs, id.vars = "Country", measure.vars = c("scale_ampusage_2015", "scale_ampusage_2016", 
                                                       "scale_ampusage_2017", "scale_ampusage_2018", "scale_ampusage_2019"))[,3]
melt_amp_pigs$N <- melt(amp_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                          "N_2017", "N_2018", "N_2019"))[,3]
melt_amp_pigs$IsolPos <- melt(amp_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                        "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_amp_pigs)[c(2,3)] <- c("Year", "Resistance")
melt_amp_pigs <- melt_amp_pigs[!is.na(melt_amp_pigs$N),]

melt_amp_pigs$lower_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[1]]))
melt_amp_pigs$upper_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[2]]))

p_amp_pigs <- ggplot(melt_amp_pigs, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp), size=0.5, width = 0) + theme_bw()  + theme(legend.position =  "bottom", plot.title = element_text(face = "bold")) +
  labs(title = "Ampicillin Usage in Fattening Pigs: 2015-2018", x = "Ampicillin Fattening Pig Sales (mg/PCU)", y = "Proportion Fattening Pigs Resistant")+
  scale_y_continuous(expand = c(0,0))  +  scale_x_continuous(expand = c(0,0.3))+
  geom_smooth(method=lm , color="red", fill="hotpink", se=TRUE) + coord_cartesian(ylim=c(0,1.1)) 

#Create the Combined Dataset - Tet Fattening Pigs 

melt_tet_pigs <- melt(tet_pigs, id.vars = "Country", measure.vars = c("2015", "2016", "2017", "2018", "2019"))
melt_tet_pigs$usage <- melt(usage_pigs, id.vars = "Country", measure.vars = c("scale_tetusage_2015", "scale_tetusage_2016", 
                                                                              "scale_tetusage_2017", "scale_tetusage_2018", "scale_tetusage_2019"))[,3]
melt_tet_pigs$N <- melt(tet_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                        "N_2017", "N_2018", "N_2019"))[,3]
melt_tet_pigs$IsolPos <- melt(tet_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                              "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_tet_pigs)[c(2,3)] <- c("Year", "Resistance")
melt_tet_pigs <- melt_tet_pigs[!is.na(melt_tet_pigs$N),]

melt_tet_pigs$lower_tet <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[1]]))
melt_tet_pigs$upper_tet <- unlist(lapply(1:nrow(melt_tet_pigs), function(i) prop.test(melt_tet_pigs$IsolPos[i],melt_tet_pigs$N[i])[[6]][[2]]))

p_tet_pigs <- ggplot(melt_tet_pigs, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet), size=0.5, width = 0) + theme_bw()  + theme(legend.position =  "bottom", plot.title = element_text(face = "bold")) +
  labs(title = "Tetracycline Usage in Fattening Pigs: 2015-2018", x = "Tetracycline Fattening Pig Sales (mg/PCU)", y = "Proportion Fattening Pigs Resistant")+
  scale_y_continuous(expand = c(0,0))  +  scale_x_continuous(expand = c(0,0.3))+ 
  geom_smooth(method=lm , color="red", fill="hotpink", se=TRUE) + coord_cartesian(ylim=c(0,1.1)) 
 
  

#Create the Combined Dataset - Amp Broilers 

melt_amp_broil <- melt(amp_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))
melt_amp_broil$usage <- melt(usage_broil, id.vars = "Country", measure.vars = c("scale_ampusage_2014", "scale_ampusage_2016", 
                                                                              "scale_ampusage_2017", "scale_ampusage_2018"))[,3]
melt_amp_broil$N <- melt(amp_broil, id.vars = "Country", measure.vars = c("N_2014", "N_2016", "N_2017", "N_2018"))[,3]
melt_amp_broil$IsolPos <- melt(amp_broil, id.vars = "Country", measure.vars = c("PosIsol_2014", "PosIsol_2016", 
                                                                              "PosIsol_2017", "PosIsol_2018"))[,3]
colnames(melt_amp_broil)[c(2,3)] <- c("Year", "Resistance")
melt_amp_broil <- melt_amp_broil[!is.na(melt_amp_broil$N),]

melt_amp_broil$lower_amp <- unlist(lapply(1:nrow(melt_amp_broil), function(i) prop.test(melt_amp_broil$IsolPos[i],melt_amp_broil$N[i])[[6]][[1]]))
melt_amp_broil$upper_amp <- unlist(lapply(1:nrow(melt_amp_broil), function(i) prop.test(melt_amp_broil$IsolPos[i],melt_amp_broil$N[i])[[6]][[2]]))

p_amp_broil <- ggplot(melt_amp_broil, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp), size=0.5, width = 0) + theme_bw() + theme(legend.position =  "bottom", plot.title = element_text(face = "bold")) +
  labs(title = "Ampicillin Usage in Broilers: 2014-2018", x = "Ampicillin Broiler Sales (mg/PCU)", y = "Proportion Broilers Resistant")+
  scale_y_continuous(expand = c(0,0))  +  scale_x_continuous(expand = c(0,0.3))+
  geom_smooth(method=lm , color="red", fill="hotpink", se=TRUE)  + coord_cartesian(ylim=c(0,1.1)) 

#Create the Combined Dataset - Tet Broilers  

melt_tet_broil <- melt(tet_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))
melt_tet_broil$usage <- melt(usage_broil, id.vars = "Country", measure.vars = c("scale_tetusage_2014", "scale_tetusage_2016", 
                                                                                "scale_tetusage_2017", "scale_tetusage_2018"))[,3]
melt_tet_broil$N <- melt(tet_broil, id.vars = "Country", measure.vars = c("N_2014", "N_2016", "N_2017", "N_2018"))[,3]
melt_tet_broil$IsolPos <- melt(tet_broil, id.vars = "Country", measure.vars = c("PosIsol_2014", "PosIsol_2016", 
                                                                                "PosIsol_2017", "PosIsol_2018"))[,3]
colnames(melt_tet_broil)[c(2,3)] <- c("Year", "Resistance")
melt_tet_broil <- melt_tet_broil[!is.na(melt_tet_broil$N),]

melt_tet_broil$lower_tet <- unlist(lapply(1:nrow(melt_tet_broil), function(i) prop.test(melt_tet_broil$IsolPos[i],melt_tet_broil$N[i])[[6]][[1]]))
melt_tet_broil$upper_tet <- unlist(lapply(1:nrow(melt_tet_broil), function(i) prop.test(melt_tet_broil$IsolPos[i],melt_tet_broil$N[i])[[6]][[2]]))

p_tet_broil <- ggplot(melt_tet_broil, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet), size=0.5, width = 0) + theme_bw()  + theme(legend.position =  "bottom", plot.title = element_text(face = "bold")) +
  labs(title = "Tetracycline Usage in Broilers: 2014-2018", x = "Tetracycline Broiler Sales (mg/PCU)", y = "Proportion Broilers Resistant") +
  scale_y_continuous(expand = c(0,0)) +  scale_x_continuous(expand = c(0,0.3)) +
  geom_smooth(method=lm , color="red", fill="hotpink", se=TRUE)

# Statistically Testing the non-agrgegated Data  -----------------------------------------

stat_amp_pigs <- lm(Resistance ~ usage, melt_amp_pigs); summary(stat_amp_pigs)
stat_tet_pigs <- lm(Resistance ~ usage, melt_tet_pigs); summary(stat_tet_pigs)

stat_amp_broil <- lm(Resistance ~ usage, melt_amp_broil); summary(stat_amp_broil)
stat_tet_broil <- lm(Resistance ~ usage, melt_tet_broil); summary(stat_tet_broil)

comb_nonaggre_plot <- ggarrange(p_amp_broil,p_tet_broil, p_amp_pigs, p_tet_pigs, labels = c("A","B","C","D"),ncol = 2,  nrow = 2, common.legend = TRUE,
                                legend = "bottom")

ggsave(comb_nonaggre_plot, filename = "nonaggreg_stat.png", dpi = 300, type = "cairo", width = 11, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

# Create aggregated data  -------------------------------------------------

broil_rel <- data.frame("country" = usage_broil$Country,
                        "usage_amp" = rowMeans(usage_broil[,22:25], na.rm = TRUE),
                        "usage_tet" = rowMeans(usage_broil[,26:29], na.rm = TRUE),
                        "N" = rowSums(amp_broil[,2:5], na.rm = TRUE),
                        "isolpos_amp" = rowSums(amp_broil[,6:9], na.rm = TRUE),
                        "isolpos_tet" = rowSums(tet_broil[,6:9], na.rm = TRUE))

broil_rel$propres_amp <- broil_rel$isolpos_amp /  broil_rel$N
broil_rel$propres_tet <- broil_rel$isolpos_tet /  broil_rel$N

broil_rel$lower_amp <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_amp[i],broil_rel$N[i])[[6]][[1]]))
broil_rel$upper_amp <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_amp[i],broil_rel$N[i])[[6]][[2]]))
broil_rel$lower_tet <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_tet[i],broil_rel$N[i])[[6]][[1]]))
broil_rel$upper_tet <- unlist(lapply(1:nrow(broil_rel), function(i) prop.test(broil_rel$isolpos_tet[i],broil_rel$N[i])[[6]][[2]]))

#Stat Test
broil_stat_amp <- lm(propres_amp ~ usage_amp, broil_rel)
broil_stat_tet <- lm(propres_tet ~ usage_tet, data = broil_rel)
summary(broil_stat_amp)
summary(broil_stat_tet)

p2_broil_amp <- ggplot(broil_rel, aes(x = usage_amp, y = propres_amp, color = country)) + geom_point() + theme_bw() +  
  geom_line(aes(y = predict(broil_stat_amp, broil_rel), x = usage_amp), colour = "red", size = 1.2) + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + theme(legend.position =  "bottom") +
  labs(title = "Ampicillin Usage in Broilers: 2014-2018", x = "Ampicillin Usage", y = "Proportion Broilers Resistant") +
  annotate("text",label = paste0("p-value: ", round(summary(broil_stat_amp)[[4]][[8]], digits = 6)), y = 0.95, 
           x = max(broil_rel$usage_amp, na.rm = TRUE)*0.75, size = 5)

p2_broil_tet <- ggplot(broil_rel, aes(x = usage_tet, y = propres_tet, color = country)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(broil_stat_tet, broil_rel), x = usage_tet), colour = "red", size = 1.2) + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet),  size=0.5, width = 0)+ theme(legend.position =  "bottom") + 
  labs(title = "Tetracycline Usage in Broilers: 2014-2018", x = "Tetracycline Usage", y = "Proportion Broilers Resistant") +
  annotate("text",label = paste0("p-value: ", round(summary(broil_stat_tet)[[4]][[8]], digits = 6)), y = 0.95, 
           x = max(broil_rel$usage_tet, na.rm = TRUE)*0.75, size = 5)

# Stat Testing - FATTENING PIGS -------------------------------------------------------------

pig_rel <- data.frame("Country" = usage_pigs$Country,
                      "usage_amp" = rowMeans(usage_pigs[,27:31], na.rm = TRUE),
                      "usage_tet" = rowMeans(usage_pigs[,32:36], na.rm = TRUE),
                      "N" = rowSums(amp_pigs[,2:6], na.rm = TRUE),
                      "isolpos_amp" = rowSums(amp_pigs[,7:11], na.rm = TRUE),
                      "isolpos_tet" = rowSums(tet_pigs[,7:11], na.rm = TRUE))

pig_rel$propres_amp <- pig_rel$isolpos_amp /  pig_rel$N
pig_rel$propres_tet <- pig_rel$isolpos_tet /  pig_rel$N

pig_rel$lower_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[1]]))
pig_rel$upper_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[2]]))
pig_rel$lower_tet <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_tet[i],pig_rel$N[i])[[6]][[1]]))
pig_rel$upper_tet <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_tet[i],pig_rel$N[i])[[6]][[2]]))

#Stat Test
pig_stat_amp <- lm(propres_amp ~ usage_amp, pig_rel)
pig_stat_tet <- lm(propres_tet ~ usage_tet, data = pig_rel)
summary(pig_stat_amp)
summary(pig_stat_tet)

p2_pig_amp <- ggplot(pig_rel, aes(x = usage_amp, y = propres_amp, color = Country)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_amp, pig_rel), x = usage_amp), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + theme(legend.position =  "bottom") + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
  labs(title = "Ampicillin Usage in Fat Pigs: 2015-2019", x = "Ampicillin Usage", y = "Proportion Pigs Resistant")  +
  annotate("text",label = paste0("p-value: ", round(summary(pig_stat_amp)[[4]][[8]], digits = 6)), y = 0.95, 
           x = max(pig_rel$usage_tet, na.rm = TRUE)*0.75, size = 5)

p2_pig_tet <- ggplot(pig_rel, aes(x = usage_tet, y = propres_tet, color = Country)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_tet, pig_rel), x = usage_tet), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + theme(legend.position =  "bottom") + 
  geom_errorbar(aes(ymin=lower_tet, ymax=upper_tet),  size=0.5, width = 0) + 
  labs(title = "Tetracycline Usage in Fat Pigs: 2015-2019", x = "Tetracycline Usage", y = "Proportion Pigs Resistant")  +
  annotate("text",label = paste0("p-value: ", round(summary(pig_stat_tet)[[4]][[8]], digits = 6)), y = 0.95, 
           x = max(pig_rel$usage_tet, na.rm = TRUE)*0.75, size = 5)

comb_aggre_plot <- ggarrange(p2_pig_tet, p2_pig_amp, p2_broil_tet, p2_broil_amp, ncol = 2,  nrow = 2)

ggsave(comb_aggre_plot, filename = "aggreg_stat.png", dpi = 300, type = "cairo", width = 14, height = 9, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")
