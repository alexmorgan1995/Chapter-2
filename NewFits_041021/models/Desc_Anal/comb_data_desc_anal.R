library(reshape2); library(ggplot2); library(ggpubr); library(betareg)
rm(list=ls())

# Ampicillin in Fattening Pigs ------------------------------------------

amp_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_FatPigs_Years.csv")
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

amp_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_Broil_Years.csv")
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

tet_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Tet_FatPigs_Years.csv")
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

tet_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Tet_Broil_Years.csv")
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

ggsave(p_fatpig, filename = "prev_fatpig.png", dpi = 300, type = "cairo", width = 12, height = 12, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

p_broil <- ggarrange(p2,p4, ncol = 2)

ggsave(p_broil, filename = "prev_broil.png", dpi = 300, type = "cairo", width = 12, height = 12, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")


# USAGE - Fattening Pigs ------------------------------------------------------

#Ampicillin
amp_use_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/pig_usage_years.csv")
amp_use_pigs <- amp_use_pigs[,c(1, 27:31 )]
colnames(amp_use_pigs)[2:6] <- as.character(seq(2015, 2019))

amp_pigsusage_melt <- melt(amp_use_pigs, id.vars = "Country", measure.vars = c("2015","2016","2017","2018","2019"))

p_amp_pigsusage <- ggplot(amp_pigsusage_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  coord_flip() + facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Ampicillin Usage in Fat Pigs (mg/PCU)") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),axis.text.y=element_blank() ,
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#Tetracycline
tet_use_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/pig_usage_years.csv")
tet_use_pigs <- tet_use_pigs[,c(1, 32:36 )]
colnames(tet_use_pigs)[2:6] <- as.character(seq(2015, 2019))

tet_pigsusage_melt <- melt(tet_use_pigs, id.vars = "Country", measure.vars = c("2015","2016","2017","2018","2019"))

p_tet_pigsusage <- ggplot(tet_pigsusage_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  coord_flip() + facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Tetracycline Usage in Fat Pigs (mg/PCU)") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),axis.text.y=element_blank() ,
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# USAGE - Broiler Poultry  ------------------------------------------------------

#Ampicillin
amp_use_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/broil_usage_years.csv")
amp_use_broil <- amp_use_broil[,c(1, 22:25)]
colnames(amp_use_broil)[2:5] <- as.character(c(2014, 2016, 2017, 2018))

amp_broilusage_melt <- melt(amp_use_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))

p_amp_broilusage <- ggplot(amp_broilusage_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  coord_flip() + facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Ampicillin Usage in Broilers (mg/PCU)") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),axis.text.y=element_blank() ,
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#Tetracycline
tet_use_broil <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/broil_usage_years.csv")
tet_use_broil <- tet_use_broil[,c(1, 26:29)]
colnames(tet_use_broil)[2:5] <- as.character(c(2014, 2016, 2017, 2018))

tet_broilusage_melt <- melt(tet_use_broil, id.vars = "Country", measure.vars = c("2014", "2016", "2017", "2018"))

p_tet_broilusage <- ggplot(tet_broilusage_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  coord_flip() + facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Tetracycline Usage in Broilers (mg/PCU)") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),axis.text.y=element_blank() ,
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Saving the Plots - Usage --------------------------------------------------------

p_fatpigusage <- ggarrange(p_amp_pigsusage, p_tet_pigsusage, ncol = 2)

ggsave(p_fatpigusage, filename = "usage_fatpig.png", dpi = 300, type = "cairo", width = 8, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")

p_broilusage <- ggarrange(p_amp_broilusage, p_tet_broilusage, ncol = 2)

ggsave(p_broilusage, filename = "usage_broil.png", dpi = 300, type = "cairo", width = 8, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")
