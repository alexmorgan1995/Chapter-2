library(reshape2)
library(ggplot2)

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

ggplot(amp_pigs_melt, aes(x = variable, y = value, group = Country, color = Country)) + geom_point(position = position_dodge(2)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color = Country),  size=0.5, width = 0) + coord_flip() +
  facet_grid(Country~., scales = "free", space = "free") + theme_bw() + labs(x = "Surveillance Year", y = "Proportion of Tetracycline Resistant Isolates") +
  theme(strip.text.y = element_text(angle = 0), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

