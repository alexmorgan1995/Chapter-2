library(reshape2); library(ggplot2); library(ggpubr); library(betareg)
rm(list=ls())


# Functions ---------------------------------------------------------------

sum_square_diff_dist <- function(data.obs, model.obs) {
  sumsquare <- abs((data.obs - model.obs)^2)
  return(sum(sumsquare))
}

# Import and Clean Data - Fattening Pigs (Ampicillin and Tetracycline) -----------
# Ampicillin Resistance
amp_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/Amp_FatPigs_Years.csv")
amp_pigs[,(2+5):(6+5)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,(2+10):(6+10)][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs[,2:6][amp_pigs[,2:6] < 10] <- NA # replace anything under a sample size of 10 with
amp_pigs <- amp_pigs[!(is.na(amp_pigs$N_2015) & is.na(amp_pigs$N_2016) & is.na(amp_pigs$N_2017) & 
                         is.na(amp_pigs$N_2018) & is.na(amp_pigs$N_2019)),]

colnames(amp_pigs)[12:16] <- as.character(seq(2015, 2019))

#Check that we have the same countries for tetracycline and ampicillin usage
amp_pigs$Country %in% amp_pigs$Country

#Broiler Usage
usage_pigs <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/puredata/pig_usage_years.csv")
usage_pigs <- usage_pigs[usage_pigs$Country %in% intersect(usage_pigs$Country, amp_pigs$Country),]

# NON-AGGREGATED - AMP PIGS  -------------------------------------------------------------

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

#NORMAL
stat_amp_pigs <- lm(Resistance ~ usage, melt_amp_pigs); summary(stat_amp_pigs)

ss_data_amppigs_nonaggre <- data.frame("model_pred" = predict(stat_amp_pigs, melt_amp_pigs),
                                       "data" = melt_amp_pigs$Resistance)
ss_data_amppigs_nonaggre <- ss_data_amppigs_nonaggre[!is.na(ss_data_amppigs_nonaggre$model_pred),]
ss_amppigs_nonaggre <- sum_square_diff_dist(ss_data_amppigs_nonaggre$model_pred, ss_data_amppigs_nonaggre$data)

p_amp_pigs <- ggplot(melt_amp_pigs, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp), size=0.5, width = 0) + theme_bw()  + theme(legend.position =  "bottom") +
  labs(title = "Ampicillin Usage in Fattening Pigs: 2015-2019 - non-Aggregated Data", x = "Ampicillin Usage", y = "Proportion Fattening Pigs Resistant") + 
  geom_line(aes(y = predict(stat_amp_pigs, melt_amp_pigs), x = usage), 
            colour = "red", size = 1.2) +
  annotate("text",label = paste0("Sum of squares: ", round(ss_amppigs_nonaggre, digits = 6)), y = 0.95, 
           x = max(melt_amp_pigs$usage, na.rm = TRUE)*0.75, size = 5)

#Force 0 
stat_amp_pigs_zero <- lm(Resistance ~ 0 + usage, melt_amp_pigs); summary(stat_amp_pigs_zero)

ss_data_amppigs_nonaggre_zero <- data.frame("model_pred" = predict(stat_amp_pigs_zero, melt_amp_pigs),
                                       "data" = melt_amp_pigs$Resistance)
ss_data_amppigs_nonaggre_zero <- ss_data_amppigs_nonaggre_zero[!is.na(ss_data_amppigs_nonaggre_zero$model_pred),]
ss_amppigs_nonaggre_zero <- sum_square_diff_dist(ss_data_amppigs_nonaggre_zero$model_pred, ss_data_amppigs_nonaggre_zero$data)

p_amp_pigs_zero <- ggplot(melt_amp_pigs, aes(usage, Resistance, color = Country)) + geom_point() + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp), size=0.5, width = 0) + theme_bw()  + theme(legend.position =  "bottom") +
  labs(title = "Ampicillin Usage in Fattening Pigs: 2015-2019 - non-Aggregated Data", x = "Ampicillin Usage", y = "Proportion Fattening Pigs Resistant") + 
  geom_line(aes(y = predict(stat_amp_pigs_zero, melt_amp_pigs), x = usage), 
            colour = "red", size = 1.2) +
  annotate("text",label = paste0("Sum of squares: ", round(ss_amppigs_nonaggre_zero, digits = 6)), y = 0.95, 
           x = max(melt_amp_pigs$usage, na.rm = TRUE)*0.75, size = 5)

# AGGREGATED - AMP PIGS  -------------------------------------------------

pig_rel <- data.frame("Country" = usage_pigs$Country,
                      "usage_amp" = rowMeans(usage_pigs[,27:31], na.rm = TRUE),
                      "N" = rowSums(amp_pigs[,2:6], na.rm = TRUE),
                      "isolpos_amp" = rowSums(amp_pigs[,7:11], na.rm = TRUE))

pig_rel$propres_amp <- pig_rel$isolpos_amp /  pig_rel$N

pig_rel$lower_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[1]]))
pig_rel$upper_amp <- unlist(lapply(1:nrow(pig_rel), function(i) prop.test(pig_rel$isolpos_amp[i],pig_rel$N[i])[[6]][[2]]))

#NORMAL
pig_stat_amp <- lm(propres_amp ~ usage_amp, pig_rel)
summary(pig_stat_amp)

ss_data_amppigs_aggre <- data.frame("model_pred" = predict(pig_stat_amp, pig_rel),
                                            "data" = pig_rel$propres_amp)
ss_data_amppigs_aggre <- ss_data_amppigs_aggre[!is.na(ss_data_amppigs_aggre$model_pred),]
ss_amppigs_aggre <- sum_square_diff_dist(ss_data_amppigs_aggre$model_pred, ss_data_amppigs_aggre$data)

p2_pig_amp <- ggplot(pig_rel, aes(x = usage_amp, y = propres_amp, color = Country)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_amp, pig_rel), x = usage_amp), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + theme(legend.position =  "bottom") + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
  labs(title = "Ampicillin Usage in Fat Pigs: 2015-2019 - Aggregated Data", x = "Ampicillin Usage", y = "Proportion Pigs Resistant")  +
  annotate("text",label = paste0("Sum of squares: ", round(ss_amppigs_aggre, digits = 6)), y = 0.95, 
           x = max(pig_rel$usage_amp, na.rm = TRUE)*0.75, size = 5)

#FORCE 0
pig_stat_amp_zero <- lm(propres_amp ~ 0 + usage_amp, pig_rel)
summary(pig_stat_amp_zero)

ss_data_amppigs_aggre_zero <- data.frame("model_pred" = predict(pig_stat_amp_zero, pig_rel),
                                    "data" = pig_rel$propres_amp)
ss_data_amppigs_aggre_zero <- ss_data_amppigs_aggre_zero[!is.na(ss_data_amppigs_aggre_zero$model_pred),]
ss_amppigs_aggre_zero <- sum_square_diff_dist(ss_data_amppigs_aggre_zero$model_pred, ss_data_amppigs_aggre_zero$data)

p2_pig_amp_zero <- ggplot(pig_rel, aes(x = usage_amp, y = propres_amp, color = Country)) + geom_point() + theme_bw() + 
  geom_line(aes(y = predict(pig_stat_amp_zero, pig_rel), x = usage_amp), colour = "red", size = 1.2) + 
  scale_y_continuous(limits = c(0,1)) + theme(legend.position =  "bottom") + 
  geom_errorbar(aes(ymin=lower_amp, ymax=upper_amp),  size=0.5, width = 0) + 
  labs(title = "Ampicillin Usage in Fat Pigs: 2015-2019 - Aggregated Data", x = "Ampicillin Usage", y = "Proportion Pigs Resistant")  +
  annotate("text",label = paste0("Sum of squares: ", round(ss_amppigs_aggre_zero, digits = 6)), y = 0.95, 
           x = max(pig_rel$usage_amp, na.rm = TRUE)*0.75, size = 5)

# AGGREGATED Plotz -------------------------------------------------------------------

ss_plots <- ggarrange(p_amp_pigs, p_amp_pigs_zero,
                      p2_pig_amp,p2_pig_amp_zero, ncol = 2, nrow = 2)

ggsave(ss_plots, filename = "ss_zerovsnzero.png", dpi = 300, type = "cairo", width = 13, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data")
