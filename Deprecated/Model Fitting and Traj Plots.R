library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("sensitivity"); library("fast")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData")

# Model Functions ----------------------------------------------------------

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model ODEs
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa  
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

#Importing in the Datasets
import <- function(id) {
  data <- data.frame(matrix(ncol = 6, nrow = 0))
  for(i in 1:5) {
    test  <- cbind(read.csv(paste0("results_ABC_SMC_gen_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

# Data Import -------------------------------------------------------------
#Import of Posterior Distributions
data <- do.call(rbind, list(import(tet), import(amp), import(broil)))

#Import of Fitting Data
dataamp <- read.csv("resistanceprofAnim_amp.csv")
dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))
dataamp <- dataamp[!dataamp$N < 10,]

datatetra <- read.csv("resistanceprofAnim_v1.csv")
datatetra$lower <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[1]]))
datatetra$upper <- unlist(lapply(1:nrow(datatetra), function(i) prop.test(datatetra$Positive.Sample[i],datatetra$N[i])[[6]][[2]]))
datatetra <- datatetra[!datatetra$N < 10,]

databroil <- read.csv("salm_broilers_2018.csv")
databroil$lower <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[1]]))
databroil$upper <- unlist(lapply(1:nrow(databroil), function(i) prop.test(databroil$Positive.Sample[i],databroil$N[i])[[6]][[2]]))
databroil <- databroil[!databroil$N < 10,]

#Obtain the MAPs for each dataset

MAPtet <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "tet")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "tet")]),
            "betaAA" <- map_estimate(data$betaAA[which(data$group == "data5" & data$fit == "tet")]),
            "alpha" <- map_estimate(data$alpha[which(data$group == "data5" & data$fit == "tet")]))

MAPamp <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "amp")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "amp")]),
            "betaAA" <- map_estimate(data$betaAA[which(data$group == "data5" & data$fit == "amp")]),
            "alpha" <- map_estimate(data$alpha[which(data$group == "data5" & data$fit == "amp")]))

MAPbroil <- c("phi" <- map_estimate(data$phi[which(data$group == "data5" & data$fit == "broil")]),
            "theta" <- map_estimate(data$theta[which(data$group == "data5" & data$fit == "broil")]),
            "betaAA" <- map_estimate(data$betaAA[which(data$group == "data5" & data$fit == "broil")]),
            "alpha" <- map_estimate(data$alpha[which(data$group == "data5" & data$fit == "broil")]))

MAP <- rbind(MAPtet, MAPamp, MAPbroil); colnames(MAP) <- c("phi", "theta", "betaAA", "alpha")

# ABC-SMC Posterior -------------------------------------------------------

plotlist1 <- list()

for(j in 1:length(unique(data$fit))) {
  plotlist2 <- list()
  
  for (i in 1:length(colnames(MAP))) { # Loop over loop.vector
    
    dens <- density(data[,colnames(MAP)[i]][data$fit == unique(data$fit)[j] & data$group == "data5"])
    dataplot <- melt(data[data$fit == unique(data$fit)[j],], id.vars = "group", measure.vars = colnames(MAP)[i])
    
    plotlist2[[i]] <- local({
      i = i
      p1 <- ggplot(data[data$fit == unique(data$fit)[j],], aes(x=get(colnames(MAP)[i]), fill=group)) + geom_density(alpha=.5) + theme_bw()  +
        scale_fill_discrete(labels = c("Generation 1", "Generation 2", "Generation 3", "Generation 4", "Generation 5"))+
        theme(legend.text=element_text(size=10), axis.text.x=element_text(size=10),axis.ticks.y=element_blank(), axis.text.y=element_blank(),
              axis.title.y=element_text(size=10), axis.title.x= element_text(size=10), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
              plot.title = element_text(size = 12, vjust = 3, hjust = 0.5, face = "bold"))
      if(colnames(MAP)[i] == "phi") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.04), expand = c(0, 0), name = expression(paste("Rate of Resistance Reversion (", phi, ")"))) +
        labs(fill = NULL, title = c("Tetracycline Sales in Fattening Pigs", "Ampicillin Sales in Fattening Pigs", "Tetracycline Sales in Boiler Poultry")[j]) + 
          scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "theta") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.5),expand = c(0, 0), name = expression(paste("Efficacy of Antibiotic-Mediated Recovery (", theta, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, 5), expand = c(0, 0), name = " ") 
      }
      if(colnames(MAP)[i] == "betaAA") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.2),expand = c(0, 0), name = expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")))+
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      if(colnames(MAP)[i] == "alpha") {
        p1 <- p1 + scale_x_continuous(limits = c(0,0.6),expand = c(0, 0), name = expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")"))) +
          labs(fill = NULL, title = "") + scale_y_continuous(limits = c(0, max(dens$y)*1.2), expand = c(0, 0), name = " ")
      }
      return(p1)
      
    })
    print(paste0("Plot Parameter: ",colnames(MAP)[i], " | Data: ", unique(data$fit)[j] ))
  }
  plotlist1[[unique(data$fit)[j]]] <- plotlist2
}

abc <- ggarrange(plotlist1[[1]][[1]], plotlist1[[2]][[1]], plotlist1[[3]][[1]],
                 plotlist1[[1]][[2]], plotlist1[[2]][[2]], plotlist1[[3]][[2]],
                 plotlist1[[1]][[3]], plotlist1[[2]][[3]], plotlist1[[3]][[3]],
                 plotlist1[[1]][[4]], plotlist1[[2]][[4]], plotlist1[[3]][[4]],
                 nrow = 4, ncol =3, 
                 labels =  c("A","B","C",
                             "","","",
                             "","","",
                             "","",""),
                 font.label = c(size = 20), common.legend = TRUE, legend = "bottom",
                 align = "hv", vjust = 1.05)

ggsave(abc, filename = "ABC_SMC_Post.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft Figures")

# Model Fit with Data -----------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:nrow(MAP)) {
  output1 <- data.frame()
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(NA, nrow = 1, ncol=7))
    
    if(rownames(MAP)[j] == "MAPbroil") {
      parms2 = c(ra = 0, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i])
    } 
    
    else {
      parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP[j,"betaAA"], betaAH = 0.00001, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = MAP[j,"phi"], theta = MAP[j,"theta"], alpha = MAP[j,"alpha"], tau = parmtau[i])
    }

    out <- ode(y = init, func = amr, times = times, parms = parms2)
    temp[1,1] <- parmtau[i]
    temp[1,2] <- rounding(out[nrow(out),5]) 
    temp[1,3] <- rounding(out[nrow(out),6]) 
    temp[1,4] <- rounding(out[nrow(out),7])
    temp[1,5] <- temp[1,3] + temp[1,4]
    temp[1,6] <- signif(as.numeric(temp[1,4]/temp[1,5]), digits = 3)
    temp[1,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[1,8] <- rownames(MAP)[j]
    print(temp[1,3])
    output1 <- rbind.data.frame(output1, temp)
  }
  icombhdata <- rbind(icombhdata, output1)
}

colnames(icombhdata)[1:8] <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")
icombhdata[,2:5] <- icombhdata[,2:5]*100000 #Scaling the prevalence (per 100,000)

amp <- ggplot(dataamp, aes(x = pig_amp_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.035),expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Ampicillin Sales (g/PCU)", y = "Ampicillin-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPtet",], aes(x = tau, y= IResRatA), col = "red", size = 1.02) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

tet <- ggplot(datatetra, aes(x = pig_tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.035), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Fattening Pig Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Fattening Pig Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPamp",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))

broil <- ggplot(databroil, aes(x = tetra_sales/1000, y= ResPropAnim))  + geom_point() + theme_bw() + 
  scale_x_continuous(limits = c(0,0.035), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Broiler Poultry Tetracycline Sales (g/PCU)", y = "Tetracycline-Resistant Broiler Poultry Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", size=1.01, inherit.aes =  TRUE) + 
  geom_line(data = icombhdata[icombhdata$group == "MAPbroil",], aes(x = tau, y= IResRatA), col = "red", size = 1.02)+
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))


fits <- ggarrange(tet, amp, broil, nrow = 1, ncol = 3, align = "h", labels = c("A","B", "C"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(fits, filename = "Model_Fits.png", dpi = 300, type = "cairo", width = 14, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft Figures")

# Tau and ICombH Base Plot ----------------------------------------------------------

icombhlist <- list()

#tet, amp, broil
averagesales <- c(0.0106, 0.009877583, 0.0062)

for(i in 1: length(unique(icombhdata$group))){
  i = i 
  
  icombhlist[[i]] <- local({
    
    plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[i],],
                     id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 
    
    p1 <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
      geom_vline(xintercept = averagesales[i], alpha = 0.3, size = 2) + 
      geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
      scale_y_continuous(limits = c(0,5), expand = c(0, 0))  + 
      geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[i]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
                position = "stack", angle = 45) +
      theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm")) + 
      scale_fill_manual(labels = c("Antibiotic-Sensitive Infection", "Antibiotic-Resistant Infection"), values = c("#F8766D", "#619CFF")) 
    
    if(unique(icombhdata$group)[i] == "MAPtet") {
      p1 <- p1 + labs(x ="Tetracycline Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)")  
    }
    if(unique(icombhdata$group)[i] == "MAPamp") {
      p1 <- p1 + labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)")  
    }
    if(unique(icombhdata$group)[i] == "MAPbroil") {
      p1 <- p1 + labs(x ="Tetracycline Sales in Broiler Poultry (g/PCU)", y = "Infected Humans (per 100,000)")  
    }
    
    return(p1)
    
  })
  
}

icombh <- ggarrange(icombhlist[[1]], icombhlist[[2]], icombhlist[[3]], nrow = 3, ncol = 1, align = "v", labels = c("A","B", "C"), font.label = c(size = 20),
                  common.legend = TRUE, legend = "bottom") 

ggsave(icombh, filename = "Icombh.png", dpi = 300, type = "cairo", width = 8, height = 11, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft Figures")

