library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/models/Analyses/Zeta_vs_nozeta/output")

#do.call works by exceuting a function and a list of arguenments passed to it
#lapply works by returning a list of all of the executed functions given a bunch of arguements passed to it

model_comp <- do.call(rbind, lapply(list.files(pattern = "comp_"), read.csv))[,2:5]

plot_list <- list()

for(i in 1:3) {
  plot_list[[i]] <- ggplot(model_comp[model_comp$scen == unique(model_comp$scen)[i],], 
                           aes(x = gen, y = value, fill = as.factor(variable))) + geom_bar(position = "fill", stat="identity") + theme_bw() +
    scale_x_continuous(breaks = seq(1,10), expand = c(0, 0), name = expression(paste("Generation"))) + 
    scale_y_continuous(expand = c(0, 0), name = bquote("Proportion of " ~ italic(m) ~ "particles accepted (n = 1000)")) +
    theme(legend.text=element_text(size=14),  axis.text = element_text(size=14),
          axis.title.y=element_blank(), axis.title.x = element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          legend.title=element_text(size=14), title= element_text(size= 15)) +
    scale_fill_manual(values = c("black","grey"),name = "Model", labels = c("1 (Zeta)", "2 (No Zeta)")) + 
    labs(title = paste0(c("Ampicillin ", "Tetracycline ", "Tetracycline ")[i] ,
                        "usage in ", c("Fattening Pigs", "Fattening Pigs", "Broiler Poultry")[i] ))
}

comb_plot <- ggarrange(plot_list[[2]],plot_list[[1]], plot_list[[3]], common.legend = TRUE, legend = "bottom" ,
          ncol = 1, nrow = 3)

comb_plot_annotate <- annotate_figure(comb_plot,
                left = text_grob(bquote("Proportion of " ~ italic(m) ~ "particles accepted (n = 1000)"), color = "black", rot = 90, size = 14))

ggsave(comb_plot_annotate, filename = paste0("compare_plot_supp.png"), path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1",
       dpi = 300, type = "cairo", width = 10, height = 10, units = "in")