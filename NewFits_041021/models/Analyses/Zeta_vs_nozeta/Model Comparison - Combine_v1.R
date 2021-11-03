library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Models/Chapter-2/NewFits_041021/data/new/compare")

#do.call works by exceuting a function and a list of arguenments passed to it
#lapply works by returning a list of all of the executed functions given a bunch of arguements passed to it


comp_list <- lapply(1:4, function(x) list.files(pattern = paste0("COMPARE_ABC_post_",c("ampbroil", "tetbroil", "amppigs", "tetpigs")[x])))
model_comp <-  lapply(1:4, function(x) do.call(rbind, lapply(comp_list[[x]], read.csv)))

for(i in 1:4) {
  list_vec <- 1:(nrow(model_comp[[i]])/1000)
  model_comp[[i]]$gen <- unlist(lapply(list_vec, function(x) rep(list_vec[x], 1000)))
}

plot_list <- list()


for(i in 1:4) {
  plot_list[[i]] <- ggplot(model_comp[[i]], aes(x = factor(gen), fill = factor(d_m))) + geom_bar(position = "fill") + theme_bw() +
    scale_x_discrete(expand = c(0, 0), name = expression(paste("Generation"))) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.text=element_text(size=14),  axis.text = element_text(size=14),
          axis.title.y=element_blank(), axis.title.x = element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          legend.title=element_text(size=14), title= element_text(size= 15)) +
    scale_fill_manual(values = c("black","grey"),name = "Model", labels = c("1 (Zeta)", "2 (No Zeta)")) + 
    labs(title = paste0(c("Ampicillin ", "Tetracycline ", "Ampicillin ", "Tetracycline ")[i] ,
                        "Usage in ", c("Broiler Poultry", "Broiler Poultry", "Fattening Pigs","Fattening Pigs")[i]))
}

comb_plot <- ggarrange(plot_list[[1]],plot_list[[2]], plot_list[[3]],plot_list[[4]], common.legend = TRUE, legend = "bottom" ,
          ncol = 1, nrow = 4)

comb_plot_annotate <- annotate_figure(comb_plot,
                left = text_grob(bquote("Proportion of " ~ italic(m) ~ "particles accepted (n = 1000)"), color = "black", rot = 90, size = 14))

ggsave(comb_plot_annotate, filename = paste0("compare_plot_supp.png"), path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Figures/comb_data",
       dpi = 300, type = "cairo", width = 10, height = 10, units = "in")


# Figuring out exact values -----------------------------------------------

comb_list <- list()

for(i in 1:4) {
  model <- model_comp[[i]]
  
  data_temp <- data.frame("gen" = 1:length(unique(model$gen)),
                          "zeta" =   sapply(1:length(unique(model$gen)), function(x) length(which(model$d_m[model$gen == x] == 1))/1000),
                          "nozeta" =   1-(sapply(1:length(unique(model$gen)), function(x) length(which(model$d_m[model$gen == x] == 1))/1000)))
  comb_list <- melt(data_temp, id.vars = "gen", measure.vars = c("zeta", "nozeta"))
}
