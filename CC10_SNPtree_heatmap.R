library(ggtree)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(reshape2)
library(dplyr)

# read tree file into variable tr
heat_tree <- read.tree("CC10_coreSNP_newick.tree")

#read data into variable data and organise
heat_dat <- read.csv('Table.S1.csv', stringsAsFactors = T)
rownames(heat_dat) <- heat_dat[,1]
heat_meta <- heat_dat[,1:12]
heat_dat[heat_dat > 0] <- "Present"
heat_dat[heat_dat == 0] <- "Absent"
heat_meta$ST <- as.factor(heat_meta$ST)

#subsets
heat_vir_cols <- c(129:247)
heat_vir_data <- heat_dat[,heat_vir_cols]
heat_res_plas_cols <- c(19:128)
heat_resplas_data <- heat_dat[,heat_res_plas_cols]

#colours for heatmap
heat_breaks <- c("Present", "Absent")
heat_colours <- c("white","#63c0af")

#tree
p <- ggtree(heat_tree, 
            layout = "rectangular", 
            branch.length = "none", 
            aes(color = ST)) %<+% heat_meta +
  theme(legend.position = "right")

#virulence heatmap
gheatmap(p, 
         data = heat_vir_data, 
         colnames_position = "top", 
         font.size =1.5, 
         width =6,
         colnames_angle = 45,
         hjust = 0,
         colnames_offset_x = -.5,
         colnames_offset_y = 0.8,
         color = "grey") +
  scale_fill_manual(breaks = heat_breaks,
                    values = heat_colours) +
  theme(legend.position = "right")

#resistance heatmap
gheatmap(p, 
         data = heat_resplas_data, 
         colnames_position = "top", 
         font.size =1.5, 
         width =6,
         colnames_angle = 45,
         hjust = 0,
         colnames_offset_x = -.4,
         colnames_offset_y = 0.8,
         color = "grey") +
  scale_fill_manual(breaks = heat_breaks,
                    values = heat_colours) +
  theme(legend.position = "right")

