library(ggtree)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(reshape2)
library(dplyr)

# read tree file into variable tr
heat_tree <- read.tree("HS_core.full.filt.gub.raxml.200818.rooted.tree")

#read data into variable data and organise
heat_dat <- read.csv('Table.S1.v2.csv', stringsAsFactors = T)
rownames(heat_dat) <- heat_dat[,1]
heat_meta <- heat_dat[,1:21]
heat_dat[heat_dat > 0] <- "Present"
heat_dat[heat_dat == 0] <- "Absent"
heat_meta$ST <- as.factor(heat_meta$ST..Achtman.)

#subsets
heat_vir_cols <- c(138:256)
heat_vir_data <- heat_dat[,heat_vir_cols]
heat_res_plas_cols <- c(28:137)
heat_resplas_data <- heat_dat[,heat_res_plas_cols]

#colours for heatmap
heat_breaks <- c("Present", "Absent")
heat_colours <- c("white","#63c0af")


#tree
p <- ggtree(heat_tree, 
            aes(color=ST),
            layout = "rectangular", 
            branch.length = "none") %<+% heat_meta +
  theme(legend.position = "right") +
  geom_tiplab(size = .85, align = T)

#virulence heatmap
gheatmap(p, 
         data = heat_vir_data, 
         offset = 2.8,
         colnames_position = "top", 
         font.size =1, 
         width = 4,
         colnames_angle = 45,
         hjust = 0,
         colnames_offset_x = -.2,
         colnames_offset_y = 0.7,
         color = "grey") +
  scale_fill_manual(breaks = heat_breaks,
                    values = heat_colours) +
  theme(legend.position = "right")

#resistance heatmap
gheatmap(p, 
         data = heat_resplas_data, 
         colnames_position = "top",
         offset =2.8,
         font.size =1, 
         width =4,
         colnames_angle = 45,
         hjust = 0,
         colnames_offset_x = -.2,
         colnames_offset_y = 0.7,
         color = "grey") +
  scale_fill_manual(breaks = heat_breaks,
                    values = heat_colours) +
  theme(legend.position = "right")

