library(ggtree)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(reshape2)
library(dplyr)

# read tree file into variable tr
tr_nwk <- read.tree("CC10_coreSNP_newick.tree")
tr <- read.tree("CC10_coreSNP.tree")

#read data into variable data
tree_meta <- read.csv('Table.S1.csv',stringsAsFactors = F)
rownames(tree_meta) <- tree_meta[,1]
tree_meta <- tree_meta[,1:12]
tree_meta$ST <- as.factor(tree_meta$ST)

#Create ggtree objects from tr_nwk and tr and bind metadata to tree
tree_nwk <- ggtree(tr_nwk, layout = "circular") %<+% tree_meta
tree <- ggtree(tr)

#Extract node confidence values from tree
tree$data[nrow(tree$data)+1,] <- NA
d <- tree$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
colnames(d)[6] <- "Node Confidence"

#Extract data from Newick tree
#tree_nwk$data[nrow(tree_nwk$data)-1,] <- NA
e <- tree_nwk$data
e <- e[!e$isTip,]
e$label <- as.numeric(e$label)

#Bind node values from Newick data to tree data
de <- left_join(x = d[,c(1,6)], y = e, by = "node")

#Define colours and break for source tree
full_source_breaks <- sort(unique(factor(tree_meta$Source)))
full_source_colours <- c("#01c199",
                         "#7f0154",
                         "#cc8700",
                         "#0041ab",
                         "#e974e7",
                         "#0b0104",
                         "#900caf",
                         "#0e570c",
                         "#01A9DB",
                         "#f6e28e",
                         "#ee3c25",
                         "#933f35",
                         "#00ffff")
#Source tree
tree_nwk + 
  geom_point(data = de, aes(fill = `Node Confidence`), shape =21, size =0.9) +
  geom_treescale(x =, y=, fontsize = 4, offset = 2, width =) +
  geom_tiplab2(aes(color = Source), 
               size = 2.8, 
               align = T,
               linetype = 2, 
               linesize = 0.01) + 
  scale_color_manual(breaks = full_source_breaks, 
                     values = full_source_colours)+
  scale_fill_continuous(low = "black", high = "white") +
  theme(legend.position = "left") 

#ST tree
tree_nwk + 
  geom_point(data = de, aes(fill = `Node Confidence`), shape =21, size =0.9) +
  geom_treescale(x =, y=, fontsize = 4, offset = 2, width =) +
  geom_tiplab2(aes(color = ST), 
               size = 2.8, 
               align = T,
               linetype = 2, 
               linesize = 0.01) + 
  scale_fill_continuous(low = "black", high = "white") +
  theme(legend.position = "left") 

#Define breaks and colours for continent tree
cont_breaks <- sort(unique(factor(tree_meta$Continent)))
cont_vals <- c("#01c199","#9a006a","#cc8700","#0041ab","#af4ccc")

#Continent tree
tree_nwk + 
geom_point(data = de, aes(fill = `Node Confidence`), shape =21, size =0.9) +
geom_treescale(x =, y=, fontsize = 4, offset = 2, width =) +
geom_tiplab2(aes(color = Continent), 
             size = 2.8, 
             align = T,
             linetype = 2, 
             linesize = 0.01) + 
scale_color_manual(breaks = cont_breaks, values = cont_vals)+
scale_fill_continuous(low = "black", high = "white") +
theme(legend.position = "left") 

#Define break and colours for origin tree
origin_breaks <- sort(unique(factor(tree_meta$Origin)))
origin_colours <- c("#01c199",
                          "#9a006a",
                          "#cc8700",
                          "#0041ab")
#Origin tree
tree_nwk + 
 geom_point(data = de, aes(fill = `Node Confidence`), shape =21, size =0.9) +
 geom_treescale(x =, y=, fontsize = 4, offset = 2, width =) +
 geom_tiplab2(aes(color = Origin), 
              size = 2.8, 
              align = T,
              linetype = 2, 
              linesize = 0.01) + 
 scale_color_manual(breaks = origin_breaks, values = origin_colours)+
 scale_fill_continuous(low = "black", high = "white") +
 theme(legend.position = "left") 