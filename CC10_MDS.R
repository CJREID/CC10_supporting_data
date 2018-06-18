library(ggplot2)
library(MASS)
library(vegan)
library(dplyr)
library(magrittr)
library(ellipse)
library(grid)

#read data into variable data
mds.data <- read.csv('Table.S1.csv', stringsAsFactors = T)
rownames(mds.data) <- mds.data[,1]
mds.metadata <- mds.data[,1:12]
mds.res_plas_cols <- c(19:128)
mds.resplas_data <- unique(mds.data[,mds.res_plas_cols])
mds.resplas_data <- mds.resplas_data[rowSums(mds.resplas_data) > 0,]
mds.vir_cols <- c(129:247)
mds.vir_data <- unique(mds.data[,mds.vir_cols])


#VIRULENCE
#create distance object 
virdist <- vegdist(wisconsin(mds.vir_data), method = "jaccard")
# Get the baseline solution: start with cmdscale
mds.null.vir <- isoMDS(virdist, tol=1e-7)
## See if you can get any better.
  repeat{
    mds.1.vir <- isoMDS(virdist, k=2, initMDS(virdist, k = 2), maxit=200, trace=FALSE, tol=1e-7)
    if(mds.1.vir$stress < mds.null.vir$stress) break
  }

# Scale solutions ("fix translation, rotation and scale")
mds.null.vir <- postMDS(mds.null.vir, virdist)
mds.1.vir <- postMDS(mds.1.vir, virdist)
# Compare solutions
plot(procrustes(mds.1.vir, mds.null.vir))
#convert to dataframe
mds.null.vir.df <- as.data.frame(mds.null.vir,
                             attr(x = mds.null.vir,
                                which = "dimnames"))
mds.1.vir.df <- as.data.frame(mds.1.vir, 
                          attr(x = mds.1.vir, 
                              which = "dimnames"))

rownames(mds.1.vir.df) -> mds.1.vir.df$Name
mds.1.vir.df <- left_join(mds.1.vir.df, mds.metadata)

mds.1.vir.df$ST <- as.factor(mds.1.vir.df$ST)


#RESISTANCE AND PLASMID
#create distance object 
resplasdist <- vegdist(wisconsin(mds.resplas_data), method = "jaccard")
# Get the baseline solution: start with cmdscale
mds.null.resplas <- isoMDS(resplasdist, tol=1e-7)
## See if you can get any better.
repeat{
  mds.1.resplas <- isoMDS(resplasdist, k=2, initMDS(resplasdist, k = 2), maxit=200, trace=FALSE, tol=1e-7)
  if(mds.1.resplas$stress < mds.null.resplas$stress) break
}

# Scale solutions ("fix translation, rotation and scale")
mds.null.resplas <- postMDS(mds.null.resplas, resplasdist)
mds.1.resplas <- postMDS(mds.1.resplas, resplasdist)
# Compare solutions
plot(procrustes(mds.1.resplas, mds.null.resplas))
#convert to dataframe
mds.null.resplas.df <- as.data.frame(mds.null.resplas,
                                 attr(x = mds.null.resplas,
                                      which = "dimnames"))
mds.1.resplas.df <- as.data.frame(mds.1.resplas, 
                              attr(x = mds.1.resplas, 
                                   which = "dimnames"))

rownames(mds.1.resplas.df) -> mds.1.resplas.df$Name
mds.1.resplas.df <- left_join(mds.1.resplas.df, mds.metadata)

mds.1.resplas.df$ST <- as.factor(mds.1.resplas.df$ST)


ST_cols <- c("#c70034",
             "#00c53a",
             "#953cec",
             "#f6e65d",
             "#0066fa",
             "#ff790f",
             "#ff34dd",
             "#844e37")
ST_breaks <- sort(unique(factor(mds.metadata$ST)))

cont_cols <- c("#c70034",
               "#ff790f",
               "#00c53a",
               "#0066fa",
               "#ff34dd")
cont_breaks <- sort(unique(factor(mds.metadata$Continent)))

origin_cols <- c("#c70034",
                 "#ff34dd",
                 "#00c53a",
                 "#0066fa")
origin_breaks <- sort(unique(factor(mds.metadata$Origin)))



#Plots
#Virulence by Origin
vir_origin <- ggplot() +
  geom_point(data=mds.1.vir.df,
             aes(x=points.1,
                 y=points.2,
                 color = Origin, shape = Origin), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.vir.df$points.1,
                   y=mds.1.vir.df$points.2, 
                   color = mds.1.vir.df$Origin)) +
  ggtitle("Nonmetric MDS Virulence Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_shape_manual(values = c(17,12,19,15))+
  scale_color_manual(values= origin_cols, breaks = origin_breaks)+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

#Virulence by Sequence Type
vir_ST <- ggplot() +
  geom_point(data=mds.1.vir.df,
             aes(x=points.1,
                 y=points.2,
                 color = ST), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.vir.df$points.1,
                   y=mds.1.vir.df$points.2, 
                   color = mds.1.vir.df$ST)) +
  scale_color_manual(breaks = ST_breaks, values = ST_cols)+
  ggtitle("Nonmetric MDS Virulence Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

#Virulence by Continent
vir_cont <- ggplot() +
  geom_point(data=mds.1.vir.df,
             aes(x=points.1,
                 y=points.2,
                 color = Continent, shape = Continent), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.vir.df$points.1,
                   y=mds.1.vir.df$points.2, 
                   color = mds.1.vir.df$Continent)) +
  ggtitle("Nonmetric MDS Virulence Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_shape_manual(values = c(18,17,15,19,12))+
  scale_color_manual(values= cont_cols, breaks = cont_breaks)+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

#Resistance and Plasmid by Origin
resplas_origin <- ggplot() +
  geom_point(data=mds.1.resplas.df,
             aes(x=points.1,
                 y=points.2,
                 color = Origin, shape = Origin), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.resplas.df$points.1,
                   y=mds.1.resplas.df$points.2, 
                   color = mds.1.resplas.df$Origin)) +
  ggtitle("Nonmetric MDS Resistance and Plasmid Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_shape_manual(values = c(17,12,19,15))+
  scale_color_manual(values= origin_cols, breaks = origin_breaks)+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

#Resistance and Plasmid by Sequence Type
resplas_ST <- ggplot() +
  geom_point(data=mds.1.resplas.df,
             aes(x=points.1,
                 y=points.2,
                 color = ST), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.resplas.df$points.1,
                   y=mds.1.resplas.df$points.2, 
                   color = mds.1.resplas.df$ST
  )) +
  scale_color_manual(breaks = ST_breaks, values = ST_cols)+
  ggtitle("Nonmetric MDS Resistance and Plasmid Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

#Resistance and Plasmid by Continent
resplas_cont <- ggplot() +
  geom_point(data=mds.1.resplas.df,
             aes(x=points.1,
                 y=points.2,
                 color = Continent, shape = Continent), 
             size = 1.5) + # add the point markers
  stat_ellipse(aes(x=mds.1.resplas.df$points.1,
                   y=mds.1.resplas.df$points.2, 
                   color = mds.1.resplas.df$Continent
  )) +
  ggtitle("Nonmetric MDS Resistance and Plasmid Profile")+
  xlab("Dimension 1")+
  ylab("Dimension 2")+
  scale_shape_manual(values = c(18,17,15,19,12))+
  scale_color_manual(values= cont_cols, breaks = cont_breaks)+
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))

