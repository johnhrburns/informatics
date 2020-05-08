# PCA of ground water data

gwd <- read.csv("gwd_all_raw.csv", row.names = 1) # Water quality

View(gwd)

gwd[is.na(gwd)] <- 0

View(gwd)

gwd.norm <- decostand(gwd, "normalize")

# PCA with plot

library(ade4)

library(vegan)

library(gclus)

library(ape)

source("panelutils.R")

gwd.pca <- rda(gwd, scale = TRUE)

gwd.pca

source("cleanplot.pca.R")

cleanplot.pca(gwd.pca, point = TRUE) 

# PCA with clusters

gwd.w <- hclust(dist(scale(gwd)), "ward")

gwd.w

head(gwd.w)

    ## Cut the dendrogram to yield 4 groups

grp <- cutree(gwd.w, k = 4)

grp1 <- levels(factor(grp))

    ## Get site scores using "Scaling 1" method

site.sc <- scores(gwd.pca, display = "wa", scaling =1)

    ## Plot the sites with cluster symbols and colors using "Scaling 1" method

par(op)

p <- plot(gwd.pca, display = "wa", scaling = 3, type = "n", main = "PCA correlation with clusters")

for (i in 1:length(grp1)) {points(site.sc[grp==i,], pch = (14 + i), cex = 2, col = i+1)}

text(site.sc, row.names(gwd), cex = 0.7, pos = 3)

ordihull(site.sc, groups = grp, draw = "polygon", lty = 1, col = "grey90") # pick this one or the ellipse

gwd.fit <- envfit(gwd.pca, gwd, permutations = 999)

plot(gwd.fit, p.max = 0.01, col = "red", cex = 0.7)


##
#
#
# Cluster analysis

    # Use Ward clustering analysis to identify four groups

gwd.ch <- vegdist(gwd.norm, "euc") # generate chord distances for normalized species data

gwd.ch.ward <- hclust(gwd.ch, method="ward") # compute ward's minimum variance clustering

gwd.ch.ward 

par(op)

plot(gwd.ch.ward)

k <- 4 # we are **choosing** to create 4 groups (*this is heuristic!*)

gwd.ward.g <- cutree(gwd.ch.ward, k)

plot(gwd.ward.g)

dend2 <- as.dendrogram(gwd.ch.ward)

heatmap(as.matrix(gwd.ch), symm = TRUE, Rowv = dend2, margin=c(9,5))

?heatmap

library(factoextra)

  # New Wards method

gwd.clus <- hclust(gwd.ch, method = "ward.D2")

gwd.grp <- cutree(gwd.clus, k = 4)

table(gwd.grp)

gwd.grp

plot(gwd.clus, cex = 0.6)
rect.hclust(gwd.clus, k = 4, border = 2:5)

write.csv(gwd.grp, "/Users/johnhrburns/Desktop/MARE_375_Spring_2020/Modules/Ecoinformatics/gwd.grp.csv")

#
#
#
###
## RDA

# Normalize

gwd.pred <- read.csv("gwd_pred.csv", row.names = 1)

gwd.pred[is.na(gwd.pred)] <- 0

gwd.norm <- decostand(gwd, "normalize")

pred.norm <- decostand(gwd.pred, "normalize")

# Run RDA on wq data

gwd.rda <- rda(gwd.norm ~ ., gwd.pred)

summary(gwd.rda)

# R2

gwdR2adj <- RsquareAdj(gwd.rda)$adj.r.squared

gwdR2adj 

# RDA plot

plot(gwd.rda, scaling = 1, main = "Triplot RDA - scaling 1")

gwd.sc <- scores(gwd.rda, choices = 1:2, scaling = 1, display = "sp")

arrows(0, 0, gwd.sc[,1], gwd.sc[,2], length = 0, lty = 1, col = "red")

# Global test of the RDA results:

anova.cca(gwd.rda, step = 1000)

# Tests of all canonical axes:

anova.cca(gwd.rda, by = "axis", step = 1000)

install.packages("remotes")
remotes::install_github("fawda123/ggord")

# New ways to plot the RDA

library(ggord)
ggord(gwd.rda) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gwdgroups <- read.csv("gwd.grp.csv")
levels(gwdgroups$Group) <- c("Group 1", "Group 2", "Group 3", "Group 4")
group <- gwdgroups$Group


gwd.sc <- scores(gwd.rda, choices = 1:2, scaling = 3, display = "sp")

bg <- c("#ff7f00","#1f78b4","#ffff33","#33a02c") # 4 nice colors for groups

plot(gwd.rda, type="n", scaling=3)
text(gwd.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the gwd data
points(gwd.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[group]) # the sites/groups
text(gwd.rda, display="sites", pch=21, cex=0.5, col="gray32", scaling=3, bg=bg[group]) # the sites/groups
text(gwd.rda, scaling=3, display="bp", col="#0868ac", cex=.8)  # the predictors
arrows(0, 0, gwd.sc[,1], gwd.sc[,2], length = 0, lty = 1, col = "#e31a1c")
legend("topleft", legend=levels(group), bty="n", col="gray32", pch=21, cex=1.0, pt.bg=bg)


#Error in match.arg(display, c("sites", "species", "wa", "lc", "bp", "cn",  : 
                                'arg' should be one of “sites”, “species”, “wa”, “lc”, “bp”, “cn”, “reg”

#
#

#
#
# PCoA - NOT WORKING - Bray curtis issues

gwd.bray <- vegdist(gwd)

gwd.pcoa <- cmdscale(gwd.bray, k=(nrow(spe)-1), eig=TRUE)

# Plot the sites and weighted average projection of species

ordiplot(scores(gwd.pcoa) [,c(1,2)], type="t", main="PCoA")

abline(h=0, lty=3)

abline(v=0, lty=3)

gwd.wa <- wascores(gwd.pcoa$points[,1:2], gwd)

text(gwd.wa, rownames(gwd.wa), cex=0.7, col="red")



#
#
# Nonmetric Multidimensional Scaling with GW data

gwd.nmds <- metaMDS(gwd, distance = "bray")

gwd.nmds

gwd.bray <- vegdist(gwd)

gwd.nmds$stress

plot(gwd.nmds, type = "t", main = paste("NMDS/Bray - Stress=", round(gwd.nmds$stress,3)))

# Stress plot and GOF

par(mfrow=c(1,2))

stressplot(gwd.nmds, main = "Shepard Plot")

gof = goodness(gwd.nmds)

plot(gwd.nmds, type = "t", main = "Goodness of fit")

points(gwd.nmds, display = "sites", cex = gof*200)

# Combine clustering with NMDS result

par(op)

gwd.bray.ward <- hclust(gwd.bray, "ward")

gwd.bw.groups <- cutree(gwd.bray.ward, k=4)

gwd.lev <- levels(factor(gwd.bw.groups))

site.sc <- scores(gwd.nmds)

# Add vectors for environmental data to the NMDS plot

# Make new plot and add sites

p <- ordiplot(site.sc, type = "n", main = "NMDS/Bray - clusters Ward/Bray")

for (i in 1:length(gwd.lev)) {points(site.sc[gwd.bw.groups==i,], pch=(3-i), cex=2, col=i+8)}

text(site.sc, row.names(gwd), pos=4, cex=0.7)

## Add vectors for environmental data to the NMDS plot

  # Add either hulls or ellipses and vectors

ordihull(gwd.nmds, groups = gwd.bw.groups, draw = "polygon", lty = 1, col = "grey90") # pick this one or the ellipse

ordiellipse(gwd.nmds, groups = gwd.bw.groups, draw = "polygon", lty = 1, col = "grey90") # pick this one or the hull

gwd.fit <- envfit(gwd.nmds, gwd, permutations = 999)

plot(gwd.fit, p.max = 0.01, col = "black", cex = 0.7)


