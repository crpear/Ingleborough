
# Ingleborough

# Started 21/11/17

# Use R 3.3.3

####################
####################

# Notes & to dos 



################################################################################
# 0. Preamble ##################################################################
################################################################################

################################################################################
# Installing the necessary packages...
#install.packages("moments")

################################################################################
# Loading the necessary packages...


install.packages("geoR")
install.packages("moments")
install.packages("mvoutlier")
install.packages("shapefiles")
install.packages("RColorBrewer")
install.packages("chemometrics")
install.packages("vegan")
install.packages("gstat")
install.packages("classInt")
install.packages("robustbase")
#install.packages("pcaMethods")
install.packages("rgdal")
install.packages("car")
install.packages("MASS")
install.packages("boot")
install.packages("lattice")
install.packages("compositions")
#install.packages("Compositional")
install.packages("ggplot2")
install.packages("outliers")
install.packages("nlme")
install.packages("xts")
install.packages("forecast")
install.packages("pls")


require(geoR)
require(moments)
require(mvoutlier)
require(shapefiles)
require(RColorBrewer)
require(chemometrics)
require(vegan)
require(gstat)
require(classInt)
require(robustbase)
#require(pcaMethods)
require(rgdal)
require(car)
require(MASS)
require(boot)
require(lattice)
require(compositions)
#require(Compositional)
require(ggplot2)
require(outliers)
require(nlme)
require(xts)
require(forecast)
require(pls)


# Setting working directory....
setwd("D:\\6 November stats visit 2017\\Analysis")


################################################################################


# Scatterplots with size proportional to the correlations...
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y, use="pairwise.complete.obs"),2)
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#    text(0.5, 0.5, txt, cex = cex.cor * abs(cor(x, y, use="pairwise.complete.obs"))) # for proportional
    text(0.5, 0.5, 'yayandwow', cex = cex.cor) # not proportional
}

################################################################################




######################################################################################################
# 2. Input data ######################################################################################
######################################################################################################

# Read in data

Data.0 <- read.csv("Data_prep_v1_PH_CP_reduced_v1_biomass.csv", header=T)
colnames(Data.0)

# [1] "Sample_ID"      "Quadrat_1m2"    "Easting_1m2"    "Northing_1m2"   "Average"        "DryAverage"     "CS_VIS2_1m2"    "CS_VIS3_1m2"    "CS_VIS4_1m2"    "CS_RedEdge"    
#[11] "CS_NIR1_1m2"    "CS_NIR2_1m2"    "CS_NIR3_1m2"    "CS_NIR4_1m2"    "CS_NIR5_1m2"    "CS_NIR6_1m2"    "CS_NIR7_1m2"    "CS_NIR8_1m2"    "CS_NIR9_1m2"    "CS_SWIR1_1m2"  
#[21] "CS_SWIR2_1m2"   "Grassland_type"

attach(Data.0)

# Read in shapefile - TO DO
#shp <- readOGR(".","Middle_ffridd")
#summary(shp)

# Set coordinates of the 1m quadrat
Coords1 <- as.matrix(cbind(Easting_1m2,Northing_1m2))

# Dimensions
dim(Data.0)

# Look at the data
View(Data.0)



######################################################################################################
# 3. Univariate EDA - ################################################################################
######################################################################################################

################################################################################
# Summary statistics
summary(Data.0[,c(5:22)])




################################################################################
# Histograms...
X11(width=20,height=10)
par(mfrow=c(2,2))
hist(Average, freq=F)
hist(DryAverage, freq=F)
hist(CS_VIS2_1m2, freq=F)
hist(CS_VIS3_1m2, freq=F)


################################################################################
# Maps of data
Data.1.spdf <- SpatialPointsDataFrame(Data.0[,3:4], Data.0)

# For legend...
mypalette.1 <- brewer.pal(10,"Spectral")

# Coordinate summaries
summary(Data.0[,4:5])

# north arrow and scale bar...
extra1 = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(265710,371310), scale = 40)
extra2 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(265970,370930), scale = 80, fill=c("transparent","black"))
extra3 = list("sp.text", c(265970,370945), "0", cex=0.8)
extra4 = list("sp.text", c(266050,370945), "80 m", cex=0.8)
#data.layout <- list(extra1,extra2,extra3,extra4,shp)
data.layout <- list(extra1,extra2,extra3,extra4)

X11(width=10,height=10)
spplot(Data.1.spdf, "Average", key.space = "right", pch=20, cex=0.75,
col.regions=rev(mypalette.1), cuts=9, par.settings=list(fontsize=list(text=11)),
main=list(label="Study sites",cex=1.25),sp.layout=data.layout)



######################################################################################################
# 4. Bivariate EDA - #################################################################################
######################################################################################################

corr.mat <- round(cor(Data.0[,c(5:21)]),2)
corr.mat


X11(width=18,height=12)
pairs(Data.0[,c(5:21)], cex=0.75, pch=19, upper.panel=panel.cor, cex.cor=1)

################################################################################
# Conditional boxplots

X11(width=14,height=4)
par(mfrow=c(1,2))
boxplot(Average~rev(factor(Grassland_type)), main="Average")
boxplot(DryAverage~rev(factor(Grassland_type)), main="Dry Average")



######################################################################################################
# 5. Multivariate EDA - REGRESSIONS ##################################################################
######################################################################################################


# Stanadrd linear model (OLS)
Average.model.1 <- lm(Average ~ CS_VIS2_1m2 + CS_VIS3_1m2 + CS_VIS4_1m2 + CS_RedEdge
+ CS_NIR1_1m2 + CS_NIR2_1m2 + CS_NIR3_1m2 + CS_NIR4_1m2 + CS_NIR5_1m2 + CS_NIR6_1m2
+ CS_NIR7_1m2 + CS_NIR8_1m2 + CS_NIR9_1m2 + CS_SWIR1_1m2 + CS_SWIR2_1m2 + factor(Grassland_type))

DryAverage.model.1 <- lm(DryAverage ~ CS_VIS2_1m2 + CS_VIS3_1m2 + CS_VIS4_1m2 + CS_RedEdge
+ CS_NIR1_1m2 + CS_NIR2_1m2 + CS_NIR3_1m2 + CS_NIR4_1m2 + CS_NIR5_1m2 + CS_NIR6_1m2
+ CS_NIR7_1m2 + CS_NIR8_1m2 + CS_NIR9_1m2 + CS_SWIR1_1m2 + CS_SWIR2_1m2 + factor(Grassland_type))

summary(Average.model.1)
summary(DryAverage.model.1)

X11(width=14,height=4)
par(mfrow=c(1,2))
plot(Average.model.1$fitted, Average)
plot(DryAverage.model.1$fitted, DryAverage)




####################################################################################################
# 6. Multivariate EDA - Basic PCA on image data only ###############################################
####################################################################################################


# IMPORTANT doing all analysis on standardised data (unit variance)
# This means that all PCA and GWPCA models should specify the covariance matrix (not correlation matrix)
# We do this now as GWPCA only has covariance option currently.

Data.1.scaled <- scale(as.matrix(Data.0[,7:21])) # standardised data... (scaling and centering)
#Data.1.scaled <- as.matrix(Data.0[,7:21])
rownames(Data.1.scaled) <- Data.1[,1]
n1 <- length(Data.1.scaled[,1])
n1
#fix(Data.1.scaled)


# Basic PCA...
cat("\n\nGlobal Principal Components Analysis\n")
cat(    "====================================\n")
#pca <- princomp(Data.1.scaled,cor=T,scores=T) # correlation matrix - don't use as we are using standardised data
pca <- princomp(Data.1.scaled,cor=F,scores=T) # use covariance matrix to match the following...
pca$sdev # the square root of the eigenvalues
pca$loadings
pca$scores

# Alternative PCA calculations using prcomp...
pca1 <- prcomp(Data.1.scaled) # use covariance matrix
pca1$sdev # the square root of the eigenvalues
pca1$rotation # loadings

# So matching it all up....
pca$sdev
pca1$sdev

# Actual eigenvalues are:
pca$sdev^2
pca1$sdev^2

# loadings match ups - e.g. ....
pca$loadings[,1]
pca1$rotation[,1]
 
# With respect to variance explained then this matches up - but this relates to standard deviation
pca_prop <- pca$sdev / sum(pca$sdev)
pca_prop
pca_prop1 <- pca1$sdev / sum(pca1$sdev)
pca_prop1

# With respect to variance explained then this matches up - this relates to variance....
# note sum(pca$sdev^2) should add up to 18 as 18 variables each with unit variance...
pca_prop <- pca$sdev^2 / sum(pca$sdev^2)
pca_prop
pca_prop1 <- pca1$sdev^2 / sum(pca1$sdev^2)
pca_prop1

# Outputting global PCA data to excel...
pca_v <- rbind(pca$sdev^2,pca_prop)
rownames(pca_v) <- c("Eigenvalues","Proportions")
print(pca_v)
#write.table(pca_v,"clipboard",sep="\t",col.names=T,row.names=T)

# Cumulative variance explained....
cpca_sum <- cumsum(pca_v[2,])
print(cpca_sum)
#write.table(cpca_sum,"clipboard",sep="\t",col.names=T,row.names=T)

print(pca$loadings)
#write.table(pca$loadings,"clipboard",sep="\t",col.names=T,row.names=T)




######################################################################################################
# 7. Multivariate EDA - PCA REGRESSIONS ##############################################################
######################################################################################################


# Stanadrd linear model (OLS)
Average.model.2 <- lm(Average ~ pca$scores[,1] + pca$scores[,2] + pca$scores[,3] + pca$scores[,4] + factor(Grassland_type))

DryAverage.model.2 <- lm(DryAverage ~ pca$scores[,1] + pca$scores[,2] + pca$scores[,3] + pca$scores[,4] + factor(Grassland_type))

summary(Average.model.2)
summary(DryAverage.model.2)

X11(width=14,height=4)
par(mfrow=c(1,2))
plot(Average.model.2$fitted, Average)
plot(DryAverage.model.2$fitted, DryAverage)




######################################################################################################
# 8. Multivariate EDA - Partial Least Squares REGRESSIONS ############################################
######################################################################################################

# PLS model
Average.model.3 <- plsr(Average ~ CS_VIS2_1m2 + CS_VIS3_1m2 + CS_VIS4_1m2 + CS_RedEdge
+ CS_NIR1_1m2 + CS_NIR2_1m2 + CS_NIR3_1m2 + CS_NIR4_1m2 + CS_NIR5_1m2 + CS_NIR6_1m2
+ CS_NIR7_1m2 + CS_NIR8_1m2 + CS_NIR9_1m2 + CS_SWIR1_1m2 + CS_SWIR2_1m2 + factor(Grassland_type), ncomp = 4, data = Data.0)

DryAverage.model.3 <- plsr(DryAverage ~ CS_VIS2_1m2 + CS_VIS3_1m2 + CS_VIS4_1m2 + CS_RedEdge
+ CS_NIR1_1m2 + CS_NIR2_1m2 + CS_NIR3_1m2 + CS_NIR4_1m2 + CS_NIR5_1m2 + CS_NIR6_1m2
+ CS_NIR7_1m2 + CS_NIR8_1m2 + CS_NIR9_1m2 + CS_SWIR1_1m2 + CS_SWIR2_1m2 + factor(Grassland_type), ncomp = 4, data = Data.0)

summary(Average.model.3)
summary(DryAverage.model.3)

# Need to sort out output




######################################################################################################
# 9. END #############################################################################################
######################################################################################################



















