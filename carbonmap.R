##===================================================================
##
## Vieilledent et al. 2016 Bioclimatic envelope models predict a
## decrease in tropical forest carbon stocks with climate change in
## Madagascar. Journal of Ecology
##
## Author: Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
##
## This R script can be used to reproduce the results of the
## study. Associated data are available on the Dryad repository:
## http://dx.doi.org/10.5061/dryad.9ph68
##
##===================================================================

##=====================
##
## Preambule
##
##=====================

##= Libraries
list.of.packages = c("rgdal","sp","raster","rdryad","curl",
                     "ggplot2","gam","randomForest")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)

##= Functions

## R2 function
r2.calc <- function (obs,pred) {
  SCR <- sum((obs-pred)^2,na.rm=TRUE)
  SCT <- sum((obs-mean(obs))^2,na.rm=TRUE)
  R.square <- 1-SCR/SCT
  return(R.square)
}

## RMSE
RMSE.calc <- function (obs,pred) {
  SCR <- sum((obs-pred)^2,na.rm=TRUE)
  RMSE <- sqrt(SCR/length(obs))
  return(RMSE)
}

## Bias
bias.calc <- function(obs,pred) {
  return (100*(pred-obs)/obs)
}

## AIC
AIC.calc <- function(obs,pred) {
  sd.res <- sd(obs-pred)
  lL <- sum(log(dnorm(obs,pred,sd.res)))
  AIC <- -2*lL
  return(AIC)
}

# Unzip
decompress_file <- function(directory, file, .file_cache = FALSE) {
  
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    
    # Set working directory for decompression
    # simplifies unzip directory location behavior
    wd <- getwd()
    setwd(directory)
    
    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)
    
    # uncomment to delete archive once decompressed
    # file.remove(file) 
    
    # Reset working directory
    setwd(wd); rm(wd)
    
    # Test for success criteria
    # change the search depending on 
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
}    

##=====================
##
## Download data
##
##=====================

if (!dir.exists("data")) {
  dir.create("data")
}
f <- dryad_files("10.5061/dryad.9ph68")
for (i in c(1,6,7,8,9)) {
  name_f <- strsplit(basename(f[i]),"?",fixed=TRUE)[[1]][1]
  if (!file.exists(paste0("data/",name_f))) {
    curl_download(f[i],paste0("data/",name_f))
  }
}

##=====================
##
## Load ACD data
##
##=====================

## Load CSV file
ACD <- read.csv("data/ACD.csv",header=TRUE)
ACD.plot <- SpatialPointsDataFrame(coords=cbind(ACD$LON,ACD$LAT),
                                   data=ACD,
                                   proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
                                   #proj4string=CRS("+init=epsg:4326"))
## Plot ACD points
if (!dir.exists("results")) {
  dir.create("results")
}
pdf("results/ACDplots.pdf")
plot(ACD.plot)
dev.off()
## Data-set
df <- ACD.plot@data
str(df)

##====================================================
##
## Maps describing ecoregions and forest plot sampling
##
##====================================================

##=====================
## Ecoregion and plots

## Reproject ACD.plot
# ACD.38s <- spTransform(ACD.plot,CRS("+init=epsg:32738"))
ACD.38s <- spTransform(ACD.plot,CRS("+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs"))
## Forest raster
Zvar <- stack("data/Zvar.tif")
names(Zvar) <- c("forest","ecoregion","Saatchi250","Baccini250")
forest <- Zvar[[1]]
## Biome shapefile
if (!dir.exists("data/ecoregions_shp")) {
  decompress_file(directory="data", file="ecoregions_shp.zip")
}
Biomes.v <- readOGR(dsn="data/ecoregions_shp/",layer="madagascar_ecoregion_tenaizy_38s")
## Colors
red.t <- adjustcolor("red",alpha.f=0.5)
blue.t <- adjustcolor("blue",alpha.f=0.5)
green.t <- adjustcolor("dark green",alpha.f=0.5)
orange.t <- adjustcolor("orange",alpha.f=0.5)
## Plot
pdf("results/Ecoregion_And_Plots.pdf",width=5,height=9)
par(mar=c(0,0,0,0))
plot(Biomes.v,border=grey(0.3),col=c(red.t,blue.t,green.t,orange.t))
plot(forest,col="black",add=TRUE,legend=FALSE)
plot(ACD.38s,add=TRUE,pch=16)
dev.off()

##====================================================
## Mean annual temp. and annual precip. by forest type

## Cells without forest cover
w.notfor10 <- which(is.na(values(forest)))
## Biomes raster
Biomes <-  Zvar[[2]]
## Bioclimatic rasters
Xvar <- stack("data/Xvar.tif")
names(Xvar) <- c("ELEV","EVI","VCF","TEMP","TSEAS","PRECIP")
bio1 <- Xvar[[4]]
bio4 <- Xvar[[5]]
bio12 <- Xvar[[6]]

##= Moist
w.not.moist <- which(values(Biomes)!=3)
bio1.moist <- bio1 
values(bio1.moist)[unique(w.not.moist,w.notfor10)] <- NA
w.notNA.moist <- which(!is.na(values(bio1.moist)))
s.moist <- sample(w.notNA.moist,size=20000,replace=FALSE)
df.moist <- data.frame(t=bio1.moist[s.moist],p=bio12[s.moist])

##= Dry
w.not.dry <- which(values(Biomes)!=4)
bio1.dry <- bio1
values(bio1.dry)[unique(w.not.dry,w.notfor10)] <- NA
w.notNA.dry <- which(!is.na(values(bio1.dry)))
s.dry <- sample(w.notNA.dry,size=20000,replace=FALSE)
df.dry <- data.frame(t=bio1.dry[s.dry],p=bio12[s.dry])

##= Spiny
w.not.spiny <- which(values(Biomes)!=1)
bio1.spiny <- bio1
values(bio1.spiny)[unique(w.not.spiny,w.notfor10)] <- NA
w.notNA.spiny <- which(!is.na(values(bio1.spiny)))
s.spiny <- sample(w.notNA.spiny,size=20000,replace=FALSE)
df.spiny <- data.frame(t=bio1.spiny[s.spiny],p=bio12[s.spiny])

##= Forest plots
df.plot <- data.frame(t=ACD.plot$TEMP,p=ACD.plot$PRECIP)

##= Plot
pdf("results/ForestZoneBiomes.pdf",width=10,height=10)
ggplot(df.moist, aes(x=p, y=t)) +
  geom_point(data=df.plot,col=grey(0.5)) + xlim(250,3100) + ylim(130,280) +
  geom_density2d(col="dark green") +
  geom_density2d(data=df.dry,col="orange") +
  geom_density2d(data=df.spiny,col="red") +
  labs(x="\nAnnual precipitation (mm.y-1)",y="Mean annual temperature (°C x 10)\n",size=2) +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
dev.off()

##=========================================================================
##
## Descriptive analysis: relationship between ACD and explicative variables
##
##=========================================================================

pdf("results/descriptive_plots.pdf")
qplot(TEMP,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
qplot(TSEAS,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
qplot(PRECIP,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
qplot(ELEV,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
qplot(EVI,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
qplot(VCF,ACDobs,data=df,geom=c("point","smooth"),method="lm",formula=y~ns(x,3))
dev.off()

##=====================================
##
## Statistical model with randomForest
##
##=====================================

##= Random Forests
varnames <- c("ELEV","EVI","VCF","TEMP","TSEAS","PRECIP")
f <- formula(paste0("ACDobs~",paste(varnames,collapse="+")))
set.seed(1408) ## To make computation reproducible (1408 = my birthday :))
mod.RF <- randomForest(f,data=df,importance=TRUE)

## Save the randomForest model
save(mod.RF,file="results/RFmodel.rda",compress=TRUE)

## Variable importance
sink(file="results/VarImportance.txt")
importance(mod.RF,type=1)
sink()

## Variable effect 
pdf(file="results/VarEffect.pdf",width=15,height=15)
par(mfrow=c(3,2),cex=1.4,cex.lab=1.2,font.lab=2)
par(mar=c(5,5,2,3))
partialPlot(mod.RF,df,x.var="ELEV",xlab="Elevation (E, in m)",ylab="ACD (Mg.ha-1)",main="",lwd=3)
legend("topright",legend="%IncMSE = 45%",bty="n",cex=1.2)
par(mar=c(5,1,2,3))
partialPlot(mod.RF,df,x.var="EVI",xlab="EVI",main="",lwd=3)
legend("topleft",legend="%IncMSE = 40%",bty="n",cex=1.2)
par(mar=c(5,5,2,3))
partialPlot(mod.RF,df,x.var="VCF",xlab="VCF (%)",ylab="ACD (Mg.ha-1)",main="",lwd=3)
legend("topleft",legend="%IncMSE = 38%",bty="n",cex=1.2)
par(mar=c(5,1,2,3))
partialPlot(mod.RF,df,x.var="TEMP",xlab="Temp. mean (T, in ° C x 10)",main="",lwd=3)
legend("topright",legend="%IncMSE = 48%",bty="n",cex=1.2)
par(mar=c(5,5,2,3))
partialPlot(mod.RF,df,x.var="TSEAS",xlab="Temp. seas. (S, in sd x 100)",ylab="ACD (Mg.ha-1)",main="",lwd=3)
legend("topright",legend="%IncMSE = 68%",bty="n",cex=1.2)
par(mar=c(5,1,2,3))
partialPlot(mod.RF,df,x.var="PRECIP",xlab="Precip. (P, in mm.y-1)",main="",lwd=3)
legend("topleft",legend="%IncMSE = 49%",bty="n",cex=1.2)
dev.off()

##=====================================
##
## Carbon maps for the year 2010
##
##=====================================

##= Import raster stack of predictive variables
Xvar <- stack("data/Xvar.tif")
names(Xvar) <- c("ELEV","EVI","VCF","TEMP","TSEAS","PRECIP")
values(Xvar[[2]]) <- values(Xvar[[2]])/10000 # rescale EVI to be between 0 and 1
## Plot of predictive variables
pdf(file="results/Xvar.pdf")
plot(Xvar)
dev.off()

##= Predictions (caution: takes some time to run, ~5min)
ACD.pred <- predict(Xvar,model=mod.RF,
                    na.rm=TRUE,
                    filename="results/ACD_RF_2010.tif",
                    datatype="INT2S",
                    options=c("COMPRESS=LWZ", "PREDICTOR=2"),
                    format="GTiff",overwrite=TRUE,progress="text")
## ACD.pred <- raster("results/ACD_RF_2010.tif")

##= Mask with forest in 2010
Zvar <- stack("data/Zvar.tif")
names(Zvar) <- c("forest","ecoregion","Saatchi250","Baccini250")
forest <- Zvar[[1]]

##= Map of *forest* carbon stock in Madagascar (TIF format)
ACD.pred <- raster("results/ACD_RF_2010.tif")
ACD.for10 <- ACD.pred
w.notfor10 <- which(is.na(values(forest)))
values(ACD.for10)[w.notfor10] <- NA
writeRaster(ACD.for10,"results/ACD_for10.tif",datatype="INT2S",format="GTiff",
            options=c("COMPRESS=LWZ", "PREDICTOR=2"),
            overwrite=TRUE)

##= Plot of the carbon map (PDF format)
## Legend arguments
a.arg <- list(at=seq(0,250,length.out=6), labels=c("0","50","100","150","200","250"))
l.arg <- list(text="Carbon stock (Mg.ha-1)",side=2, line=0.5, cex=1.4)
## Colors
ter.col <- c(rev(terrain.colors(255))[-c(1:10)])
ColorRamp <- colorRampPalette(c(ter.col[length(ter.col)],"black"))
Colors <- c(rev(terrain.colors(255))[-c(1:10)],ColorRamp(50))
## Plot
pdf("results/Carbon_Map_Final.pdf",width=6.5,height=10)
par(mar=c(0,0,0,0),cex=1.4)
plot(ACD.for10,col=Colors,
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,
     axis.args=a.arg,zlim=c(0,250),legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()

##==========================================================
##
## Resampling predictive maps at 1km resolution
##
##==========================================================

## Region
Extent <- "298447.1 7155919 1100819 8682411"
Res <- "1000"
nodat <- "-9999"

## Resample with gdalwarp (GDAL must be installed on the computer)
Input <- "results/ACD_RF_2010.tif"
Output <- "results/ACD_RF_2010_resamp_1km.tif"
system(paste0("gdalwarp -ot Int16 -dstnodata ",nodat," -te ",Extent," -tr ",Res," ",Res," -r bilinear -overwrite ",Input," ",Output))

## Import resample raster
ACD.resamp <- raster("results/ACD_RF_2010_resamp_1km.tif")

## Extract predictions
df$ACDpredRF <- mod.RF$predicted ## randomForests predictions
df$ACDpred250m <- extract(ACD.pred,ACD.plot) ## Predictions at 250m
df$ACDpred1km <- extract(ACD.resamp,ACD.plot) ## Predictions at 1km

##==========================================================
##
## Model performance
##
##==========================================================

##= Obs Baccini
s.obs <- 1:nrow(df)
s.obs.Bac <- which(!is.na(df$BacciniOrig500m))

##= Performance indexes

## Random Forests
obs <- df$ACDobs[s.obs]
pred <- df$ACDpredRF[s.obs]
R2.RF <- r2.calc(obs,pred)
RMSE.RF <- RMSE.calc(obs,pred)
bias.RF <- bias.calc(obs,pred)
Bias.RF <- mean(bias.RF)
AIC.RF <- AIC.calc(obs,pred)
n.RF <- length(s.obs)

## Our map at 250m
obs <- df$ACDobs[s.obs]
pred <- df$ACDpred250m[s.obs]
R2.ACD.map <- r2.calc(obs,pred)
RMSE.ACD.map <- RMSE.calc(obs,pred)
bias.ACD.map <- bias.calc(obs,pred)
Bias.ACD.map <- mean(bias.ACD.map)
AIC.ACD.map <- AIC.calc(obs,pred)
n.ACD.map <- length(s.obs)

## Our map at 1km
obs <- df$ACDobs[s.obs]
pred <- df$ACDpred1km[s.obs]
R2.ACD.resamp <- r2.calc(obs,pred)
RMSE.ACD.resamp <- RMSE.calc(obs,pred)
bias.ACD.resamp <- bias.calc(obs,pred)
Bias.ACD.resamp <- mean(bias.ACD.resamp)
AIC.ACD.resamp <- AIC.calc(obs,pred)
n.ACD.resamp <- length(s.obs)

## Baccini at 1km
obs <- df$ACDobs[s.obs.Bac]
pred <- df$Baccini1km[s.obs.Bac]
R2.Bac <- r2.calc(obs,pred)
RMSE.Bac <- RMSE.calc(obs,pred)
bias.Bac <- bias.calc(obs,pred)
Bias.Bac <- mean(bias.Bac)
AIC.Bac <- AIC.calc(obs,pred)
n.Bac <- length(s.obs.Bac)

## Saatchi at 1km
obs <- df$ACDobs[s.obs]
pred <- df$Saatchi1km[s.obs]
R2.Saa <- r2.calc(obs,pred)
RMSE.Saa <- RMSE.calc(obs,pred)
bias.Saa <- bias.calc(obs,pred)
Bias.Saa <- mean(bias.Saa)
AIC.Saa <- AIC.calc(obs,pred)
n.Saa <- length(s.obs)

## Table of results
Perf <- data.frame(matrix(NA,nrow=5,ncol=6))
names(Perf) <- c("Model","R2","RMSE","Bias","AIC","n")
Perf$Model <- c("RF","ACD250m","ACD1km","Baccini1km","Saatchi1km")
Perf$R2 <- c(R2.RF,R2.ACD.map,R2.ACD.resamp,R2.Bac,R2.Saa)
Perf$RMSE <- c(RMSE.RF,RMSE.ACD.map,RMSE.ACD.resamp,RMSE.Bac,RMSE.Saa)
Perf$Bias <- c(Bias.RF,Bias.ACD.map,Bias.ACD.resamp,Bias.Bac,Bias.Saa)
Perf$AIC <- c(AIC.RF,AIC.ACD.map,AIC.ACD.resamp,AIC.Bac,AIC.Saa)
Perf$n <- c(n.RF,n.ACD.map,n.ACD.resamp,n.Bac,n.Saa)

## Backup
sink(file="results/Performance_Index_Internal.txt")
print(Perf)
sink()

##= Plot
pdf(file="results/Predicted_observed.pdf",width=10,height=10)
par(mfrow=c(2,2),cex.axis=1.4,cex.lab=1.5,cex.main=1.5)
obs <- df$ACDobs[s.obs]
## Mada 250m
pred <- df$ACDpred250m[s.obs]
par(mar=c(3,5,2,1))
plot(obs,pred,
     col=grey(0.5),
     main="Madagascar 250m",
     xlab="",
     ylab="Predicted ACD (Mg.ha-1)",
     axes=FALSE,
     xlim=c(0,350),
     ylim=c(0,350))
abline(a=0,b=1)
axis(1,at=seq(0,350,by=50),labels=seq(0,350,by=50))
axis(2,at=seq(0,350,by=50),labels=seq(0,350,by=50))
## Mada 1km
pred <- df$ACDpred1km[s.obs]
par(mar=c(3,2,2,1))
plot(obs,pred,
     col=grey(0.5),
     main="Madagascar 1km",
     xlab="",
     ylab="",
     axes=FALSE,
     xlim=c(0,350),
     ylim=c(0,350))
abline(a=0,b=1)
axis(1,at=seq(0,350,by=50),labels=seq(0,350,by=50))
axis(2,at=seq(0,350,by=50),labels=seq(0,350,by=50))
## Saatchi 2011
pred <- df$Saatchi1km[s.obs]
par(mar=c(5,5,2,1))
plot(obs,pred,
     col=grey(0.5),
     main="Saatchi 1km",
     xlab="Observed ACD (Mg.ha-1)",
     ylab="Predicted ACD (Mg.ha-1)",
     axes=FALSE,
     xlim=c(0,350),
     ylim=c(0,350))
abline(a=0,b=1)
axis(1,at=seq(0,350,by=50),labels=seq(0,350,by=50))
axis(2,at=seq(0,350,by=50),labels=seq(0,350,by=50))
## Baccini 2012
pred <- df$Baccini1km[s.obs]
par(mar=c(5,2,2,1))
plot(obs,pred,
     col=grey(0.5),
     main="Baccini 1km",
     xlab="Observed ACD (Mg.ha-1)",
     ylab="",
     axes=FALSE,
     xlim=c(0,350),
     ylim=c(0,350))
abline(a=0,b=1)
axis(1,at=seq(0,350,by=50),labels=seq(0,350,by=50))
axis(2,at=seq(0,350,by=50),labels=seq(0,350,by=50))
dev.off()

##= Plot of the bias
b.ACD.resamp <- lowess(df$ACDobs[s.obs],bias.ACD.resamp)
b.Bac <- lowess(df$ACDobs[s.obs.Bac],bias.Bac)
b.Saa <- lowess(df$ACDobs[s.obs],bias.Saa)

min.Bias <- min(b.ACD.resamp$y,b.Bac$y,b.Saa$y)
max.Bias <- max(b.ACD.resamp$y,b.Bac$y,b.Saa$y)

pdf(file="results/bias.pdf")
par(mar=c(4,4,1,1),cex=1.4)
plot(b.ACD.resamp$x,b.ACD.resamp$y,
     bty="n",
     axes=FALSE,
     col="black",
     type="l",
     xlim=c(0,350),
     ylim=c(-90,120),
     xlab="ACD (Mg.ha-1)",
     ylab="Bias (% of ACD)",
     lwd=2,
     lty=1)
axis(1,at=seq(from=0,to=350,by=50),labels=seq(from=0,to=350,by=50))
axis(2,at=seq(from=-90,to=120,by=30),labels=seq(from=-90,to=120,by=30))
segments(x0=0,y0=0,x1=350,y1=0,lwd=1,col=grey(0.7),lty=1)

lines(b.Bac$x,b.Bac$y,lty=2,lwd=2,col=grey(0.5))
lines(b.Saa$x,b.Saa$y,lty=4,lwd=2,col=grey(0.5))

legend(x=0,y=-50,legend=c("Madagascar 1km","Baccini 1km","Saatchi 1km"),
       lty=c(1,2,4),
       lwd=c(2,2,2),
       col=c("black",rep(grey(0.5),2)),
       bty="n",cex=0.8)
dev.off()

##=====================
##
## Cross-validation
##
##=====================

ncross <- 10
R2.RF <- RMSE.RF <- Bias.RF <- AIC.RF <- vector()
R2.Bac <- RMSE.Bac <- Bias.Bac <- AIC.Bac <- vector()
R2.Saa <- RMSE.Saa <- Bias.Saa <- AIC.Saa <- vector()

for (i in 1:ncross) {
  
  ##= Training and test df (subset of df 70/30)
  nobs <- nrow(df)
  train.obs <- sort(sample(nobs,size=floor(nobs*70/100),replace=FALSE))
  n.train <- length(train.obs)
  test.obs <- c(1:nobs)[-sort(train.obs)]
  n.test <- length(test.obs)
  df.train <- df[train.obs,]
  df.test <- df[test.obs,] 
  
  ##= Random Forests
  mod.RF <- randomForest(f,data=df.train)
  df.test$ACDpredRF <- predict(mod.RF,df.test)
  
  ##= Obs Baccini
  s.obs <- 1:nrow(df.test)
  s.obs.Bac <- which(!is.na(df.test$BacciniOrig500m))
  
  ##= Performance indexes
  
  ## Random Forests
  obs <- df.test$ACDobs[s.obs]
  pred <- df.test$ACDpredRF[s.obs]
  R2.RF[i] <- r2.calc(obs,pred)
  RMSE.RF[i] <- RMSE.calc(obs,pred)
  bias.RF <- bias.calc(obs,pred)
  Bias.RF[i] <- mean(bias.RF)
  AIC.RF[i] <- AIC.calc(obs,pred)
  
  ## Baccini
  obs <- df.test$ACDobs[s.obs.Bac]
  pred <- df.test$BacciniOrig500m[s.obs.Bac]
  R2.Bac[i] <- r2.calc(obs,pred)
  RMSE.Bac[i] <- RMSE.calc(obs,pred)
  bias.Bac <- bias.calc(obs,pred)
  Bias.Bac[i] <- mean(bias.Bac)
  AIC.Bac[i] <- AIC.calc(obs,pred)
  
  ## Saatchi
  obs <- df.test$ACDobs[s.obs]
  pred <- df.test$SaatchiOrig1km[s.obs]
  R2.Saa[i] <- r2.calc(obs,pred)
  RMSE.Saa[i] <- RMSE.calc(obs,pred)
  bias.Saa <- bias.calc(obs,pred)
  Bias.Saa[i] <- mean(bias.Saa)
  AIC.Saa[i] <- AIC.calc(obs,pred)
  
} ## End of bootstrap

## Table of results
Perf.mean <- data.frame(matrix(NA,nrow=3,ncol=5))
names(Perf.mean) <- c("Model","R2","RMSE","Bias","AIC")
Perf.mean$Model <- c("RF","Baccini","Saatchi")
Perf.mean$R2 <- apply(cbind(R2.RF,R2.Bac,R2.Saa),2,mean)
Perf.mean$RMSE <- apply(cbind(RMSE.RF,RMSE.Bac,RMSE.Saa),2,mean)
Perf.mean$Bias <- apply(cbind(Bias.RF,Bias.Bac,Bias.Saa),2,mean)
Perf.mean$AIC <- apply(cbind(AIC.RF,AIC.Bac,AIC.Saa),2,mean)

## Backup
sink(file="results/Performance_Index_CV_mean.txt")
print(Perf.mean)
sink()

## Table of results
Perf.sd <- data.frame(matrix(NA,nrow=3,ncol=5))
names(Perf.sd) <- c("Model","R2","RMSE","Bias","AIC")
Perf.sd$Model <- c("RF","Baccini","Saatchi")
Perf.sd$R2 <- apply(cbind(R2.RF,R2.Bac,R2.Saa),2,sd)
Perf.sd$RMSE <- apply(cbind(RMSE.RF,RMSE.Bac,RMSE.Saa),2,sd)
Perf.sd$Bias <- apply(cbind(Bias.RF,Bias.Bac,Bias.Saa),2,sd)
Perf.sd$AIC <- apply(cbind(AIC.RF,AIC.Bac,AIC.Saa),2,sd)

## Backup
sink(file="results/Performance_Index_CV_sd.txt")
print(Perf.sd)
sink()

##=============================================
##
## Forest carbon stock by forest type
##
##=============================================

##= Create matrix for export
df.ACD.biome <- data.frame(forest=c("Moist","Dry","Spiny"),
                           mean.ACD.plot=rep(0,3),sd.ACD.plot=rep(0,3),
                           mean.ACD.map=rep(0,3),sd.ACD.map=rep(0,3),
                           mean.ACD.Saatchi1km=rep(0,3),sd.ACD.Saatchi1km=rep(0,3),
                           mean.ACD.Baccini500m=rep(0,3),sd.ACD.Baccini500m=rep(0,3))

##= Observed plots
df.ACD.biome$mean.ACD.plot <- tapply(df$ACDobs,df$Ecoregion,mean)[c(2,1,3)]
df.ACD.biome$sd.ACD.plot <- tapply(df$ACDobs,df$Ecoregion,sd)[c(2,1,3)]

##= Load Biomes raster
Biomes <- Zvar[[2]]

##== Compute mean ACD by forest type
coef.surf <- 228.1411*228.1411/10000
surf.for10 <- round(sum(!is.na(ACD.for10[]))*coef.surf) ## in ha
ACD.moist <- ACD.dry <- ACD.spiny <- ACD.for10
##= Moist
w.not.moist <- which(values(Biomes)!=3)
ACD.moist[w.not.moist] <- NA
mean.ACD.moist <- mean(ACD.moist[],na.rm=TRUE)
sd.ACD.moist <- sd(ACD.moist[],na.rm=TRUE)
## Surface
surf.moist <- round(sum(!is.na(ACD.moist[]))*coef.surf) ## in ha
##= Dry
w.not.dry <- which(values(Biomes)!=4)
ACD.dry[w.not.dry] <- NA
mean.ACD.dry <- mean(ACD.dry[],na.rm=TRUE)
sd.ACD.dry <- sd(ACD.dry[],na.rm=TRUE)
## Surface
surf.dry <- round(sum(!is.na(ACD.dry[]))*coef.surf) ## in ha
##= Spiny
w.not.spiny <- which(values(Biomes)!=1)
ACD.spiny[w.not.spiny] <- NA
mean.ACD.spiny <- mean(ACD.spiny[],na.rm=TRUE)
sd.ACD.spiny <- sd(ACD.spiny[],na.rm=TRUE)
## Surface
surf.spiny <- round(sum(!is.na(ACD.spiny[]))*coef.surf) ## in ha

##= Total amount of C
TotACD.10 <- round(cellStats(ACD.for10,stat="sum",na.rm=TRUE)*coef.surf) ## 873,086,364 T of forest C

## Save results
df.ACD.biome$mean.ACD.map <- c(mean.ACD.moist,mean.ACD.dry,mean.ACD.spiny)
df.ACD.biome$sd.ACD.map <- c(sd.ACD.moist,sd.ACD.dry,sd.ACD.spiny)
df.ACD.biome$surf <- c(surf.moist,surf.dry,surf.spiny)

##=========================
##= For Saatchi

## Import
Saatchi250 <- Zvar[[3]]
values(Saatchi250)[w.notfor10] <- NA
ACD.moist <- ACD.dry <- ACD.spiny <- Saatchi250
## Moist
values(ACD.moist)[w.not.moist] <- NA
mean.ACD.moist <- mean(ACD.moist[],na.rm=TRUE)
sd.ACD.moist <- sd(ACD.moist[],na.rm=TRUE)
## Dry
values(ACD.dry)[w.not.dry] <- NA
mean.ACD.dry <- mean(ACD.dry[],na.rm=TRUE)
sd.ACD.dry <- sd(ACD.dry[],na.rm=TRUE)
## Spiny
values(ACD.spiny)[w.not.spiny] <- NA
mean.ACD.spiny <- mean(ACD.spiny[],na.rm=TRUE)
sd.ACD.spiny <- sd(ACD.spiny[],na.rm=TRUE)

## Save results
df.ACD.biome$mean.ACD.Saatchi1km <- c(mean.ACD.moist,mean.ACD.dry,mean.ACD.spiny)
df.ACD.biome$sd.ACD.Saatchi1km <- c(sd.ACD.moist,sd.ACD.dry,sd.ACD.spiny)

##= Total amount of C
TotACD.10.Saatchi <- df.ACD.biome$mean.ACD.Saatchi1km %*% df.ACD.biome$surf ## 764,167,76 T of forest C
## Another option:
## TotACD.10.Saatchi <- round(cellStats(Saatchi250,stat="sum",na.rm=TRUE)*coef.surf) ## 769,853,842 T of forest C

##=========================
##= For Baccini

## Import
Baccini250 <- Zvar[[4]]
values(Baccini250)[w.notfor10] <- NA
Baccini250[] <- 0.47*Baccini250[] ## Baccini's map is biomass !!
ACD.moist <- ACD.dry <- ACD.spiny <- Baccini250
## Moist
values(ACD.moist)[w.not.moist] <- NA
mean.ACD.moist <- mean(ACD.moist[],na.rm=TRUE)
sd.ACD.moist <- sd(ACD.moist[],na.rm=TRUE)
## Dry
values(ACD.dry)[w.not.dry] <- NA
mean.ACD.dry <- mean(ACD.dry[],na.rm=TRUE)
sd.ACD.dry <- sd(ACD.dry[],na.rm=TRUE)
## Spiny
values(ACD.spiny)[w.not.spiny] <- NA
mean.ACD.spiny <- mean(ACD.spiny[],na.rm=TRUE)
sd.ACD.spiny <- sd(ACD.spiny[],na.rm=TRUE)

## Save results
df.ACD.biome$mean.ACD.Baccini500m <- c(mean.ACD.moist,mean.ACD.dry,mean.ACD.spiny)
df.ACD.biome$sd.ACD.Baccini500m <- c(sd.ACD.moist,sd.ACD.dry,sd.ACD.spiny)

##= Total amount of C
TotACD.10.Baccini <- df.ACD.biome$mean.ACD.Baccini500m %*% df.ACD.biome$surf ## 652,895,945 T of forest C
TotACD <- data.frame(Source=c("ACD.map","Saatchi1km","Baccini500m"),TotACD=c(TotACD.10,TotACD.10.Saatchi,TotACD.10.Baccini))

##===============
## Export results
sink("results/mean_ACD_forest.txt")
print(df.ACD.biome)
sink()

sink("results/Tot_ACD_forest.txt")
print(TotACD)
sink()

##=====================================================================================
##
## Effect of climate change on ACD
##
##=====================================================================================

##= Unzip climatic projection data (takes some time to run: ~2min)
##= Climatic projections are unzipped in the ./climproj folder in the working directory
if (!dir.exists("data/climproj")) {
  decompress_file(directory="data",file="climproj.zip")
}

##= Loop indices
yr <- c("2050","2080") ## For 2050, 2080
mod <- c("ac","gs","ip","mc","he","cc","no") ## For global climate models (GCMs): CNRM-CM5, GISS-E2-R, HadGEM2-ES 
rcp <- c("45","85") ## For representative concentration pathways (RCPs): RCP 4.5, RCP 8.5
biovar <- c("bio1","bio4","bio12") ## BIO1=Annual Mean Temperature, BIO4=Temperature Seasonality (standard deviation *100)
## BIO12=Annual Precipitation

##= Total number of models
n.mod <- length(mod)*length(rcp)*length(yr)

##= Mask with forest in 2010
forest <- Zvar[[1]]
w.notfor10 <- which(is.na(forest[]))
coef.surf <- 228.1411*228.1411/10000

##= Data-frame of results
df.acd <- data.frame(MOD=rep(mod,each=4),RCP=rep(c("45","45","85","85"),7),YR=rep(yr,14),ACD=NA)

##= Loop (takes a long time to run: 28 runs in total, ~5min by run)
for (i in 1:length(mod)) {
  for (j in 1:length(rcp)) {
    for (l in 1:length(yr)) {
      ##= Message
      i.mod <- (i-1)*length(rcp)*length(yr) + (j-1)*length(yr) + l
      cat(paste0("\n","Model ",i.mod,"/",n.mod,": ",yr[l],mod[i],rcp[j],"\n"))
      ##= Explicative bioclimatic variables
      bio1.region <- raster(paste0("data/climproj/bio1_",mod[i],"_",rcp[j],"_",yr[l],"_region.tif"))
      bio4.region <- raster(paste0("data/climproj/bio4_",mod[i],"_",rcp[j],"_",yr[l],"_region.tif"))
      bio12.region <- raster(paste0("data/climproj/bio12_",mod[i],"_",rcp[j],"_",yr[l],"_region.tif"))
      Xvar.proj <- stack(list(ELEV=Xvar[[1]],EVI=Xvar[[2]],VCF=Xvar[[3]],
                              TEMP=bio1.region,TSEAS=bio4.region,PRECIP=bio12.region))
      ##= Predict with randomForest model
      ACD.proj <- predict(Xvar.proj, model=mod.RF,
                          na.rm=TRUE,
                          filename=paste0("results/ACD_",mod[i],"_",rcp[j],"_",yr[l],".tif"),
                          datatype="INT2S",
                          options=c("COMPRESS=LWZ", "PREDICTOR=2"),
                          format="GTiff",overwrite=TRUE,progress="text")
      ##= National forest ACD
      ACD.proj.for10 <- ACD.proj
      values(ACD.proj.for10)[w.notfor10] <- NA
      df.acd$ACD[i.mod] <- round(cellStats(ACD.proj.for10,stat="sum",na.rm=TRUE)*coef.surf)         
    }
  }
}

##== Loss
TotACD.10 <- TotACD.10
df.acd$loss.acd <- -(TotACD.10-df.acd$ACD)
df.acd$loss.perc <- -round(100*(TotACD.10-df.acd$ACD)/TotACD.10,2)
df.acd$ACD.2010 <- TotACD.10

##== Backup df.acd
sink(file="results/df.acd.txt")
df.acd
sink()

##== For RCP 85
## 2050
Mean <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2050,4:6],2,mean)
Max <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2050,4:6],2,min) ## Inverse min -> Maximal loss
Min <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2050,4:6],2,max) ## Inverse max -> Minimal loss
df.acd50.rcp85 <- rbind(Mean,Max,Min) 
sink(file="results/df.acd50.rcp85.txt")
df.acd50.rcp85
sink()
## 2080
Mean <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2080,4:6],2,mean)
Max <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2080,4:6],2,min) ## Inverse min -> Maximal loss
Min <- apply(df.acd[df.acd$RCP==85 & df.acd$YR==2080,4:6],2,max) ## Inverse max -> Minimal loss
df.acd80.rcp85 <- rbind(Mean,Max,Min) 
sink(file="results/df.acd80.rcp85.txt")
df.acd80.rcp85
sink()

##== ACD raster in 2050 and 2080 for RCP 8 from model averaging
## 2050
ACD.50ac85 <- raster("results/ACD_ac_85_2050.tif")
ACD.50cc85 <- raster("results/ACD_cc_85_2050.tif")
ACD.50gs85 <- raster("results/ACD_gs_85_2050.tif")
ACD.50he85 <- raster("results/ACD_he_85_2050.tif")
ACD.50ip85 <- raster("results/ACD_ip_85_2050.tif")
ACD.50mc85 <- raster("results/ACD_mc_85_2050.tif")
ACD.50no85 <- raster("results/ACD_no_85_2050.tif")
Stack.2050.85 <- stack(ACD.50ac85,ACD.50cc85,ACD.50gs85,ACD.50he85,ACD.50ip85,ACD.50mc85,ACD.50no85)
## 2080
ACD.80ac85 <- raster("results/ACD_ac_85_2080.tif")
ACD.80cc85 <- raster("results/ACD_cc_85_2080.tif")
ACD.80gs85 <- raster("results/ACD_gs_85_2080.tif")
ACD.80he85 <- raster("results/ACD_he_85_2080.tif")
ACD.80ip85 <- raster("results/ACD_ip_85_2080.tif")
ACD.80mc85 <- raster("results/ACD_mc_85_2080.tif")
ACD.80no85 <- raster("results/ACD_no_85_2080.tif")
Stack.2080.85 <- stack(ACD.80ac85,ACD.80cc85,ACD.80gs85,ACD.80he85,ACD.80ip85,ACD.80mc85,ACD.80no85)

##= Mean
ACD.2050.85.mean <- mean(Stack.2050.85)
ACD.2080.85.mean <- mean(Stack.2080.85)

##= Export
## Remove non forest
ACD.2050.85.for10 <- ACD.2050.85.mean
values(ACD.2050.85.for10)[w.notfor10] <- NA
ACD.2080.85.for10 <- ACD.2080.85.mean
values(ACD.2080.85.for10)[w.notfor10] <- NA
## Write o disk
writeRaster(ACD.2050.85.for10,"results/ACD.2050.85.for10.tif",datatype="INT2S",
            options=c("COMPRESS=LWZ", "PREDICTOR=2"),
            format="GTiff",overwrite=TRUE)
writeRaster(ACD.2080.85.for10,"results/ACD.2080.85.for10.tif",datatype="INT2S",
            options=c("COMPRESS=LWZ", "PREDICTOR=2"),
            format="GTiff",overwrite=TRUE)

##============================================================================
##
## Plot maps of future forest carbon stocks under the effect of climate change
##
##============================================================================

##==================
##= Load ACD rasters

ACD.for10 <- raster("results/ACD_for10.tif")
ACD.2050.85.for10 <- raster("results/ACD.2050.85.for10.tif")
ACD.2080.85.for10 <- raster("results/ACD.2080.85.for10.tif")

##========================
##= Plots of carbon stocks

##pdf("results/Evolution_ACD_Climate.pdf",width=10,height=7)
png("results/Evolution_ACD_Climate.png",width=300*10,height=300*7,res=300)
par(mfrow=c(1,3))
## 2010
par(cex=1.2,mar=c(0,0,4,0))
plot(ACD.for10,col=Colors,
     maxpixels=50000,
     axes=FALSE,box=FALSE,legend=FALSE,
     zlim=c(0,250),
     main="2010")
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
mtext(text="C stock: 873 086",side=3,line=0,cex=1.2)
## 2050
par(cex=1.2,mar=c(0,0,4,0))
plot(ACD.2050.85.for10,col=Colors,
     maxpixels=50000,
     axes=FALSE,box=FALSE,legend=FALSE,
     zlim=c(0,250),
     main="2050")
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
mtext(text="C stock: 799 097 (-8%)",side=3,line=0,cex=1.2)
## 2080
par(cex=1.2,mar=c(0,0,4,0))
plot(ACD.2080.85.for10,col=Colors,
     maxpixels=50000,
     axes=FALSE,box=FALSE,legend=FALSE,
     zlim=c(0,250),
     main="2080")
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
mtext(text="C stock: 720 944 (-17%)",side=3,line=0,cex=1.2)
dev.off()

## Legend
l.arg <- list(text="Carbon stock (Mg.ha-1)",side=2, line=0.5, cex=1.4)
png("results/Evolution_ACD_Climate_legend.png",height=300*10,width=300*2,res=300)
par(mar=c(0,0,0,0))
plot(ACD.for10,col=Colors,
     maxpixels=50000,
     legend.only=TRUE, horizontal=FALSE,
     axis.args=a.arg,zlim=c(0,250),legend.args=l.arg)
dev.off()

##======================
##= Carbon stock changes

## Changes
ACD.2010.2050 <- ACD.2050.85.for10-ACD.for10
ACD.2010.2080 <- ACD.2080.85.for10-ACD.for10

##= Histograms
png("results/Hist_carbon_stock_change.png",height=300*10,width=300*5,res=300)
par(mfrow=c(2,1),cex=1.2)
## 2050
par(mar=c(3,4,2,1))
hist(values(ACD.2010.2050),
     xlim=c(-150,50),
     col=grey(0.5),
     main="2010-2050",
     xlab="")
## 2080
par(mar=c(5,4,2,1))
hist(values(ACD.2010.2080),
     xlim=c(-150,50), 
     col=grey(0.5),
     main="2010-2080",
     xlab="Carbon stock change (Mg.ha-1)")
dev.off()

##= Legend arguments
a.arg <- list(at=seq(-150,50,length.out=5), labels=c("-150","-100","-50","0","50"))
l.arg <- list(text="Carbon stock change (Mg.ha-1)",side=2, line=0.5, cex=1.4)
## Colors
ColorRamp <- colorRampPalette(c("black","#d73027","#ffffbf","#1a9850","black"))
Colors <- ColorRamp(301)

## Plots
png("results/Carbon_stock_change.png",width=300*7,height=300*7,res=300)
par(mfrow=c(1,2))
## 2050
par(cex=1.2,mar=c(0,0,0,0))
plot(ACD.2010.2050,col=Colors[1:201],
     axes=FALSE,box=FALSE,legend=FALSE,
     main="",
     ##maxpixels=50000,
     axis.args=a.arg,zlim=c(-150,50),legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
## 2080
par(cex=1.2,mar=c(0,0,0,0))
plot(ACD.2010.2080,col=Colors[1:201],
     axes=FALSE,box=FALSE,legend=FALSE,
     main="",
     ##maxpixels=50000,
     axis.args=a.arg,zlim=c(-150,50),legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()

## Legend
l.arg <- list(text="Carbon stock change (Mg.ha-1)",side=2, line=0.5, cex=1.4)
png("results/Carbon_stock_change_legend.png",height=300*10,width=300*2,res=300)
par(mar=c(0,0,0,0))
plot(ACD.2010.2080,col=Colors[1:201],
     maxpixels=50000,
     legend.only=TRUE, horizontal=FALSE,
     axis.args=a.arg,zlim=c(-150,50),legend.args=l.arg)
dev.off()

##============================================================================
##
## Percentage of change
##
##============================================================================

# Load ACD rasters
ACD.for10 <- raster("results/ACD_for10.tif")
ACD.2050.85.for10 <- raster("results/ACD.2050.85.for10.tif")
ACD.2080.85.for10 <- raster("results/ACD.2080.85.for10.tif")

# Percentage
ACD.2050.pch <- 100*(ACD.2050.85.for10-ACD.for10)/ACD.for10

# Refuge area
RA.2050 <- ACD.2050.pch
RA.2050[values(ACD.2050.pch) > -15 | values(ACD.2050.pch) < 15] <- 1
RA.2050[values(ACD.2050.pch) <= -15 | values(ACD.2050.pch) >= 15] <- 2

# Plot refuge areas
png("results/Refuge_areas.png",height=300*7,width=300*7,res=300)
par(mar=c(0,0,0,0))
plot(RA.2050,col=c("darkgreen","orange"),
     axes=FALSE,box=FALSE,legend=FALSE,
     main="")
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()

##============================================================================
##
## Maps of predicted climatic anomalies
##
##============================================================================

##= Temperature seasonality
## Actual
bio4.region <- Xvar[[5]]
## Mean in 2080
bio4.ac.2080 <- raster("data/climproj/bio4_ac_85_2080_region.tif")
bio4.cc.2080 <- raster("data/climproj/bio4_cc_85_2080_region.tif")
bio4.gs.2080 <- raster("data/climproj/bio4_gs_85_2080_region.tif")
bio4.he.2080 <- raster("data/climproj/bio4_he_85_2080_region.tif")
bio4.ip.2080 <- raster("data/climproj/bio4_ip_85_2080_region.tif")
bio4.mc.2080 <- raster("data/climproj/bio4_mc_85_2080_region.tif")
bio4.no.2080 <- raster("data/climproj/bio4_no_85_2080_region.tif")
Stack.bio4.2080 <- stack(c(bio4.ac.2080,bio4.cc.2080,bio4.gs.2080,bio4.he.2080,
                           bio4.ip.2080,bio4.mc.2080,bio4.no.2080))
bio4.2080 <- mean(Stack.bio4.2080)
## Anomalies
bio4.anomalies <- bio4.region
bio4.anomalies[] <- bio4.2080[]-bio4.region[]
values(bio4.anomalies)[w.notfor10] <- NA
range(bio4.anomalies[],na.rm=TRUE)
## Legend
a.arg <- list(at=c(-250,-125,0,125,250), labels=c("-250","-125","0","+125","+250"))
l.arg <- list(text="Temp. seas. anomalies (sd x 100)",side=3, line=0.5, cex=1.4)
## Colors
Colors <- cm.colors(255)
## Plot
##pdf("results/TS_anomalies.pdf",width=6.5,height=10)
png("results/TS_anomalies.png",width=100*6,height=100*10,res=100)
par(cex=1.2,mar=c(2,0,0,0))
plot(bio4.anomalies,col=Colors,
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=TRUE,horizontal=TRUE,
     axis.args=a.arg, zlim=c(-250,250), legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()
## Summary
summary(bio4.anomalies) ## +138 sd x 100

##= Precipitation
## Actual
bio12.region <- Xvar[[6]]
## Mean in 2080
bio12.ac.2080 <- raster("data/climproj/bio12_ac_85_2080_region.tif")
bio12.cc.2080 <- raster("data/climproj/bio12_cc_85_2080_region.tif")
bio12.gs.2080 <- raster("data/climproj/bio12_gs_85_2080_region.tif")
bio12.he.2080 <- raster("data/climproj/bio12_he_85_2080_region.tif")
bio12.ip.2080 <- raster("data/climproj/bio12_ip_85_2080_region.tif")
bio12.mc.2080 <- raster("data/climproj/bio12_mc_85_2080_region.tif")
bio12.no.2080 <- raster("data/climproj/bio12_no_85_2080_region.tif")
Stack.bio12.2080 <- stack(c(bio12.ac.2080,bio12.cc.2080,bio12.gs.2080,bio12.he.2080,
                            bio12.ip.2080,bio12.mc.2080,bio12.no.2080))
bio12.2080 <- mean(Stack.bio12.2080)
## Anomalies
bio12.anomalies <- bio12.region
bio12.anomalies[] <- bio12.2080[]-bio12.region[]
values(bio12.anomalies)[w.notfor10] <- NA
range(bio12.anomalies[],na.rm=TRUE)
## Legen
a.arg <- list(at=seq(-250,250,length.out=5), labels=c("-250","-125","0","+125","+250"))
l.arg <- list(text="Precip. anomalies (mm.y-1)",side=3, line=0.5, cex=1.4)
## Colors
Colors <- cm.colors(255)
## Plot
##pdf("results/P_anomalies.pdf",width=6.5,height=10)
png("results/P_anomalies.png",width=100*6,height=100*10,res=100)
par(cex=1.2,mar=c(2,0,0,0))
plot(bio12.anomalies,col=Colors,
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=TRUE,horizontal=TRUE,
     axis.args=a.arg, zlim=c(-250,+250), legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()
## Summary
summary(bio12.anomalies) ## -107 mm.y-1

##= Temp. mean.
## Actual
bio1.region <- Xvar[[4]]
## Mean in 2080
bio1.ac.2080 <- raster("data/climproj/bio1_ac_85_2080_region.tif")
bio1.cc.2080 <- raster("data/climproj/bio1_cc_85_2080_region.tif")
bio1.gs.2080 <- raster("data/climproj/bio1_gs_85_2080_region.tif")
bio1.he.2080 <- raster("data/climproj/bio1_he_85_2080_region.tif")
bio1.ip.2080 <- raster("data/climproj/bio1_ip_85_2080_region.tif")
bio1.mc.2080 <- raster("data/climproj/bio1_mc_85_2080_region.tif")
bio1.no.2080 <- raster("data/climproj/bio1_no_85_2080_region.tif")
Stack.bio1.2080 <- stack(c(bio1.ac.2080,bio1.cc.2080,bio1.gs.2080,bio1.he.2080,
                           bio1.ip.2080,bio1.mc.2080,bio1.no.2080))
bio1.2080 <- mean(Stack.bio1.2080)
## Anomalies
bio1.anomalies <- bio1.region
bio1.anomalies[] <- bio1.2080[]-bio1.region[]
values(bio1.anomalies)[w.notfor10] <- NA
range(bio1.anomalies[],na.rm=TRUE)
## Legend
a.arg <- list(at=seq(30,45,length.out=4), labels=c("+30","+35","+40","+45"))
l.arg <- list(text="Temp. mean anomalies (°C x 10)",side=3, line=0.5, cex=1.4)
## Colors
heat.col <- c(rev(heat.colors(255))[-c(1:10)])
ColorRamp <- colorRampPalette(c(heat.col[length(heat.col)],"black"))
Colors <- c(rev(heat.colors(255))[-c(1:10)],ColorRamp(50))
## Plot
##pdf("results/T_anomalies.pdf",width=6.5,height=10)
png("results/T_anomalies.png",width=100*6,height=100*10,res=100)
par(cex=1.2,mar=c(2,0,0,0))
plot(bio1.anomalies,col=Colors,
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=TRUE,horizontal=TRUE,
     axis.args=a.arg, zlim=c(30,45), legend.args=l.arg)
plot(Biomes.v,border=grey(0.3),col="transparent",add=TRUE)
dev.off()
## Summary
summary(bio1.anomalies) ## +3.7°C

##=====================================================================
##
## Amout of carbon loss that would be associated to future
## deforestation on the period 2010-2080
##
##=====================================================================

## Amout of carbon loss due to deforestation from 2010 to 2080
TotACD.10 <- TotACD.10
theta <- 0.5 ## %.yr-1
y <- 2080-2010
theta.prim <- 1-(1-theta/100)^y
theta.prim ## 
loss.defo <- TotACD.10*theta.prim
loss.defo ## 258,372,713 ; 258 Tg (1 Tg= 10^12 g)

## Save workspace
save.image(file="results/Results.RData")
##load("results/Results.RData")

##====================================================================================================
##============================================ END ===================================================
