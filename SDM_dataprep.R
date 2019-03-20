################################################################
###               SDM Preparation                            ###
###        Environmental data: crop, resampling              ###
################################################################
# author: yjyan


----UTF8 -----
library(raster)
library(rgdal)
library(sp)
library(spatial)

root <- ""
setwd(root)
datapath <- paste(root, "data/", sep="") 
outputpath <- paste(root,"GPRespCoor/",sep="")

GPRespCoor <- read.csv(paste(datapath,"GPRespCoor.csv",sep=""), header=T,stringsAsFactors=FALSE)

QTPbou <- shapefile(paste(datapath,"/QTPbou/DBATP_P_Proj.shp",sep=""))

cropQTP <- function(r){
  cr <- crop(r,QTPbou)
  fr <- rasterize(QTPbou,cr)   
  r <- mask(x=cr, mask=fr)
}

#crop, project and resample current climate (temperature divided by 10)
e <- extent(c(72, 110, 23, 44)) 
bio_c <- list.files("Climate/Current/bio_5m/",full.names = T) #current bio
bio_r <- lapply(bio_c[1:19],FUN=raster)
bio_r <- lapply(bio_r,crop,e)
stand <- raster("bio1proj")

pro_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
name = c("01","02","03","04","05","06","07","08","09","10","11", "12", "13", "14", "15", "16", "17", "18", "19")
resamp_c <- list()
for (i in 1:19){
  proj4string(bio_r[[i]])<- pro_wgs 
  resamp_c[[i]] <- projectRaster(from=bio_r[[i]] , to=stand, method="bilinear")
  writeRaster(resamp_c[[i]],file=paste("Climate/Current/bio/bio",name[i],".asc",sep=""),format="ascii")
}

#temperature data, divide by 10
bio_c <- list.files("Climate/Current/bio/",full.names = T) #current bio
bio_t <- c(bio_c[1:2],bio_c[5:11])
bio_r <- lapply(bio_t,FUN=raster)
plot(bio_r[[6]])
div <- function(x)x/10
bio_t10 <- lapply(bio_r,div)
for (i in 1:9){
  writeRaster(bio_t10[[i]],file=bio_t[[i]],overwrite=T)
}

# Future climate data
bio_f <- list.files("Climate/Future/bio5min/",full.names = T) #current bio
bio_fname <- list.files("Climate/Future/bio5min/",full.names = F)
bio_rf <- lapply(bio_f,FUN=raster)
bio_rf <- lapply(bio_rf,crop,e)

pro_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
for (i in 1:180){
  bio_fname[i] <- unlist(strsplit(bio_fname[i],".tif"))
}

resamp_c <- list()
for (i in 1:180){
  proj4string(bio_rf[[i]])<- pro_wgs 
  resamp_c[[i]] <- projectRaster(from=bio_rf[[i]] , to=stand, method="bilinear")
}

#temperature data, divide by 10
bio_fr <- list.files("Climate/Future/bio10km/",full.names = T) 
bio_ft <- c(bio_fr[1:2],bio_c[5:11])
bio_r <- lapply(bio_t,FUN=raster)
plot(bio_r[[6]])
div <- function(x)x/10
bio_t10 <- lapply(bio_r,div)
bio_ft10 <- resamp_c
for (i in 1:20){
  resamp_c[[(i-1)*9+1]] <- resamp_c[[(i-1)*9+1]]/10
  resamp_c[[(i-1)*9+2]] <- resamp_c[[(i-1)*9+2]]/10
  resamp_c[[(i-1)*9+3]] <- resamp_c[[(i-1)*9+3]]/10
  resamp_c[[(i-1)*9+8]] <- resamp_c[[(i-1)*9+8]]/10
  #writeRaster(resamp_c[[i]],file=paste("Climate/Future/bio10km/",bio_fname[i],".asc",sep=""),format="ascii")
}

# Correlation of climate data
env <- stack(lapply(bio_c,FUN=raster))
allenvm <- rasterToPoints(envc)
cor_env <- cor(allenvm,method="spearman")
write.csv(cor_env,file="cor_env.csv")
