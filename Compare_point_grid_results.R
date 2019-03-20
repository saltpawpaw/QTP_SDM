### Boyce analysis ###
### Calculate the range overlap ###
# author: yjyan

----UTF8 -----

library(ecospat)
library(raster)

### 1. Caculate Boyce Index ####
root <- ""
outputpath <- paste(root,"/GPRespCoor/",sep="")
boyce_alsp <- splist993$SpList
boyce_uncalsp <- boyce_uncalsp$X1
LAllgpcoor <- LAllgpcoor_final

boyceindex <- function(sp.n){
  name <- sp.n
  if (sp.n %in% LAllgpcoor$SpID){
    tmp_spdisc <- raster(paste(outputpath,sp.n,"/proj_current/proj_current_",
                               sp.n,"_ensemble.grd",
                               sep=""),band=1)
    tmp_Resp <- which(LAllgpcoor$SpID==sp.n)
    tmp_RespCoord = LAllgpcoor[tmp_Resp,c('X','Y')]
   # get NAs id
    #na.id <- which(is.na(tmp_Resp))
    # remove NAs to enforce PA sampling to be done on explanatory rasters
    #tmp_Resp <- tmp_Resp[-na.id]
    #tmp_RespCoord = LAllgpcoor[-na.id,c('X','Y')]
    tmp_boyce <- ecospat.boyce(tmp_spdisc,tmp_RespCoord,nclass=0)[[2]]
  }else{
    tmp_boyce <- NA
  }
  return (list(spname = name, boyce.value = tmp_boyce))
}
boyce_list1 <- lapply(boyce_alsp,boyceindex)
boyceindex(boyce_alsp[218])
boyce1 <- list()
for (i in 1:993){
  boyce1[i] <- boyce_list1[[i]]$boyce.value
}
boyce1_dataframe <- cbind(boyce_alsp,do.call(rbind,boyce1))
write.csv(boyce1_dataframe,file="boyce1.csv")

#### 2. Plot point distribution & county distribution ####
PRespCoor <- read_excel("splist.xlsx", 
                        sheet = "XY_pcoor_alb")
PRespCoor[,10] <- paste("X",PRespCoor$SpID,sep="") 

plot.pCoor_gCoor <- function(sp.n){
  cat("\n",sp.n)
  tmp_Resp <- as.numeric(LAllgcoor[,sp.n])
  na.id <- which(is.na(tmp_Resp))
  tmp_Resp <- tmp_Resp[-na.id]
  tmp_RespCoord <- LAllgcoor[-na.id,c('X','Y')]
  if (length(which(PRespCoor[,10]==sp.n))!=0){
    tmp_pcor <- PRespCoor[which(PRespCoor[,10]==sp.n),8:9]
    png(filename = paste("PRespCoor_plot/",sp.n,".jpg",sep=""),width = 900,height=600)
    plot(QTPbou,main = nrow(tmp_pcor))
    points(tmp_RespCoord$X,tmp_RespCoord$Y,col="#FFC0CB7F",pch=15) ##grid based distribution
    points(tmp_pcor$X_alb,tmp_pcor$Y_alb,col="green",pch=20)  ##points based distribution
    dev.off()
  }
  # }else{
  #   png(filename = paste("PRespCoor_plot/",sp.n,".jpg",sep=""),width = 900,height=600)
  #   plot(QTPbou,main = "0")
  #   points(tmp_RespCoord$X,tmp_RespCoord$Y,col="#FFC0CB7F",pch=15) ##grid based distribution
  #   dev.off()
  # }
}

lapply(colnames(GPRespCoor1)[4:1498],plot.pCoor_gCoor)

#### 3. Calculate spatial congruence ####
library(adehabitatHR)
library(rgeos)
splist_point_sub <- subset.data.frame(splist_point,RecordNum >= 5, select = c("SpID","X_alb","Y_alb"))
colnames(splist_point_sub) <- c("SpID","X","Y")
splistp <- unique(splist_point_sub$SpID)

calPgCon <- function(sp.n){
  gcoorlist <- LAllgcoor[which(LAllgcoor$SpID==sp.n),3:4] #get grid coordinate
  if (nrow(gcoorlist)>5){
    pcoorlist <- splist_point_sub[which(splist_point_sub$SpID==sp.n),2:3] #get point coordinate
    pgcoorlist <- rbind(gcoorlist,pcoorlist)
    pcoorsp <- SpatialPoints(pcoorlist)
    gcoorsp <- SpatialPoints(gcoorlist)
    pgcoorsp <- SpatialPoints(pgcoorlist)
    prange <- mcp(pcoorsp,percent=95,unin = "km",unout = "km2")
    grange <- mcp(gcoorsp,percent=95,unin = "km",unout = "km2")
    pgrange <- mcp(pgcoorsp, percent=95,unin = "km",unout = "km2")
    pgrange_area <- pgrange@data$area
    pgintersect <- intersect(prange,grange)
    if (is.null(pgintersect)!=T){
      pgintersect_area <- pgintersect@polygons[[1]]@area
    }else{
      pgintersect_area <- 0
    }
    pgcon <- c(sp.n,pgintersect_area,pgrange_area,pgintersect_area/pgrange_area)
    return(pgcon)
  }
}
splistp_pgcon <- lapply(splistp,calPgCon)
splistp_pgcon_db <- do.call(rbind,splistp_pgcon)
colnames(splistp_pgcon_db) <- c("spID","pgintersect_area","pgrange_area",
                                "pgcon")
write.csv(splistp_pgcon_db,file="splistp_pgcon.csv")

#### 4. Comparing grid&point modeled results ####
changearea_g <- cbind(paste(changearea_g[,1],changearea_g[,2],sep="_"),
                      changearea_g)
names(changearea_g)[1] <- "ID"
changearea_p <- cbind(paste(changearea_p[,1],changearea_p[,2],sep="_"),
                      changearea_p)
names(changearea_p)[1] <- "ID"
#combine grid & point data frame
changearea <- merge(changearea_g,changearea_p,by="ID",all=F)
plot(changearea$CurrentRangeSize.x,changearea$CurrentRangeSize.y)
write.csv(changearea,file="changearea_cp_gd.csv")
changearea_per <- changearea[,16:25]/changearea[,4:13]

#### 5.Calculate niche overlap ####
library(phyloclim)
# Convert raster to grid
calNicheOverP <- function(sp.n){
  spdisc_g <- cropQTP(raster(paste(outputpath,"model_g/",sp.n,"/proj_current/proj_current_",
                                   sp.n,"_ensemble.grd",
                                   sep=""),band=1))
  spdisc_p <- cropQTP(raster(paste(outputpath,"model_p/",sp.n,"/proj_current/proj_current_",
                                   sp.n,"_ensemble.grd",
                                   sep=""),band=1))
  writeRaster(spdisc_g,file=paste(outputpath,"currentprob/",sp.n,"_g.asc",sep=""),overwrite=T)
  writeRaster(spdisc_p,file=paste(outputpath,"currentprob/",sp.n,"_p.asc",sep=""),overwrite=T)
  tmp_niche_over <- niche.overlap(c(paste(sp.n,"_g.asc",sep=""),paste(sp.n,"_p.asc",sep="")))
  tmp_niche <- c(tmp_niche_over[2,1],tmp_niche_over[1,2])#first is I, seconde is D
  return(tmp_niche)
}
nicheover_list <- lapply(splist,calNicheOverP)
nicheoverp_db <- do.call(rbind,nicheover_list)
write.csv(nicheoverp_db,file="nicheoverp_db.csv") 
