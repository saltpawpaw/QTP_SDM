################################################################
###               SDM result analysis                        ###
### Richness pattern, species num change, elevation change   ###
### Species range size change, summarize model performance   ###
################################################################
# author: yjyan

----UTF8 -----

library(biomod2)
library(maptools)
library(raster)
library(rgdal)
library(sp)
library(spatial)
library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)

# input: species list, path to SDM result
# input: DEM of study area, study time period/model name
data_dir <- ""
setwd(data_dir)
result_dir <- paste(data_dir,"result",sep="")
if (file.exists(result_dir)==F){
  dir.create(file.path(data_dir,result_dir))
}
library(readr)
spinfo <- read_delim("splist.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
splist <- spinfo$SpList

### 1. Model performace ####
Extract_Eva <- function (x) {  
  test <- data.frame(matrix(nrow = 27, ncol = 6))
  evalname <- rownames(x)
  mod_eval_names <- c("name")
  #i=1
  for (i in 1: (nrow(x)/4)){                        
    test[i,2:3] <- x[1+4*(i-1),1:2]
    mod_eval_names[i] <- evalname[1+4*(i-1)]
  }
  #test <- data.frame(test)
  evalname1 <- data.frame(matrix(ncol = 5, nrow = nrow(test)))
  for (i in 1:nrow(test)){
    evalname1[i,] <- unlist(strsplit(mod_eval_names[i],"[.]"))
  }
  test[,1] <- evalname1[,3]
  test[,4:5] <- evalname1[,4:5]
  colnames(test) <- c('Model', 'TSS', 'ROC','RUN','PA',"SpID")
  return(test)
}

mMpList <- list.files(paste(data_dir,"ModelProperties/",sep=""),full.names = T)
mEvaList <- mMpList[grepl("models",mMpList)==T] 

AllEva <- list()
for (sp.n in splist){
  print(sp.n)
  ROC_sp <- matrix(ncol=1 ,nrow=27)
  n <- which(splist[]==sp.n)
  eva_sp <- data.frame(read.table(mEvaList[[n]],sep="",row.names=1,stringsAsFactors=F,header=T))
  model_eva <- Extract_Eva(eva_sp)
  #model_eva[,6] <- rep(sp.n,times=27)
  model_eva <- melt(model_eva[,1:3], id.vars = "Model",      
                    measure.vars = c("TSS", "ROC"))
  model_eva <- cbind(rep(sp.n,times=length(model_eva)),model_eva)
  colnames(model_eva) <- c("SpID","Model", "Eval", "Value")
  AllEva[[n]] <- model_eva
}

AllEva_df <- AllEva[[1]]
for (i in 2:length(AllEva)){
  AllEva_df <- rbind(AllEva_df,AllEva[[i]])
}
names(AllEva_df)[1] <- "spNo"
names(spinfo)[1] <- "spNo"
AllEva_df[,1] <- as.character(AllEva_df[,1])
AllEvasp <- inner_join(AllEva_df,SDMSP,by="spNo") #add species information
write.csv(AllEvasp,file="result/Evaluation_Metrics_AllSp.csv")

### 2. Factor contribution ####
mVimList <- mMpList[grepl("importance",mMpList)==T] 

varimps_GPCX <- matrix(ncol=length(splist),nrow= 8)
for (sp.n in splist){
  n <- which(splist[]==sp.n)
  varimp <- data.frame(read.table(paste(data_dir,"ModelProperties/",splist[n],"_variables_importance.txt",sep=""),sep="",stringsAsFactors=FALSE,header=T))
  varimps_GPCX[,n] <- varimp[,1]
}
varimps_GPCX <- t(varimps_GPCX)
climname <- rownames(varimp)
colnames(varimps_GPCX)<-climname
rownames(varimps_GPCX)<- splist
varimps_GPCX <- cbind(spinfo,as.data.frame(varimps_GPCX))
write.csv(varimps_GPCX,file=paste(result_dir,"/factor_contribution.csv",sep=""))

# Plot boxplot for different factors 
oldpar <- par()
png(filename="result/varimps-1.png",width=1200,height=900,bg="white")
par(mfrow = c(1,1),cex.axis = 3.5, cex.lab=5,mar = c(15,10,5,3))#,mar=c(3,3,3,3)
boxplot(varimps_GPCX[,7:14],las=2,xaxt="n")
axis(1,at=1:8,labels=F)
#climname <- c("MAT","MDR","TSN","MTWQ","MTCQ","AP","PSN","PWEQ")
text(1:8, par("usr")[3]-0.035, srt = 0, 
     labels = climname, xpd = TRUE,cex=3)
dev.off()

##Get the factor that contribute most to each sp
maxBio <- list()
for (i in 1:nrow(varimps_GPCX)){
  maxBio[i] <- colnames(sort(varimps_GPCX[i,7:14],decreasing = T))[1]
}
maxBio <- do.call(rbind,maxBio)
varimps_GPCX <- cbind(varimps_GPCX,maxBio)
write.csv(varimps_GPCX,file="varimps_GPCX.csv")
#varimps <- melt(varimps_GPCX[,2:9], id.vars = "splist",      
#                measure.vars = climname)
#varimps[,3]<- as.numeric(levels(varimps[ ,3]))[varimps[ ,3]]
#colnames(varimps) <- c("SpID","Bio","Imp")

### Table 1.range size change & conservation status change ####
splist_r <- spinfo[which(spinfo$SpRecal==4),1]
root <- ""
changeareaAll <- matrix(nrow=length(splist_r)*4,ncol=10)
for (sp.n in splist_r){
  #sp.n <- splist[1]
  #time.n <- futuretime[2]
  n <- which(splist_r[]==sp.n)
  spdisc <- raster(paste(data_dir,sp.n,"/proj_current/proj_current_",
                         sp.n,"_ensemble_ROCbin.grd",
                         sep=""),band=1)
  for (time.n in futuretime){
    #time.n <- futuretime[1]
    layername <- paste("/proj_",time.n,sep="")
    cat('\n',layername)
    m <- which(futuretime[]==time.n)
    spdisf <- raster(paste(data_dir,sp.n,layername,layername,"_",
                           sp.n,"_ensemble_ROCbin.grd",
                           sep=""),band=1)
    changecf <- BIOMOD_RangeSize(spdisc,spdisf,SpChange.Save=paste("change",sp.n,time.n,sep=""))
    changeareaAll[4*(n-1)+m,] <- changecf$Compt.By.Models
    rownamechange[4*(n-1)+m] <- paste(sp.n,time.n,sep="_")   
    writeRaster(changecf$Diff.By.Pixel[[1]],paste(root,"binarychangeAll/",sp.n,"_average",time.n,".asc",sep=""),format="ascii",overwrite=T)
  }
}
rownames(changeareaAll) <- rownamechange
colnames(changeareaAll) <- colnames(changecf$Compt.By.Models)
write.csv(changeareaAll,file="result/changeareaAll_rsp.csv")

changeareaAll <- read.csv(file="result/changeareaAll_rsp.csv",header=T,stringsAsFactors = F)
changeareaAll <- changeareaAll[,-3]

cutC <- c(0,30,50,80,100)
lossB <- cut(changearea_s$PercLoss,breaks=cutC,labels=c("-30 - 0","-50 - -30","-80 - -50","-100 - -80"))
cutC_m <- c(-100,-80,-50,-30,0,600)
lossB_fulldis <- cut(changearea_s$SpeciesRangeChange,breaks=cutC_m,labels=c("-100 - -80","-80 - -50","-50 - -30","-30 - 0","> 0"))
changearea_s <- cbind(changearea_s,lossB,lossB_fulldis)
write.csv(changearea_s,file="result/changeareaAll.csv")


### Fig1. Range center change ####
WeightedCenter <- function(dat,polys){ 
  proj4string(dat)  <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  dat <- crop(dat,extent(polys))
  ## Compute weighted x and y coordinates within each rasterized region
  z <- rasterize(polys, dat)## Convert polygons to a raster layer
  xx <- zonal(init(dat, v="x")*dat, z) / zonal(dat,z)
  yy <- zonal(init(dat, v="y")*dat, z) / zonal(dat,z)
  centerXY <- c(xx[1,2],yy[1,2])
  return (centerXY)
}

spProCenter <- function(spname,time.n){
  layername <- paste("/proj_",time.n,sep="")
  dat <- raster(paste(data_dir,spname,layername,layername,"_",
                      spname,"_ensemble.grd",
                      sep=""),band=1)
  return(WeightedCenter(dat,QTPbou))
}

spProCenterXY2000 <- lapply(splist,spProCenter,times[1])
spProCenterXY20001 <- data.frame(t(as.data.frame(spProCenterXY2000)))
rownames(spProCenterXY20001) <- splist
spProCenterXY20001 <- cbind(spname=rownames(spProCenterXY20001),spProCenterXY20001)

spProCenterXY5026_1 <- data.frame(t(as.data.frame(lapply(splist_o,spProCenter,times[2]))))
spProCenterXY5026_2 <- data.frame(t(as.data.frame(lapply(splist_r,spProCenter,futuretime[1]))))
rownames(spProCenterXY5026_1) <- splist_o
rownames(spProCenterXY5026_2) <- splist_r
spProCenterXY5026_S <- rbind(spProCenterXY5026_1,spProCenterXY5026_2)
spProCenterXY5026_S <- cbind(spname = rownames(spProCenterXY5026_S),spProCenterXY5026_S)

spProCenter_S <- inner_join(spProCenterXY20001,spProCenterXY5026_S,by="spname")

spProCenterXY5085_1 <- data.frame(t(as.data.frame(lapply(splist_o,spProCenter,times[3]))))
spProCenterXY5085_2 <- data.frame(t(as.data.frame(lapply(splist_r,spProCenter,futuretime[2]))))
rownames(spProCenterXY5085_1) <- splist_o
rownames(spProCenterXY5085_2) <- splist_r
spProCenterXY5085_S <- rbind(spProCenterXY5085_1,spProCenterXY5085_2)
spProCenterXY5085_S <- cbind(spname = rownames(spProCenterXY5085_S),spProCenterXY5085_S)
spProCenter_S <- inner_join(spProCenter_S,spProCenterXY5085_S,by="spname")

spProCenterXY7026_1 <- data.frame(t(as.data.frame(lapply(splist_o,spProCenter,times[4]))))
spProCenterXY7026_2 <- data.frame(t(as.data.frame(lapply(splist_r,spProCenter,futuretime[3]))))
rownames(spProCenterXY7026_1) <- splist_o
rownames(spProCenterXY7026_2) <- splist_r
spProCenterXY7026_S <- rbind(spProCenterXY7026_1,spProCenterXY7026_2)
spProCenterXY7026_S <- cbind(spname = rownames(spProCenterXY7026_S),spProCenterXY7026_S)
spProCenter_S <- inner_join(spProCenter_S,spProCenterXY7026_S,by="spname")

spProCenterXY7085_1 <- data.frame(t(as.data.frame(lapply(splist_o,spProCenter,times[5]))))
spProCenterXY7085_2 <- data.frame(t(as.data.frame(lapply(splist_r,spProCenter,futuretime[4]))))
rownames(spProCenterXY7085_1) <- splist_o
rownames(spProCenterXY7085_2) <- splist_r
spProCenterXY7085_S <- rbind(spProCenterXY7085_1,spProCenterXY7085_2)
spProCenterXY7085_S <- cbind(spname = rownames(spProCenterXY7085_S),spProCenterXY7085_S)
spProCenter_S <- inner_join(spProCenter_S,spProCenterXY7085_S,by="spname")

ls(spProCenter_S)
spProCenter_S <- spProCenter_S[,-2]
colnames(spProCenter_S) <- c("spname","X_current","Y_current",
                             "X_5026","Y_5026","X_5085","Y_5085",
                             "X_7026","Y_7026","X_7085","Y_7085")

# angle change
calAngle <- function(i,data1,data2){ # data1起始点, data2终点坐标
  x1 = data1[i,1]
  x2 = data2[i,1]
  y1 = data1[i,2]
  y2 = data2[i,2]
  z = atan((abs(y2-y1))/(abs(x2-x1)))/(2*pi)*360
  if (y2>y1&&x2>x1){
    zn = 90-z
  }else if (y2<y1&&x2>x1){
    zn = 90+z
  }else if (y2<y1&&x2<x1){
    zn = 270-z
  }else if (y2>y1&&x2<x1){
    zn = 270+z
  }
  zn
}
anglef <- matrix(ncol=4,nrow=nrow(spProCenter_S))
for (j in 1:nrow(spProCenter_S)){
  anglef[j,1] <- calAngle(j,spProCenter_S[,2:3],spProCenter_S[,4:5])
  anglef[j,2] <- calAngle(j,spProCenter_S[,2:3],spProCenter_S[,6:7])
  anglef[j,3] <- calAngle(j,spProCenter_S[,2:3],spProCenter_S[,8:9])
  anglef[j,4] <- calAngle(j,spProCenter_S[,2:3],spProCenter_S[,10:11])
}
#dir <- cut_interval(c(0,anglef[,1],360),n=16)

##dispersal distance
colnames(spProCenter_S) <- c("spname","X","Y",
                             "X","Y","X","Y",
                             "X","Y","X","Y")

spProDist26 <- dist(rbind(spProCenter_S[,2:3],spProCenter_S[,4:5],spProCenter_S[,8:9]),method="euclidean") 
spProDist85 <- dist(rbind(spProCenter_S[,2:3],spProCenter_S[,6:7],spProCenter_S[,10:11]),method="euclidean")
spProDist26 <- as.matrix(spProDist26)
spProDist85 <- as.matrix(spProDist85)
spProDist851 <- spProDist85[(length(splist)+1):(length(splist)*2),1:length(splist)]
spProDist852 <- spProDist85[(length(splist)*2+1):(length(splist)*3),1:length(splist)]
spProDist851 <- diag(spProDist851)
spProDist852 <- diag(spProDist852)

spProDist261 <- spProDist26[(length(splist)+1):(length(splist)*2),1:length(splist)]
spProDist262 <- spProDist26[(length(splist)*2+1):(length(splist)*3),1:length(splist)]
spProDist261 <- diag(spProDist261)
spProDist262 <- diag(spProDist262)
spProDist <- data.frame(Dist=c(spProDist261,spProDist262,spProDist851,spProDist852),
                        Year=rep(c("2050","2070","2050","2070"),each=length(splist)),
                        RCP=rep(c("RCP 2.6","RCP 8.5"),each=length(splist)*2))
write.csv(spProDist,file="result/spProDist.csv")
# Summarize
cutA <- seq(from=0,to=360,by=22.5)
dirA <- list(cut(anglef[,1],breaks=cutA,labels=colnames(windfr)),
          cut(anglef[,2],breaks=cutA,labels=colnames(windfr)),
          cut(anglef[,3],breaks=cutA,labels=colnames(windfr)),
          cut(anglef[,4],breaks=cutA,labels=colnames(windfr)))
library(climatol)
dirA <- dirA[[1]]
distA <- distA5026
CompWindfr <- function(dirA,distA,splist){
  tmp_windfr <- data.frame(cbind(dirA,distA,rep(1,times=length(splist))))
  tmp_windfr5 <- dcast(tmp_windfr,distA~dirA,sum)
  colnames(tmp_windfr5) <- c("N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW" 
                         ,"WSW","W","WNW","NW","NNW")
  rownames(tmp_windfr5) <- levels(distA)
  return (tmp_windfr5)
}
cutD1 <- c(0,50,100,150,305)
distA5026 <- cut(spProDist261,breaks=cutD1,labels=c("<50","50 - 100","100-150",">150"))
cutD1 <- c(0,50,100,150,400)
distA5085 <- cut(spProDist851,breaks=cutD1,labels=c("<50","50 - 100","100-150",">150"))
cutD1 <- c(0,50,100,150,300)
distA7026 <- cut(spProDist262,breaks=cutD1,labels=c("<50","50 - 100","100-150",">150"))
cutD1 <- c(0,100,200,300,500)
distA7085 <- cut(spProDist852,breaks=cutD1,labels=c("<100","100 - 200","200-300",">300"))

windfr5026 <- CompWindfr(dirA[[1]],distA5026,splist)
windfr5085 <- CompWindfr(dirA[[2]],distA5085,splist)
windfr7026 <- CompWindfr(dirA[[3]],distA7026,splist)
windfr7085 <- CompWindfr(dirA[[4]],distA7085,splist)
write.csv(rbind(windfr5085,windfr7026,windfr7085),
          file="result/windfr.csv")

#### Fig.2 Protection ratio evaluation ####
#recognize whether the species is in a conserve area or not
con_areag <- shapefile(paste(datapath,"reservebou/Reserve_Bou_Alb.shp",sep=""))## national nature reserve
proj4string(con_areag) <- CRS("+proj=aea +lat_1=25 +lat_2=47 +lat_0=30 +lon_0=110 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

sumCon <- function(spname,time.n,conarea){ #conarea为保护区的polygon
  layername <- paste("/proj_",time.n,sep="")
  if (file.exists(paste(data_dir,spname,layername,layername,"_",
                        spname,"_ensemble_ROCbin.grd",
                        sep=""),band=1) == TRUE){
    dat <- raster(paste(data_dir,spname,layername,layername,"_",
                        spname,"_ensemble_ROCbin.grd",
                        sep=""),band=1)
    dat <- cropQTP(dat)
    dat[is.na(dat[])] <- 0 
    extra_rst <- extract(dat,conarea,fun=sum)#sum the value
  }
}
times_r <- c(times[1],futuretime)
library(snowfall)
sfInit(parallel=TRUE, cpus=3)
sfLibrary('raster', character.only=TRUE)
sfExport('splist')
sfExport('sumCon')
sfExport('times_r')
sfExport('QTPbou')
sfExport('cropQTP')
sfExport('con_areag')
sfExport('data_dir')

spconser2000 <- sfLapply(splist,sumCon,times_r[[1]],con_areag)
spconser5026 <- sfLapply(splist,sumCon,times_r[[2]],con_areag)
spconser5085 <- sfLapply(splist,sumCon,times_r[[3]],con_areag)
spconser7026 <- sfLapply(splist,sumCon,times_r[[4]],con_areag)
spconser7085 <- sfLapply(splist,sumCon,times_r[[5]],con_areag)

spconser5026 <- spconser5026[!sapply(spconser5026, is.null)] 
spconser5085 <- spconser5085[!sapply(spconser5085, is.null)] 
spconser7026 <- spconser7026[!sapply(spconser7026, is.null)] 
spconser7085 <- spconser7085[!sapply(spconser7085, is.null)] 

spconser20001 <- rapply(spconser2000,c)
write.csv(spconser20001,"result/spconser20001.csv")
spconser50261 <- rapply(spconser5026,c)
write.csv(spconser50261,"result/spconser50261.csv")
spconser50851 <- rapply(spconser5085,c)
write.csv(spconser50851,"result/spconser50851.csv")
spconser70261 <- rapply(spconser7026,c)
write.csv(spconser70261,"result/spconser70261.csv")
spconser70851 <- rapply(spconser7085,c)
write.csv(spconser70851,"result/spconser70851.csv")

### Analysis for each nature reserve 
sumCon2 <- function(spconser,spconser1,spnames,namec){ 
  n <- length(spconser1)/length(spconser) #n: number of reserves
  spconser <- data.frame(cbind(as.numeric(do.call(rbind,spconser)),rep(1:n,times=length(spconser)),rep(spnames,each=n)))
  spconser[,1] <- as.numeric(as.character(spconser[,1]))
  spconser[,2] <- as.numeric(as.character(spconser[,2]))
  spconser[,3] <- as.character(spconser[,3])
  spconser <- subset(spconser,spconser[,1]>0)
  #spconser <- inner_join(spconser,spnames2,by=c("X3"="spnames"))
  write.csv(spconser,file=paste(data_dir,namec,sep=""))
  return(spconser)
}

spconser20002 <- sumCon2(spconser2000,spconser20001,splist,"result/spconser20002.csv")
spconser50262 <- sumCon2(spconser5026,spconser50261,splist_r,"result/spconser50262.csv")
spconser50852 <- sumCon2(spconser5085,spconser50851,splist_r,"result/spconser50852.csv")
spconser70262 <- sumCon2(spconser7026,spconser70261,splist_r,"result/spconser70262.csv")
spconser70852 <- sumCon2(spconser7085,spconser70851,splist_r,"result/spconser70852.csv")

spconser50852_2 <- read.csv("result/spconser50852_2.csv",header=T,stringsAsFactors = F)
spconser70262_2 <- read.csv("result/spconser70262_2.csv",header=T,stringsAsFactors = F)
spconser70852_2 <- read.csv("result/spconser70852_2.csv",header=T,stringsAsFactors = F)

selectSp <- function(spconser_2,spconser){
  spconser_2 <- spconser_2[,-1]
  spconser_2[,2] <- as.numeric(spconser_2[,2])
  spconser <- rbind(spconser_2,spconser)
  spconser_s <- spconser[which(spconser[,3]%in%splist),1:3] 
  return (spconser_s)
}

spconser50262_s <- selectSp(spconser50262_2,spconser50262)
spconser50852_s <- selectSp(spconser50852_2,spconser50852)
spconser70262_s <- selectSp(spconser70262_2,spconser70262)
spconser70852_s <- selectSp(spconser70852_2,spconser70852)

write.csv(spconser50262_s,file="result/spconser50262_s.csv")
write.csv(spconser50852_s,file="result/spconser50852_s.csv")
write.csv(spconser70262_s,file="result/spconser70262_s.csv")
write.csv(spconser70852_s,file="result/spconser70852_s.csv")

confreq <- list()
confreq[[1]] <- data.frame(cbind(table(spconser20002[,2])))
confreq[[1]] <- data.frame(NNR = rownames(confreq[[1]]),confreq[[1]])
confreq[[2]] <- data.frame(cbind(table(spconser50262_s[,2])))
confreq[[2]] <- data.frame(NNR = rownames(confreq[[2]]),spNum = confreq[[2]])
confreq[[3]] <- data.frame(cbind(table(spconser50852_s[,2])))
confreq[[3]] <- data.frame(NNR = rownames(confreq[[3]]),spNum = confreq[[3]])
confreq[[4]] <- data.frame(cbind(table(spconser70262_s[,2])))
confreq[[4]] <- data.frame(NNR = rownames(confreq[[4]]),spNum = confreq[[4]])
confreq[[5]] <- data.frame(cbind(table(spconser70852_s[,2])))
confreq[[5]] <- data.frame(NNR = rownames(confreq[[5]]),spNum = confreq[[5]])

confreq_df <- left_join(confreq[[1]],confreq[[2]],by="NNR") #combine results of different period
confreq_df <- left_join(confreq_df,confreq[[3]],by="NNR")
confreq_df <- left_join(confreq_df,confreq[[4]],by="NNR")
confreq_df <- left_join(confreq_df,confreq[[5]],by="NNR")
confreq_df[,1] <- as.numeric(confreq_df[,1])
confreq_df <- confreq_df[order(confreq_df$NNR),]
con_areag_data <- con_areag@data #extract the dataframe from conservation area polygon
confreq_df <- cbind(con_areag_data,confreq_df)
write.csv(confreq_df,file=paste(data_dir,"result/confreq_df.csv",sep=""))

Changearea_1 <- Changearea[which(Changearea$SpID%in%splist),1:13]
Changearea_1[,11] <- as.numeric(as.character(Changearea_1[,11])) 
Changearea_1[,13] <- as.character(Changearea_1[,13])
changearea_s <- rbind(Changearea_1,changeareaAll)
write.csv(changearea_s,"result/changearea_s.csv")

# get current range size unduplicated
ind <- duplicated(changearea_s$SpID)
sprangesize <- cbind(changearea_s$SpID,changearea_s$CurrentRangeSize)
sprangesize <- sprangesize[!ind,]
sprangesize[,2] <- as.numeric(sprangesize[,2])
colnames(sprangesize) <- c("spNo","crangesize")
#sprangesize[,1] <- as.character(sprangesize[,1])
sprangesize <- cbind(sprangesize,changearea_s[which(changearea_s$Period=="5026"),10],changearea_s[which(changearea_s$Period=="5085"),10],
                      changearea_s[which(changearea_s$Period=="7026"),10],changearea_s[which(changearea_s$Period=="7085"),10])
colnames(sprangesize)<- c("spNo","current","5026","5085","7026","7085") 
write.csv(sprangesize,file=paste(data_dir,"sprangesize.csv",sep=""))
sprangesize <- as.data.frame(sprangesize,stringsAsFactors=F)
sprangesize[,2:6] <- sapply(sprangesize[,2:6],as.numeric)
# get conserved range for each species
spSumCon <- function(spconser){###input: spconser20002
  spsumcon <- tapply(spconser[,1],spconser[,3],sum)
  spsumcon <- data.frame(spNo=as.character(rownames(spsumcon)),conarea=spsumcon)
  spsumcon[,1] <- as.character(spsumcon[,1])
  return(spsumcon)
}

spsumcon2000 <- spSumCon(spconser20002)
spsumcon5026 <- spSumCon(spconser50262_s)
spsumcon5085 <- spSumCon(spconser50852_s)
spsumcon7026 <- spSumCon(spconser70262_s)
spsumcon7085 <- spSumCon(spconser70852_s)

spsumconp <- left_join(sprangesize,spsumcon2000,by="spNo")
spsumconp <- left_join(spsumconp,spsumcon5026,by="spNo")
spsumconp <- left_join(spsumconp,spsumcon5085,by="spNo")
spsumconp <- left_join(spsumconp,spsumcon7026,by="spNo")
spsumconp <- left_join(spsumconp,spsumcon7085,by="spNo")
colnames(spsumconp) <- c(colnames(spsumconp)[1:6],"con2000","con5026","con5085","con7026","con7085")
cprorate <- data.frame(spsumconp,pcurrent=spsumconp$con2000/spsumconp$current,
                p5026=(spsumconp$con5026)/(spsumconp$`5026`),
                p5085=spsumconp$con5085/spsumconp$`5085`,
                p7026=spsumconp$con7026/spsumconp$`7026`,
                p7085=spsumconp$con7085/spsumconp$`7085`)

# cut protection ratio into different categories
cut_ratio <- c(seq(from=0,to=0.8,by=0.1),1)
label_ratio <- c("<10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%",">80%")
cprorate_cut <- data.frame(current_cut=cbind(table(cut(cprorate$pcurrent,breaks=cut_ratio,labels=label_ratio))),
                     p5026_cut=cbind(table(cut(cprorate$p5026,breaks=cut_ratio,labels=label_ratio))),
                     p5085_cut=cbind(table(cut(cprorate$p5085,breaks=cut_ratio,labels=label_ratio))),
                     p7026_cut=cbind(table(cut(cprorate$p7026,breaks=cut_ratio,labels=label_ratio))),
                     p7085_cut=cbind(table(cut(cprorate$p7085,breaks=cut_ratio,labels=label_ratio))))
write.csv(cprorate_cut,file="result/cprorate_cut.csv")

### Calculate the protection ratio change for different lifeforms
cproratebox <- data.frame(prorate=c(cprorate$pcurrent,cprorate$p5026,cprorate$p5085,
                                    cprorate$p7026,cprorate$p7085),
                          spNo=rep(cprorate$spNo,times=5),
                          CCS=rep(cprorate$Status,times=5),
                          Lifeform=rep(cprorate$Lifeform,times=5),
                          Time=rep(c("a2000","a2650","a8550","a2670","a8570"),each=nrow(cprorate)),
                          Year=c(rep("2000",times=nrow(cprorate)),rep("2050",times=nrow(cprorate)*2),rep("2070",times=nrow(cprorate)*2)),
                          RCP = rep(c("Current","RCP 2.6","RCP 8.5","RCP 2.6","RCP 8.5"),each=nrow(cprorate)))

ccproratebox <- data.frame(prorate=c(ccprorate$currentc,ccprorate$C5026,ccprorate$C5085,
                                     ccprorate$C7026,ccprorate$C7085),
                           spNo=rep(ccprorate$spNo,times=5),
                           CCS=rep(ccprorate$CCS2,times=5),
                           Lifeform=rep(ccprorate$Lifeform,times=5),
                           Time=rep(c("a2000","a2650","a8550","a2670","a8570"),each=1095),
                           Year=c(rep("2000",times=1095),rep("2050",times=1095*2),rep("2070",times=1095*2)),
                           RCP = rep(c("Current","RCP 2.6","RCP 8.5","RCP 2.6","RCP 8.5"),each=1095))
ccproratebox_sub <- subset(ccproratebox,ccproratebox$prorate<=1) 

mycolor <- rainbow(5, alpha=0.2)
mycolor = brewer.pal(5,"Set2")

ggplot(cproratebox,aes(x=CCS,y=prorate))+ 
  geom_boxplot(aes(col=Time),size=0.3)+
  scale_x_discrete(breaks=c("EN", "LC", "DD"),
                   labels=c("Threatened", "Low Risk", "Data Deficient"))+
  ylab("Protection ratio")+
  xlab("Conservation Status")+
  facet_grid(.~Lifeform)+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.15,"cm"),
        axis.text.x = element_text(margin=unit(c(0.2,0.2,0.2,0.2),"cm"),angle=30, vjust=0.5),#legend.position="top",
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25),"cm")),
        legend.title=element_blank(),legend.text = element_text(size = 10),
        strip.text = element_text(size = 14),
        panel.grid=element_blank()
        #panel.border=element_blank(),
        #axis.line=element_line(size=1,color="black")
  )+
  scale_color_manual(values=mycolor,name = "",
                     labels = c("Current", "RCP 2.6, Y2050", "RCP 2.6, Y2070","RCP 8.5, Y2050","RCP 8.5, Y2070"))

ggsave("result/plot/cprorate1.png",width=7,height=5)
write.csv(ccproratebox,file="ccproratebox.csv")

tmp_pro <- ccproratebox[which(RCP=="Current"),]
attach(tmp_pro)
aggregate(tmp_pro,by=list(CCS),FUN=mean)

### chi-sq test
consp_en <- rbind(cprorate[which(cprorate$Status=="EN"),13:17],
                  cprorate[which(cprorate$Status=="VU"),13:17],
                  cprorate[which(cprorate$Status=="CR"),13:17])
consp_lc <- rbind(cprorate[which(cprorate$Status=="LC"),13:17],
                  cprorate[which(cprorate$Status=="NT"),13:17])
consp_w <- cprorate[which(spsumconp$Lifeform=="woody"),13:17]
consp_h <- cprorate[which(spsumconp$Lifeform=="herb"),13:17]


result_t <- t.test(consp_en$pcurrent,consp_lc$pcurrent,alternative="greater")

lapply(consp_lc[,2:5],t.test,consp_lc$pcurrent,paired=T,alternative="greater")
lapply(consp_en[,2:5],t.test,consp_en$pcurrent,paired=T,alternative="less")

lapply(cprorate[,14:17],t.test,cprorate$pcurrent,paired=T,alternative="greater")
t.test(cprorate$p5026,cprorate$p7026,alternative="greater",paired=T)
t.test(cprorate$p5085,cprorate$p7085,alternative="greater",paired=T)

lapply(consp_w[,2:5],t.test,consp_w$pcurrent,paired=T,alternative="less")
lapply(consp_h[,2:5],t.test,consp_h$pcurrent,paired=T,alternative="less")
