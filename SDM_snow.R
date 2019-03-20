################################################################
###               Parallel SDM                               ###
################################################################
# author: yjyan

----UTF8 -----
library(biomod2)
library(raster)
library(rgdal)
library(sp)

outputpath1 <- "" 

setwd(outputpath1)
datapath <- paste(outputpath1,"Data/",sep="")
envpath <- paste(datapath,"Climate/",sep="")

bio_c <- list.files(paste(envpath,"Current/",sep=""),full.names = T) 
env <- lapply(bio_c,FUN=raster)
#bouQTP <- env[[1]]
env <- stack(env)

GPRespCoor <- read.csv(paste(datapath,"LAllgpcoor.csv",sep=""), header=T,stringsAsFactors=FALSE)
GPRespCoor <- GPRespCoor[,-1]
SDMSp <- read.table(paste(datapath,"recalsp.txt",sep=""), header=F,stringsAsFactors=FALSE)
spnames <- SDMSp[,1]
sp.n <- spnames[1]
mySpeciesOcc <- GPRespCoor
MyBiomodSF <- function(sp.n,mySpeciesOcc,myExpl,noRunList){
  ### definition of data for this run
  myRespName = sp.n
  cat('\n',myRespName,'modelling...')
  n <- which(spnames[]==myRespName) 
  myResp <- which(mySpeciesOcc$SpID==sp.n)
  # myResp <- as.numeric(mySpeciesOcc[,myRespName])
  # na.id <- which(is.na(myResp))
  # myResp <- myResp[-na.id]
  if (length(myResp)>=20){
    myRespCoord = mySpeciesOcc[myResp,c('X','Y')]
    ### Initialisation
    myRespCoord <- cbind(as.numeric(tmp_spxy[,1]),as.numeric(tmp_spxy[,2]))
    sptest <- SpatialPoints(myRespCoord,proj4string = CRS(proj4string(QTPbou)))
    spbou <- sptest[QTPbou]
    myRespCoord <- spbou@coords
    myResp <- rep(1,times=nrow(myRespCoord))
    if (length(myResp)>=20){
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.xy = myRespCoord,
                                           resp.name = myRespName,
                                           PA.nb.rep = 3,
                                           PA.nb.absences = length(myResp),
                                           PA.strategy = 'random')
      ### Modelling
      myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                          models = c('GLM','MAXENT.Phillips','RF'),
                                          #models.options = myBiomodOption,
                                          NbRunEval=3,
                                          DataSplit=85,
                                          Yweights=NULL,
                                          VarImport=3,
                                          models.eval.meth = c('TSS','ROC'),
                                          SaveObj = TRUE,
                                          rescal.all.models = T,
                                          do.full.models=F,
                                          modeling.id="test")
      ### save models evaluation scores and variables importance on hard drive
      eva <- t(as.matrix(data.frame(myBiomodModelOut@models.evaluation@val)))
      var_imp <- as.matrix(data.frame(myBiomodModelOut@variables.importances@val))
      var_mean <- data.frame(apply(var_imp,1,mean))
      capture.output(eva,
                     file=file.path(paste("ModelProperties/",myRespName,"_models_evaluation.txt", sep="")))
      capture.output(var_mean,
                     file=file.path(paste("ModelProperties/",myRespName,"_variables_importance.txt", sep="")))
      # subclim <- order(var_mean[,1], decreasing=T)
      # subclim <- subclim[1:5]
      # myExpl <- unstack(myExpl)
      # myExpl_sub <- stack(myExpl[subclim])
      evaname <- rownames(eva)
      n_AUC <- evaname[grep("Testing",evaname)]
      model_AUC <- eva[n_AUC,2]
      if (sum(model_AUC>0.9)!= 0){
        ### Building ensemble-models
        myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                              chosen.models = 'all',
                                              em.by = "all",
                                              eval.metric = c('ROC'),
                                              eval.metric.quality.threshold = c(0.9),
                                              prob.mean = F,
                                              prob.cv = F,
                                              prob.ci = F,
                                              #prob.ci.alpha = 0.05,
                                              prob.median = F,
                                              committee.averaging = F,
                                              prob.mean.weight = T)
        #prob.mean.weight.decay = 'proportional' )
        ### Do projections on current varaiable
        myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = get_needed_models(myBiomodEM),
                                          binary.meth= 'ROC',
                                          clamping.mask = F)
        ### Do ensemble-models projections on current varaiable
        myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                                 EM.output = myBiomodEM,
                                                 binary.meth = 'ROC',
                                                 total.consensus = TRUE)
      }
    }
  }
  capture.output(noRunList,
                 file=file.path(paste("ModelProperties/",sp.n,".txt", sep="")))
}

#### Modeling current distribution ####
GPnoRunList <- c(rep("1",times=length(spnames)))
library(snowfall)
sfInit(parallel=TRUE, cpus=8) 
## Export packages
sfLibrary('biomod2', character.only=TRUE)
## Export variables
sfExport('GPRespCoor')
sfExport('GPnoRunList')
sfExport('env')
sfExport('spnames')
## Do the run
sfLapply(spnames,MyBiomodSF,GPRespCoor,env,GPnoRunList)
## stop snowfall
sfStop( nostop=FALSE)
MyBiomodSF()

#### Modeling future suitability ####
setwd("")
myBiomodEMList <- list()
spnamesb <- list.files(getwd())
spnamesb <- spnamesb[5:length(spnamesb)]

for (sp.n in spnamesb){
  n <- which(spnamesb[]==sp.n)
  listOut <- list.files(paste(sp.n,"/",sep=""),full.names=T)
  emOut <- listOut[grepl("ensemble",listOut)==T]
  if (identical(emOut,character(0))==F){
    emModelOut <- load(emOut)
    myBiomodEMList[[n]] <- get(emModelOut)
    rm(list = c(emModelOut, 'emModelOut'))
  }
}

MyBiomodSMF <- function(sp.n,Explf){
  #futuretime <- c("5026")
  futuretime <- c("5026","5085","7026","7085")
  #modelname <- "average"
  #model.n <- "average"
  n <- which(spnamesb[]==sp.n)
  myBiomodEM <- myBiomodEMList[[n]]
  #return (myBiomodEM)
  if(is.null(myBiomodEM)==F){
    for (time.n in futuretime){
        projname <- paste(model.n,time.n,sep="")
        #i <- which(modelname[]==model.n)
        i <- which(futuretime[]==time.n)
        myExplf <- Explf[[i]]
        BiomodEF <- BIOMOD_EnsembleForecasting(
          EM.output = myBiomodEM,
          new.env = myExplf,
          proj.name = projname,
          binary.meth = 'ROC',
          total.consensus = TRUE)
    } 
  }
}

modelname <- c("bc","he","ip","mg","no")

library(snowfall)
for (time.n in futuretime){
 #time.n = futuretime[4]
   cat("\n",time.n)
   for (model.n in modelname){
     cat("\n",model.n)
     pathf <- paste(envpath,time.n,"/",model.n,"/",sep="")
     projname <- paste(model.n,time.n,sep="")
     Listf <- list.files(pathf,full.names=T)
     environ <- lapply(Listf,raster)
     myExplf <- stack(environ)
     #MyBiomodSMF(spnames[1])
     sfInit(parallel=TRUE, cpus=3) ## we select 2 CPUs. If you have 8 CPUs, put 8.
     sfLibrary('biomod2', character.only=TRUE)
     sfExport('spnames')
     sfExport('projname')
     sfExport('myExplf')
     sfExport('myBiomodEMList')
     sfLapply(spnames, MyBiomodSMF)
     sfStop(nostop=FALSE)
   }
 } 
