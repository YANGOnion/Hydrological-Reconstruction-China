
#-------------------------------------------------------------------------------#

#' @author Wencong Yang
#' @description calibration of HBV snow module

#-------------------------------------------------------------------------------#

#' @description calibration of all grids in a 1 degree tile 
tileCalib=function(n,e,path){
  require(data.table)
  require(raster)
  require(Rcpp)
  require(RcppDE)
  require(hydroGOF)
  require(tiff)
  sourceCpp('cpp/hbv_snow.cpp')
  ## MOD10 metric
  biclass=function(pred,obs){# balanced
    CP=sum(obs)
    CN=sum(!obs)
    TPR=sum(pred & obs)/CP
    TNR=sum((!pred) & (!obs))/CN
    if(CP==0){
      return(TNR)
    }else if(CN==0){
      return(TPR)
    }else{
      return((TPR+TNR)/2)
    }
  }
  ## calib on 1 grid
  gridCalib=function(i,NP=20,itermax=20){
    rowind=ind[i,1]
    colind=ind[i,2]
    dt_prcp=data.table(date=prcp_date,prcp=r_prcp[rowind,colind,])
    dt_tair=data.table(date=tair_date,tair=r_tair[rowind,colind,])
    dt_cloud=data.table(date=mod_date,cloud=r_cloud[rowind,colind,])
    dt_cover=data.table(date=mod_date,cover=r_cover[rowind,colind,])
    dt_depth=data.table(date=depth_date,depth=r_depth[rowind,colind,])
    dt=Reduce(function(a,b){merge(a,b,by='date',all=T)},list(dt_prcp,dt_tair,dt_cloud,dt_cover,dt_depth))
    dt=dt[!is.na(prcp)][!is.na(tair)] # # filter non-land grid
    if(nrow(dt)==0) return(rep(NA,8)) # no forcing data: NA
    if(dt_tair[tair<0,.N]/dt_tair[,.N]<0.05) return(c(0,0,1,5,0.1,5,NA,NA)) # no snow: default param & NA metric
    ## obs
    mod_ind=dt[367:.N][,which((month(date)%in%c(1:4,11:12))&(cloud<40))]
    obs_cover=dt[367:.N][mod_ind,cover/(100-cloud)>=0.1]
    depth_ind=dt[367:.N][,which((month(date)%in%c(1:4,11:12)))]
    obs_depth=dt[367:.N][depth_ind,depth]
    calib_ind=367:nrow(dt)
    useDepth=(!any(is.na(obs_depth)))&(sd(obs_depth)!=0)
    useCover=(length(mod_ind)!=0)
    ## calibration
    lower_par=c(TP=-3,TM=-3,SCF=0.9,CFX=0.5,CWH=0,Tcover=2)
    upper_par=c(TP=3,TM=3,SCF=1.5,CFX=10,CWH=0.2,Tcover=10)
    if(useDepth&useCover){
      hbvSnow=function(param){
        dt_snow=runModel(df=dt,I_Snowpack0=10,I_Meltwater0=10,
                         P_TP=param[1],P_TM=param[2],P_SCF=param[3],P_CFX=param[4],P_CFR=0.05,P_CWH=param[5])
        kg=KGE(dt_snow[calib_ind,]$Snowpack[depth_ind],obs_depth)
        if(is.na(kg)){
          kg=0
        }
        bi=biclass(dt_snow[calib_ind,]$Snowpack[mod_ind]>param[6],obs_cover)
        return(1-1/4*bi-3/4*kg)
      }
    }else if(useDepth){
      hbvSnow=function(param){
        dt_snow=runModel(df=dt,I_Snowpack0=10,I_Meltwater0=10,
                         P_TP=param[1],P_TM=param[2],P_SCF=param[3],P_CFX=param[4],P_CFR=0.05,P_CWH=param[5])
        kg=KGE(dt_snow[calib_ind,]$Snowpack[depth_ind],obs_depth)
        if(is.na(kg)){
          kg=0
        }
        return(1-kg)
      }
    }else if(useCover){
      hbvSnow=function(param){
        dt_snow=runModel(df=dt,I_Snowpack0=10,I_Meltwater0=10,
                         P_TP=param[1],P_TM=param[2],P_SCF=param[3],P_CFX=param[4],P_CFR=0.05,P_CWH=param[5])
        bi=biclass(dt_snow[calib_ind,]$Snowpack[mod_ind]>param[6],obs_cover)
        return(1-bi)
      }
    }else{
      return(c(0,0,1,5,0.1,5,NA,NA)) # no benchmark: default param & NA metric
    }
    set.seed(1234)
    outDEoptim <- DEoptim(hbvSnow,lower_par,upper_par,DEoptim.control(NP=NP,itermax=itermax))
    ## metric
    param=outDEoptim$optim$bestmem
    if(useDepth&useCover){
      dt_snow=runModel(df=dt,I_Snowpack0=10,I_Meltwater0=10,
                       P_TP=param[1],P_TM=param[2],P_SCF=param[3],P_CFX=param[4],P_CFR=0.05,P_CWH=param[5])
      kg=KGE(dt_snow[calib_ind,]$Snowpack[depth_ind],obs_depth)
      bi=biclass(dt_snow[calib_ind,]$Snowpack[mod_ind]>param[6],obs_cover)
    }else if(useDepth){
      kg=1-outDEoptim$optim$bestval
      bi=NA
    }else if(useCover){
      kg=NA
      bi=1-outDEoptim$optim$bestval
    }
    return(c(param,kg,bi))
  }
  ## load tile data
  prcp_date=seq(as.Date('2000-01-01'),as.Date('2017-12-31'),1)
  tair_date=seq(as.Date('2000-01-01'),as.Date('2019-12-31'),1)
  mod_date=fread(paste0(path,'mod_date.csv'))[,as.Date(mod_date)]
  depth_date=seq(as.Date('2000-01-01'),as.Date('2019-12-31'),1)
  r_prcp=readTIFF(paste0(path,'tile_prcp_cgdpav2/n',n,'e',e,'.tif'))
  r_prcp[r_prcp<0]=NA
  r_tair=readTIFF(paste0(path,'tile_tem_stationERA5/n',n,'e',e,'.tif'))
  r_tair[r_tair<(-999)]=NA
  r_cloud=readTIFF(paste0(path,'tile_cloud_mod10/n',n,'e',e,'.tif'))
  r_cloud[r_cloud<0]=NA
  r_cover=readTIFF(paste0(path,'tile_snowCover_mod10/n',n,'e',e,'.tif'))
  r_cover[r_cover<0]=NA
  r_depth=readTIFF(paste0(path,'tile_snowDepth_china/n',n,'e',e,'.tif'))*1.8
  r_depth[r_depth<0]=NA
  ind=expand.grid(1:dim(r_prcp)[1],1:dim(r_prcp)[2])
  ref=raster(paste0(path,'tile_prcp_cgdpav2/n',n,'e',e,'.tif'))
  ## parallel
  results=sapply(1:nrow(ind),function(i){
    tryCatch({
      return(gridCalib(i,20,20))
    }, error = function(e) {
      return(rep(-999,8)) # error: -999
    })
  })
  r_result=array(dim=c(nrow(ref),ncol(ref),8))
  for(i in 1:nrow(ind)){
    r_result[ind[i,1],ind[i,2],]=results[,i]
  }
  r_result=brick(r_result,xmn=xmin(ref),xmx=xmax(ref),ymn=ymin(ref),ymx=ymax(ref),crs=crs(ref))
  writeRaster(r_result,paste0(path,'tile_snowParam/n',n,'e',e,'.tif'),overwrite=TRUE)
  return(0)
}

### parallel calibration
path='./HBV_China/'
library(stringr)
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
foreach(x=setdiff(dir(paste0(path,'tile_prcp_cgdpav2')),dir(paste0(path,'tile_snowParam'))),.packages=c('stringr','data.table','tiff','hydroGOF','raster','Rcpp','RcppDE')) %dopar% {
  n=as.numeric(str_extract_all(x,'[0-9]+')[[1]][1])
  e=as.numeric(str_extract_all(x,'[0-9]+')[[1]][2])
  tileCalib(n,e,path)
}
stopCluster(cl)
registerDoSEQ()





