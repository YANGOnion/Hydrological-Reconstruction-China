
#-------------------------------------------------------------------------------#

#' @author Wencong Yang
#' @description calibration of HBV soil water module

#-------------------------------------------------------------------------------#

#' @description calibration of all grids in a 1 degree tile 
tileCalib=function(n,e,path){
  require(data.table)
  require(lubridate)
  require(raster)
  require(Rcpp)
  require(RcppDE)
  require(tiff)
  require(hydroGOF)
  sourceCpp('cpp/hbv_soil.cpp')
  ## calib on 1 grid
  gridCalib=function(i,NP=20,itermax=20){
    rowind=ind[i,1]
    colind=ind[i,2]
    dt_prcp=data.table(date=prcp_date,prcp=r_prcp[rowind,colind,])
    dt_tair=data.table(date=tair_date,tair=r_tair[rowind,colind,])
    dt_pet=data.table(date=pet_date,pet=r_pet[rowind,colind,])
    dt_lai=data.table(date=lai_date,lai=r_lai[rowind,colind,])
    dt_sw=data.table(date=sw_date,sw=r_sw[rowind,colind,])
    dt=Reduce(function(a,b){merge(a,b,by='date',all=T)},list(dt_prcp,dt_tair,dt_pet,dt_lai,dt_sw))
    dt=dt[!is.na(prcp)][!is.na(tair)][!is.na(pet)][!is.na(lai)] # filter non-land grid
    if(nrow(dt)==0) return(rep(NA,5)) # no forcing data: NA
    if(dt[,all(is.na(sw))]) return(c(0.3,300,2,0.7,NA)) # no benchmark: default param & NA metric
    ## obs
    dt=dt[year(date)>=2014]
    calib_ind=dt[,which(date>=as.Date('2015-04-01'))]
    sw=dt[calib_ind,sw]
    ## calibration
    snowParam=r_snowParam[rowind,colind,1:5]
    if(all(is.na(snowParam))) return(rep(NA,5)) # no snow params: NA
    lower_par=c(mc=0.1,FC=50,beta=0.1,LP=0.2)
    upper_par=c(mc=0.5,FC=800,beta=6,LP=1)
    hbvSoil=function(param){
      dt_soil=runModel(df=dt,I_Snowpack0=10,I_Meltwater0=10,I_Si0=0,I_Su0p=0.5,
                       P_TP=snowParam[1],P_TM=snowParam[2],P_SCF=snowParam[3],P_CFX=snowParam[4],P_CFR=0.05,P_CWH=snowParam[5],
                       P_mc=param[1],P_FC=param[2],P_beta=param[3],P_LP=param[4])
      kg=KGE(dt_soil[calib_ind,]$Sup*param[2]/1000,sw)
      if(is.na(kg)){
        kg=-10
      }
      return(1-kg)
    }
    set.seed(1234)
    outDEoptim <- DEoptim(hbvSoil,lower_par,upper_par,DEoptim.control(NP=NP,itermax=itermax))
    ## metric
    return(c(outDEoptim$optim$bestmem,1-outDEoptim$optim$bestval))
  }
  ## load tile data
  prcp_date=seq(as.Date('2000-01-01'),as.Date('2017-12-31'),1)
  tair_date=seq(as.Date('2000-01-01'),as.Date('2019-12-31'),1)
  pet_date=seq(as.Date('2000-01-01'),as.Date('2019-12-31'),1)
  lai_date=seq(as.Date('2000-01-01'),as.Date('2017-12-31'),1)
  sw_date=seq(as.Date('2015-04-01'),as.Date('2019-12-31'),1)
  r_prcp=readTIFF(paste0(path,'tile_prcp_cgdpav2/n',n,'e',e,'.tif'))
  r_prcp[r_prcp<0]=NA
  r_tair=readTIFF(paste0(path,'tile_tem_stationERA5/n',n,'e',e,'.tif'))
  r_tair[r_tair<(-999)]=NA
  r_pet=readTIFF(paste0(path,'tile_pet_stationERA5PT/n',n,'e',e,'.tif'))
  r_pet[r_pet<(-999)]=NA
  r_pet[r_pet<0]=0 # negative pet
  r_lai=readTIFF(paste0(path,'tile_lai_glass/n',n,'e',e,'.tif'))
  r_lai[r_lai<(-999)]=0
  r_sw=readTIFF(paste0(path,'tile_soilWater_smap/n',n,'e',e,'.tif'))
  r_sw[r_sw<(-999)]=NA
  r_snowParam=as.array(brick(paste0(path,'tile_snowParam/n',n,'e',e,'.tif')))
  ind=expand.grid(1:dim(r_prcp)[1],1:dim(r_prcp)[2])
  ref=raster(paste0(path,'tile_prcp_cgdpav2/n',n,'e',e,'.tif'))
  ## parallel
  results=sapply(1:nrow(ind),function(i){
    tryCatch({
      return(gridCalib(i,20,20))
    }, error = function(e) {
      return(rep(-999,5)) # error: -999
    })
  })
  r_result=array(dim=c(nrow(ref),ncol(ref),5))
  for(i in 1:nrow(ind)){
    r_result[ind[i,1],ind[i,2],]=results[,i]
  }
  r_result=brick(r_result,xmn=xmin(ref),xmx=xmax(ref),ymn=ymin(ref),ymx=ymax(ref),crs=crs(ref))
  writeRaster(r_result,paste0(path,'tile_soilWaterParam/n',n,'e',e,'.tif'),overwrite=TRUE)
  return(0)
}


### parallel calibration
library(stringr)
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
foreach(x=setdiff(dir(paste0(path,'tile_prcp_cgdpav2')),dir(paste0(path,'tile_soilWaterParam'))),.packages=c('stringr','data.table','lubridate','tiff','hydroGOF','raster','Rcpp','RcppDE')) %dopar% {
  n=as.numeric(str_extract_all(x,'[0-9]+')[[1]][1])
  e=as.numeric(str_extract_all(x,'[0-9]+')[[1]][2])
  tileCalib(n,e,path)
}
stopCluster(cl)
registerDoSEQ()

