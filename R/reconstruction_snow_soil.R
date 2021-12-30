
#-------------------------------------------------------------------------------#

#' @author Wencong Yang
#' @description reconstruction historical snow and soil moisture

#-------------------------------------------------------------------------------#

path_P='./sim_raster_year/'
path_PET='./year_PET_PT_mergeTEM/'
path_T='./year_tem_ele_era_nmax50/'
path_LAI='./year_China01/'
path_out='./cgdpav2_driven/'

require(tiff)
paramSnow=readTIFF('./snowParam.tif')[,,1:5]
paramSnow[paramSnow<(-999)]=NA
paramSoilWater=readTIFF('./soilWaterParam.tif')[,,1:4]
paramSoilWater[paramSoilWater<(-999)]=NA

#' @description historical simulation of snow and soil moisture variables
#' @param yr simulation year
#' @param path_P path for precipitation
#' @param path_PET path for potential evatranspiration
#' @param path_T path for air temperature
#' @param path_LAI path for leaf area index
#' @param path_out output path
#' @param paramSnow snow module parameter arrays
#' @param paramSoilWater soil water module parameter arrays
#' @param Snowpack_init initial snowpack array
#' @param Meltwater_init initial melt water storage array
#' @param Si_init initial interception water storage array
#' @param Si_init initial soil water storage array
hbv_year=function(yr,path_P,path_PET,path_T,path_LAI,path_out,paramSnow,paramSoilWater,
                  Snowpack_init,Meltwater_init,Si_init,Sup_init){
  require(data.table)
  require(raster)
  require(Rcpp)
  require(tiff)
  require(Matrix)
  require(abind)
  sourceCpp('cpp/hbv_soil.cpp')
  ### load forcing
  r_prcp=readTIFF(paste0(path_P,yr,'.tif'))
  r_prcp[r_prcp<0]=NA
  r_tair=readTIFF(paste0(path_T,yr,'.tif'))
  r_tair[r_tair<(-999)]=NA
  r_pet=readTIFF(paste0(path_PET,yr,'.tif'))
  r_pet[r_pet<(-999)]=NA
  r_pet[r_pet<0]=0 # negative pet
  r_lai=readTIFF(paste0(path_LAI,yr,'.tif'))
  r_lai[r_lai<(-999)]=0
  ### grid sim
  ind=expand.grid(1:dim(r_prcp)[1],1:dim(r_prcp)[2])
  dates=seq(as.Date(paste0(yr,'-01-01')),as.Date(paste0(yr,'-12-31')),1)
  dt_model=rbindlist(lapply(1:nrow(ind),function(irow){
    i=ind[irow,1]
    j=ind[irow,2]
    dt=data.table(prcp=r_prcp[i,j,],pet=r_pet[i,j,],tair=r_tair[i,j,],lai=r_lai[i,j,])
    if(nrow(na.omit(dt))!=0){
      dt_model=runModel(df=dt,I_Snowpack0=Snowpack_init[i,j],I_Meltwater0=Meltwater_init[i,j],
                        I_Si0=Si_init[i,j],I_Su0p=Sup_init[i,j],
                        P_TP=paramSnow[i,j,1],P_TM=paramSnow[i,j,2],P_SCF=paramSnow[i,j,3],
                        P_CFX=paramSnow[i,j,4],P_CFR=0.05,P_CWH=paramSnow[i,j,5],
                        P_mc=paramSoilWater[i,j,1],P_FC=paramSoilWater[i,j,2],
                        P_beta=paramSoilWater[i,j,3],P_LP=paramSoilWater[i,j,4])
      return(data.table(dt_model)[,`:=`(i=i,j=j,date=dates)])
    }else{
      return(data.table(Snowpack=NA,Meltwater=NA,Si=NA,Sup=NA,Rain=NA,Snow=NA,Melt=NA,Tosoil=NA,Ei=NA,Ea=NA,i=i,j=j,date=dates))
    }
  }))
  rm(r_prcp,r_pet,r_tair,r_lai)
  gc()
  ### output annual data
  crs_out=crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  # RAIN
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Rain))
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'RAIN/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  # SNOW
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Snow))
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'SNOW/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  # SNOWPACK
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Snowpack))
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'SNOWPACK/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  # MELT
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Melt))
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'MELT/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  # SM
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Sup))*paramSoilWater[,,2]/1000
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'SM/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  # ET
  r=abind(lapply(dates,function(d){
    dt=dt_model[date==d]
    as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt[,Ei+Ea]))
  }),along=3)
  r=brick(r,crs=crs_out,xmn=72,xmx=140,ymn=18,ymx=56)
  writeRaster(r,paste0(path_out,'ET/',yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
  ### update init
  dt=dt_model[date==dates[length(dates)]]
  State_init=list(as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Snowpack)),
                  as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Meltwater)),
                  as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Si)),
                  as.matrix(sparseMatrix(i=dt$i,j=dt$j,x=dt$Sup)))
  rm(dt_model)
  gc()
  return(State_init)
}

### initial states
State_init=list(array(10,dim=c(380,680)),array(10,dim=c(380,680)),array(0,dim=c(380,680)),array(0.5,dim=c(380,680)))

### simulation with a 2-year warm-up
for(yr in c(1981,1981,1981:2017)){
  State_init=hbv_year(yr,path_P,path_PET,path_T,path_LAI,path_out,paramSnow,paramSoilWater,
                      State_init[[1]],State_init[[2]],State_init[[3]],State_init[[4]])
  gc()
}
