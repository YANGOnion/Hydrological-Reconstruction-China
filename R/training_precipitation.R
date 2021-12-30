
#-------------------------------------------------------------------------------#

#' @author Wencong Yang
#' @description tile-by-tile training for precipitation merging

#-------------------------------------------------------------------------------#

library(stringr)
path_training='./training/'
path_valid='./valid/'
path_historical='./historical/'
path_historical_sim='./historical_sim/'
grids=str_extract(dir(path_training),'[0-9_]+(?=\\.)')


#' @description get lon & lat for each grid in a 1 degree tile
#' @param n lower lat of the tile
#' @param e lower lon of the tile
#' @param gridlen size of the tile in degree
getLonLat=function(n,e,gridlen){
  require(raster)
  r=raster(ext=extent(e,e+gridlen,n,n+gridlen),res=0.1,
           crs=crs('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
  tile_lon=as.matrix(init(r, 'x'))
  tile_lat=as.matrix(init(r, 'y'))
  ind=expand.grid(1:dim(r)[1],1:dim(r)[2])
  dt=data.table(t(sapply(1:nrow(ind),function(i){
    rowind=ind[i,1]
    colind=ind[i,2]
    lon=tile_lon[rowind,colind]
    lat=tile_lat[rowind,colind]
    return(c(i,lon,lat))
  })))
  names(dt)=c('id','lon','lat')
  return(dt)
}


#' @description training a precipitation prediction model in a tile
#' @param grid_name grid name, e.g., "35_110" for lat=35 and lon=110
#' @param path_training path for inputing training data files
#' @param path_valid path for outputing validation data files
#' @param path_historical path for inputing historical data files
#' @param path_historical path for outputing historical prediction files
#' @param winter whether to use cold season benchmarks
predGrid=function(grid_name,path_training,path_valid,path_historical,path_historical_sim,winter=TRUE){
  require(caret)
  require(caretEnsemble)
  require(data.table)
  require(stringr)
  dt_grid=fread(paste0(path_training,grid_name,'.csv'))
  n=as.numeric(str_split(grid_name,'_')[[1]][1])
  e=as.numeric(str_split(grid_name,'_')[[1]][2])
  dt_lonlat=getLonLat(n,e,1)
  dt_grid=merge(dt_grid,dt_lonlat,by='id')
  features=c('prcp_rs','lai','lai_d','ele','ele_d','prcp_aux','lon','lat')
  if(!winter){
    features=setdiff(features,'lai')
    dt_grid=dt_grid[month(date)%in%(5:9)]
  }
  ### classification
  X=as.matrix(dt_grid[,features,with=F])
  y=dt_grid[,as.character(isRain)]
  set.seed(1)
  model=train(x=X,y=y,method='ranger',
              trControl=trainControl(method='cv',number=5,savePredictions='final',
                                     returnData=FALSE,trim=TRUE),
              tuneGrid=expand.grid(mtry=c(5),splitrule='gini',min.node.size=5),num.trees=100,num.threads=1)
  fwrite(cbind(dt_grid[,.(date,id)],data.table(model$pred)[order(rowIndex),.(pred,obs,Resample)]),
         paste0(path_valid,grid_name,'_cls.csv'))
  ### regression
  X=as.matrix(dt_grid[isRain==1,features,with=F])
  y=dt_grid[isRain==1,prcp]
  set.seed(1)
  trainIndex=createFolds(y,5,returnTrain=TRUE)
  model_list <- caretList(
    x=X,y=y,
    trControl=trainControl(method='cv',number=5,savePredictions='final',index=trainIndex,
                           returnData=FALSE,trim=TRUE,predictionBounds=c(0,NA)),
    tuneList=list(
      ranger=caretModelSpec(method="ranger",num.trees=100,num.threads=1,
                            tuneGrid=data.frame(.mtry=5,.splitrule='variance',.min.node.size=3)),
      lm=caretModelSpec(method="lm")
    )
  )
  set.seed(1)
  model_stack <- caretStack(
    model_list,
    method="brnn",
    trControl=trainControl(method='cv',number=5,savePredictions='final',
                           returnData=FALSE,trim=TRUE,predictionBounds=c(0,NA)),
    tuneGrid=expand.grid(neurons=2)
  )
  fwrite(cbind(dt_grid[isRain==1,.(date,id)],data.table(model_stack$ens_model$pred)[order(rowIndex),.(pred,obs,Resample)]),
         paste0(path_valid,grid_name,'_reg.csv'))
  ### sim
  dt_his=fread(paste0(path_historical,grid_name,'.csv'))
  dt_his=merge(dt_his,dt_lonlat,by='id')
  X=as.matrix(dt_his[,features,with=F])
  pred_cls=predict(model,X)
  pred_reg=predict(model_stack,X)
  fwrite(dt_his[,.(date,id,isRain=pred_cls,prcp=pred_reg)],
         paste0(path_historical_sim,grid_name,'.csv'))
  return(0)
}


### parallel training
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
foreach(grid_name=setdiff(grids,str_extract(dir(path_historical_sim),'.+(?=.csv)')),.packages=c('data.table','caret','caretEnsemble','raster','stringr'),.errorhandling='remove') %dopar% {
  n=as.numeric(str_split(grid_name,'_')[[1]][1])
  e=as.numeric(str_split(grid_name,'_')[[1]][2])
  if((n>=40)|((e<=100)&(n>=27))){
    winter=FALSE
  }else{
    winter=TRUE
  }
  predGrid(grid_name,path_training,path_valid,path_historical,path_historical_sim,winter)
}
stopCluster(cl)
registerDoSEQ()



