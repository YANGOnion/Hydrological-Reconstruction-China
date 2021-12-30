
#-------------------------------------------------------------------------------#

#' @author Wencong Yang
#' @description reconstruction historical precipitation

#-------------------------------------------------------------------------------#

#' @description get lon & lat for each grid in a 1 degree tile
#' @param n lower lat of the tile
#' @param e lower lon of the tile
#' @param gridlen size of the tile in degree
getLonLat=function(n,e,gridlen){
  require(raster)
  r=raster(ext=extent(e,e+gridlen,n,n+gridlen),res=0.1,crs=crs('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
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

#------------------------ historical sim csv to tif ----------------------------# 

library(stringr)
path_historical_sim='./historical_sim/'
path_sim_raster='./sim_raster/'
path_sim_raster_year='./sim_raster_year/'
grids=str_extract(dir(path_historical_sim),'[0-9_]+(?=\\.)')

### historical predictions per tile per year
library(raster)
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
foreach(grid_name=grids,.packages=c('data.table','raster','stringr'),.errorhandling='remove') %dopar% {
  n=as.numeric(str_split(grid_name,'_')[[1]][1])
  e=as.numeric(str_split(grid_name,'_')[[1]][2])
  dt_lonlat=getLonLat(n,e,1)
  dt_sim=fread(paste0(path_historical_sim,n,'_',e,'.csv'))
  dt_sim[isRain==0,prcp:=0]
  dt_sim[prcp<0,prcp:=0]
  lapply(1981:2017,function(yr){
    dates=seq(as.Date(paste0(yr,'-01-01')),as.Date(paste0(yr,'-12-31')),1)
    dt=dt_sim[year(date)==yr]
    dt=merge(merge(data.table(expand.grid(id=1:100,date=dates)),
                   dt_lonlat,by='id'),dt,by=c('id','date'),all.x=T)
    r=brick(array(data = dt[order(lat,decreasing=T)][order(date,lon)]$prcp,dim=c(10, 10, length(dates))),
            xmn=e,xmx=e+1,ymn=n,ymx=n+1,crs=crs('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'))
    writeRaster(r,paste0(path_sim_raster,yr,'_',n,'_',e,'.tif'),format='GTiff',overwrite=TRUE)
  })
}
stopCluster(cl)
registerDoSEQ()

### merge all tiles per year
library(doParallel)
cl <- makeCluster(13)
registerDoParallel(cl)
foreach(yr=1981:2017,.packages=c('data.table','raster')) %dopar% {
  r=lapply(grids,function(grid_name){
    brick(paste0(path_sim_raster,yr,'_',grid_name,'.tif'))
  })
  r=do.call(merge,r)
  writeRaster(extend(r,extent(72,140,18,56)),paste0(path_sim_raster_year,yr,'.tif'),format='GTiff',overwrite=TRUE)
  rm(r)
  gc()
}
stopCluster(cl)
registerDoSEQ()

#---------------------------- fill no CGDPA days -------------------------------#

library(data.table)
library(raster)
library(lubridate)
library(tiff)

dt_miss=fread('./CGDPA/miss.txt') # missing days of CGDPA
dt_miss[,miss:=as.Date(as.character(miss),format='%Y%m%d')-1]
r_ref=raster(paste0(path_sim_raster_year,'2015.tif'))
r_ref[!is.na(r_ref)]=1

### fill with MSWEP: <2017-11-01
lapply(2015:2017, function(yr){
  r=readTIFF(paste0(path_sim_raster_year,yr,'.tif'))
  r[r<0]=NA
  r_aux=brick(paste0('./MSWEPv2_China/year01/',yr,'.tif'))
  for(x in dt_miss[year(miss)==yr][miss<=as.Date('2017-10-31'),yday(miss)]){
    r[,,x]=as.array(r_aux[[x]]*r_ref)
  } 
  r=brick(r,xmn=72,xmx=140,ymn=18,ymx=56,crs=crs(r_ref))
  names(r)=paste0('X',yr,'.',1:yday(paste0(yr,'-12-31')))
  writeRaster(r,paste0(path_sim_raster_year,yr,'.tif'),format='GTiff',overwrite=TRUE)
})

### fill with CMPA:>=2017-11-01
yr=2017
r=readTIFF(paste0(path_sim_raster_year,yr,'.tif'))
r[r<0]=NA
r_cmpa=brick(paste0('./CMPAv2/year01/',yr,'.tif'))
for(x in yday(seq(as.Date('2017-11-01'),as.Date('2017-12-31'),1))){
  r_prcp=r_cmpa[[x]]
  r_prcp[is.na(r_prcp)]=0
  r[,,x]=as.array(r_prcp*r_ref)
} 
r=brick(r,xmn=72,xmx=140,ymn=18,ymx=56,crs=crs(r_ref))
names(r)=paste0('X',yr,'.',1:yday(paste0(yr,'-12-31')))
writeRaster(r,paste0(path_sim_raster_year,yr,'.tif'),format='GTiff',overwrite=TRUE)



