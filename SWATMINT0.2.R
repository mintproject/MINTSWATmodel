# Add NSE and Observed for cases with observed flow
# 
library(parallel)
cl <- makeCluster(4)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(SWATmodel,RSQLite,argparse,stringi,stringr,rgdal,ggplot2,rgeos,rnoaa,moments,sf,readr,tools,
               diffobj,png,grid,gridExtra,ncdfgeom,purrr,raster)
source("./MINTSWATcalib.R")
source("./SWAToutput.R")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/get_grdc_gage.R?root=ecohydrology")
source("https://raw.githubusercontent.com/cran/sampSurf/master/R/spCircle.R")
setwd("~")
basedir=getwd()
outbasedir=paste0(basedir,"/MINTSWATmodel_output")
inbasedir=paste0(basedir,"/MINTSWATmodel_input")
print(detectCores())
dir.create(outbasedir)
dir.create(inbasedir)
setwd(inbasedir)
Sys.setenv(R_USER_CACHE_DIR=inbasedir)
# If a parameter change scenario, we use --swatscen
parser <- ArgumentParser()
parser$add_argument("-p","--swatparam", action="append", metavar="param:val[:regex_file]",
                    help = "Add in SWAT parameters that need to be modified")
parser$add_argument("-s","--swatscen", metavar="calib01",
                    help = "Scenario folder name, 'calib' for calibration, 'scen' for scenario")
parser$add_argument("-d","--swatiniturl", metavar="url to ArcSWAT init or GRDC format dataset",
                    help = "Scenario folder name")
# Examples:
# geojson example 
exampleargs=c("-d https://data.mint.isi.edu/files/files/geojson/guder.json")
# GRDC Calibration example 
exampleargs=c("-s calib01","-p deiter:10","-p rch:3","-d https://portal.grdc.bafg.de/grdcdownload/external/7ce24ffd-3c99-407f-84e8-bd4a99417c06/2022-07-08_00-05.zip")
# ArcSWAT example 
#exampleargs=c("-d https://raw.githubusercontent.com/vtdrfuka/MINTSWATmodel/main/tb_s2.zip")
#
args <- parser$parse_args()
if(is.null(args$swatiniturl)){
  args <- parser$parse_args(c(exampleargs))
}
print(paste0("This run's args: ",args))
dlfilename=basename(args$swatiniturl)
dlurl=trimws(args$swatiniturl)

paramloc=grep("deiter",args$swatparam)
if(length(paramloc)>0){
  deiter=as.numeric(strsplit(args$swatparam[paramloc],split = ":")[[1]][2])
}else{
  deiter=200
}
paramloc=grep("rch",args$swatparam)
if(length(paramloc)>0){
  rch=as.numeric(strsplit(args$swatparam[paramloc],split = ":")[[1]][2])
}else{
  rch=3
}

# *** download
dlfiletype=file_ext(dlfilename)
if(dlfiletype=="json"){
  print("geojson single run")
  download.file(dlurl,paste0("data.",dlfiletype))
  swatrun="basic"
} else {
  print("different")
  dlfiletype="zip"
  download.file(dlurl,paste0("data.",dlfiletype))
  if(grepl("Q_Day",unzip("data.zip", list=T)[1])){
    swatrun="GRDC"
  }    
}

if(swatrun=="GRDC"){
  print("GRDC Format Uninitialized")
  dir.create("GRDCstns")
  setwd("GRDCstns")
  currentdir=getwd()
  unzip("../data.zip")
  tryCatch({
    stationbasins_shp=readOGR("stationbasins.geojson")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  for(grdcfilename in list.files(pattern = "_Q_Day")){
    print(grdcfilename)    
    setwd(currentdir)
    flowgage=get_grdc_gage(lfilename =  grdcfilename)
    if(is.null(flowgage)){print("Not enough Gage Info");next()}
    basinid=strsplit(grdcfilename,"_")[[1]][1]
    basinname=basinid
    if(is.character(flowgage)){next()}
    GRDC_mindate=min(flowgage$flowdata$mdate)
    GRDC_maxdate=max(flowgage$flowdata$mdate)
    # Depends on: rnoaa, lubridate::month,ggplot2
    declat=flowgage$declat
    declon=flowgage$declon
    proj4_ll = "+proj=longlat"
    proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")
    crs_ll=CRS(proj4_ll)
    crs_utm=CRS(proj4_utm)
    basin_area=flowgage$area
    # Building 3 basin Feature for NetCDF based on basin shape if available
    # or a virtual circular basins on area and outlet. 
    if(exists("stationbasins_shp")&&any(stationbasins_shp@data$grdc_no==basinid)){
      subs1_shp=subset(stationbasins_shp,grdc_no==basinid)
      proj4_utm = paste0("+proj=utm +zone=", trunc((180+gCentroid(subs1_shp)$x)/6+1), " +datum=WGS84 +units=m +no_defs")
      proj4_ll = "+proj=longlat"
      crs_ll=CRS(proj4_ll)
      crs_utm=CRS(proj4_utm)
      subs1_shp_utm=spTransform(subs1_shp,crs_utm)
      initsizeguess=-sqrt(gArea(subs1_shp_utm))/6
      f <- function (x,a) {(gArea(subs1_shp_utm)*a-
                              gArea(gBuffer(subs1_shp_utm,width=x)))^2}
      hru3scale <- optimize(f, c(initsizeguess, 0), tol = 0.0001,a=1/3)$minimum
      hru2scale <- optimize(f, c(initsizeguess, 0), tol = 0.0001,a=2/3)$minimum
      hru3_utm=gBuffer(subs1_shp_utm,width=hru3scale)
      hru2_utm=gBuffer(subs1_shp_utm,width=hru2scale)
      gArea(hru2_utm)/gArea(subs1_shp_utm)
      gArea(hru3_utm)/gArea(subs1_shp_utm)
      hru1_utm=gDifference(subs1_shp_utm,hru2_utm)
      hru2_utm=gDifference(hru2_utm,hru3_utm)
      
      combined_hrus=list(c(hru1_utm,hru2_utm,hru3_utm))
      list(combined_hrus, makeUniqueIDs = T) %>% 
        flatten() %>% 
        do.call(rbind, .)
      subs1_shp_utm <- do.call(bind, combined_hrus)
      subs1_shp_ll=spTransform(subs1_shp_utm,crs_ll)
    } else{
      latlon <- cbind(declon,declat)
      gagepoint_ll <- SpatialPoints(latlon)
      proj4string(gagepoint_ll)=proj4_ll
      gagepoint_utm=spTransform(gagepoint_ll,crs_utm)
      hru1_utm=spCircle(sqrt(basin_area*1000^2/pi), spUnits = crs_utm,
                        centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                        nptsPerimeter = 30,spID = 1)$spCircle
      hru2_utm=spCircle(sqrt(basin_area*2/3*1000^2/pi), spUnits = crs_utm,
                        centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                        nptsPerimeter = 30,spID = 1 )$spCircle
      hru3_utm=spCircle(sqrt(basin_area/3*1000^2/pi), spUnits = crs_utm,
                        centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                        nptsPerimeter = 30,spID = 1 )$spCircle
      hru1_utm=gDifference(hru1_utm,hru2_utm)
      hru2_utm=gDifference(hru2_utm,hru3_utm)
      
      combined_hrus=list(c(hru1_utm,hru2_utm,hru3_utm))
      list(combined_hrus, makeUniqueIDs = T) %>% 
        flatten() %>% 
        do.call(rbind, .)
      subs1_shp_utm <- do.call(bind, combined_hrus)
      proj4string(subs1_shp_utm)=proj4_utm
      subs1_shp_ll=spTransform(subs1_shp_utm,crs_ll)
      
    }
    
    if(exists("stationbasins_shp")&&length(try(which(stationbasins_shp$grdc_no==as.numeric(flowgage$id))))>0){
      basinloc=which(stationbasins_shp$grdc_no==as.numeric(flowgage$id))
      basin=stationbasins_shp[basinloc,]
      basinutm=spTransform(basin,CRS(proj4_utm))
      wxlat=gCentroid(basin)$y
      wxlon=gCentroid(basin)$x
    } else {
      wxlat=declat
      wxlon=declon
    }
    stradius=20;minstns=30
    station_data=ghcnd_stations()
    while(stradius<2000){
      print(paste0("Looking for WX Stations within: ",stradius,"km"))
      junk=meteo_distance(
        station_data=station_data,
        lat=wxlat, long=wxlon,
        units = "deg",
        radius = stradius,
        limit = NULL
      )
      if(length(unique(junk$id))>minstns){break()}
      stradius=stradius*1.2
    }
    basinoutdir=paste0(outbasedir,"/",basinid);dir.create(basinoutdir)
    dir.create(basinoutdir,recursive = T)
    setwd(basinoutdir)
    WXData=FillMissWX(wxlat,wxlon,date_min = "1979-01-01",date_max = "2022-01-01", StnRadius = stradius,method = "IDW",alfa = 2)
    
    GRDC_mindate=min(WXData$date)
    GRDC_maxdate=max(WXData$date)
    AllDays=data.frame(date=seq(GRDC_mindate, by = "day", length.out = GRDC_maxdate-GRDC_mindate))
    WXData=merge(AllDays,WXData,all=T)
    
    WXData$PRECIP=WXData$P
    WXData$PRECIP[is.na(WXData$PRECIP)]=-99
    WXData$TMX=WXData$MaxTemp
    WXData$TMX[is.na(WXData$TMX)]=-99
    WXData$TMN=WXData$MinTemp
    WXData$TMN[is.na(WXData$TMN)]=-99
    WXData$DATE=WXData$date
    build_swat_basic(dirname=basinoutdir, iyr=min(year(WXData$DATE),na.rm=T),    ###***basin name!
                     nbyr=(max(year(WXData$DATE),na.rm=T)-min(year(WXData$DATE),na.rm=T)), 
                     wsarea=basin_area, elev=mean(WXData$prcpElevation,na.rm=T), 
                     declat=declat, declon=declon, hist_wx=WXData)
    build_wgn_file(metdata_df=WXData,declat=declat,declon=declon)
    
    if(!is.null(args$swatscen) && 
       substr(trimws(args$swatscen),1,5)=="calib"){
      MINTSWATcalib()
    }
    
    runSWAT2012()
    SWAToutput()
  }
}

if(dlfiletype=="json"){
  basinname=strsplit(basename(args$swatiniturl),split = "\\.")[[1]][1]
  basinoutdir=paste0(outbasedir,"/",basinname);dir.create(basinoutdir)
  basin=readOGR("data.json")
  declat=gCentroid(basin)$y
  declon=gCentroid(basin)$x
  proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")
  proj4_ll = "+proj=longlat"
  crs_ll=CRS(proj4_ll)
  crs_utm=CRS(proj4_utm)
  basinutm=spTransform(basin,CRS(proj4_utm))
  basin_area=gArea(basinutm)/10^6
  
  # Replace with conversion of geojson
  latlon <- cbind(declon,declat)
  gagepoint_ll <- SpatialPoints(latlon)
  proj4string(gagepoint_ll)=proj4_ll
  gagepoint_utm=spTransform(gagepoint_ll,crs_utm)
  hru1_utm=spCircle(sqrt(basin_area*1000^2/pi), spUnits = crs_utm,
                    centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                    nptsPerimeter = 30,spID = 1)$spCircle
  hru2_utm=spCircle(sqrt(basin_area*2/3*1000^2/pi), spUnits = crs_utm,
                    centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                    nptsPerimeter = 30,spID = 1 )$spCircle
  hru3_utm=spCircle(sqrt(basin_area/3*1000^2/pi), spUnits = crs_utm,
                    centerPoint = c(x = gagepoint_utm@coords[1], y = gagepoint_utm@coords[2]),
                    nptsPerimeter = 30,spID = 1 )$spCircle
  hru1_utm=gDifference(hru1_utm,hru2_utm)
  hru2_utm=gDifference(hru2_utm,hru3_utm)
  
  combined_hrus=list(c(hru1_utm,hru2_utm,hru3_utm))
  list(combined_hrus, makeUniqueIDs = T) %>% 
    flatten() %>% 
    do.call(rbind, .)
  subs1_shp_utm <- do.call(bind, combined_hrus)
  proj4string(subs1_shp_utm)=proj4_utm
  subs1_shp_ll=spTransform(subs1_shp_utm,crs_ll)
  # End replace with GeoJSON conversion
  stradius=20;minstns=30
  station_data=ghcnd_stations()
  while(stradius<2000){
    print(stradius)
    junk=meteo_distance(
      station_data=station_data,
      lat=gCentroid(basin)$y, long=gCentroid(basin)$x,
      units = "deg",
      radius = stradius,
      limit = NULL
    )
    if(length(unique(junk$id))>minstns){break()}
    stradius=stradius*1.2
  }
  setwd(basinoutdir)
  WXData=FillMissWX(gCentroid(basin)$y,gCentroid(basin)$x,date_min = "1979-01-01",date_max = "2022-01-01", StnRadius = stradius,method = "IDW",alfa = 2)
  GRDC_mindate=min(WXData$date)
  GRDC_maxdate=max(WXData$date)
  AllDays=data.frame(date=seq(GRDC_mindate, by = "day", length.out = GRDC_maxdate-GRDC_mindate))
  WXData=merge(AllDays,WXData,all=T)
  
  WXData$PRECIP=WXData$P
  WXData$PRECIP[is.na(WXData$PRECIP)]=-99
  WXData$TMX=WXData$MaxTemp
  WXData$TMX[is.na(WXData$TMX)]=-99
  WXData$TMN=WXData$MinTemp
  WXData$TMN[is.na(WXData$TMN)]=-99
  WXData$DATE=WXData$date
  build_swat_basic(dirname=basinoutdir, iyr=min(year(WXData$DATE),na.rm=T), 
                   nbyr=(max(year(WXData$DATE),na.rm=T)-min(year(WXData$DATE),na.rm=T)), 
                   wsarea=basin_area, elev=mean(WXData$prcpElevation,na.rm=T), 
                   declat=declat, declon=declon, hist_wx=WXData)
  build_wgn_file(metdata_df=WXData,declat=declat,declon=declon)
  runSWAT2012()
  SWAToutput()
}

setwd(outbasedir)
unlink(list.files(pattern = "output.*",recursive = TRUE))
unlink(list.files(pattern = "*.out",recursive = TRUE))
unlink(list.files(pattern = "pcp1.pcp",recursive = TRUE))
unlink(list.files(pattern = "tmp1.tmp",recursive = TRUE))
quit()
