SWAToutput=function(outnetcdf_cf=TRUE,outsqlite=TRUE,outplot=TRUE,cleanswatoutput=FALSE){
# outnetcdf_cf=TRUE;outsqlite=TRUE;outplot=TRUE;cleanswatoutput=FALSE
# Requires being in a SWAT TxtInOut Dir and WXData dataframe
  output_hru=readSWAT("hru",".")
  output_sub=readSWAT("sub",".")
  output_rch=readSWAT("rch",".")
# SQLite output
  if(outsqlite){
  sqlitefile=paste0("./",args$swatscen,"MINTSWATtables.sqlite")
  con <- dbConnect(RSQLite::SQLite(), sqlitefile)
  dbWriteTable(con, "output_hru", output_hru,overwrite = TRUE)
  dbWriteTable(con, "output_rch", output_rch,overwrite = TRUE)
  dbWriteTable(con, "output_sub", output_sub,overwrite = TRUE)
  dbListTables(con)
  }
  if(outplot){
  # Graphical Summary
  output_plot=merge(output_rch,WXData,by.x="mdate",by.y="date")
  output_plot$Qpredmm=output_plot$FLOW_OUTcms/(basin_area*10^6)*3600*24*1000
  output_plot$Qmm=output_plot$Qm3ps/(basin_area*10^6)*3600*24/10
  
  maxRange <- 1.1*(max(output_plot$P,na.rm = T) + max(output_plot$Qpredmm,na.rm = T))
  
  p1<- ggplot() +
    # Use geom_tile to create the inverted hyetograph. geom_tile has a bug that displays a warning message for height and width, you can ignore it.
    geom_tile(data = output_plot, aes(x=mdate,y = -1*(P/2-maxRange), # y = the center point of each bar
                                      height = P,width = 1),
              fill = "black",
              color = "black") +
    # Plot your discharge data
    geom_line(data=output_plot,aes(x=mdate, y = Qpredmm, colour= "Qpred"), size=1) +
    scale_colour_manual("",breaks = c("Qmm", "Qpred"),values = c("red", "blue")) +
    # Create a second axis with sec_axis() and format the labels to display the original precipitation units.
    scale_y_continuous(name = "Discharge (mm/day)",
                       sec.axis = sec_axis(trans = ~-1*(.-maxRange),
                                           name = "Precipitation (mm/day)"))+
    scale_x_continuous(name = NULL,labels = NULL)+
    ggtitle(toupper(basinname))
  pdf(file = paste0(basinoutdir,"/","HydroSummary.pdf"),width = 6,height = 4)
  plot(p1)
  dev.off()
  }
  # Creating a NetCDF-CF of the output_rch dataframe
  #
  if(outnetcdf_cf){
    print(subs1_shp_ll)
    subs1_shp_ll_sf=st_as_sf(subs1_shp_ll)
    d = data.frame(HRU = as.character(1:3),HRU_NAME=paste0("HRU_",1:3))
    d$geom = subs1_shp_ll_sf$geom
    subs1_shp_ll_sf=st_as_sf(d)
    subs1_centroids <- subs1_shp_ll_sf %>%
      st_transform(5070) %>% # Albers Equal Area
      st_set_agr("constant") %>%
      st_centroid() %>%
      st_transform(4269) %>% #NAD83 Lat/Lon
      st_coordinates() %>%
      as.data.frame()
    subs1_poly <- st_sf(st_cast(subs1_shp_ll_sf, "MULTIPOLYGON"))
    
    nc_file <- "swat_rch.nc"
    unlink(nc_file)
    
    for(rch_var in colnames(dplyr::select(output_rch,-INFO,-mdate,-RCH))){
      print(rch_var)
      var_data=reshape(dplyr::select(output_rch,mdate,RCH,eval(rch_var)),
                       v.names=eval(rch_var),idvar="mdate",
                       timevar="RCH",direction="wide")
      
      colnames(var_data)=c("date",str_split(
        colnames(var_data[,2:length(colnames(var_data))]),
        "\\.",simplify=T)[,2])
      var_dates <- var_data$date
      var_data <- dplyr::select(var_data, -date)
      var_meta <- list(name = gsub("[a-z]","",rch_var,ignore.case = FALSE), 
                       long_name = rch_var)
      
      write_timeseries_dsg(nc_file = nc_file, 
                           instance_names = subs1_shp_ll_sf$HRU, 
                           lats = subs1_centroids$Y, 
                           lons = subs1_centroids$X, 
                           times = var_dates, 
                           data = var_data, 
                           data_unit = rep(gsub("[A-Z]|_","",rch_var,ignore.case = FALSE), (ncol(var_data) - 1)), 
                           data_prec = "float", 
                           data_metadata = var_meta, 
                           attributes = list(title = "SWAT Output RCH using ncdfgeom"), 
                           add_to_existing = TRUE)
      
    }
    
    ncvariables=list(gsub("[a-z]","",colnames(dplyr::select(
      output_rch,-INFO,-mdate,-RCH)),ignore.case = FALSE))
    write_geometry(nc_file = nc_file, 
                   geom_data = subs1_poly,
                   variables = ncvariables)
  }  
  if(cleanswatoutput){
    unlink("output.*")
    unlink("watout.dat")
  }
}
