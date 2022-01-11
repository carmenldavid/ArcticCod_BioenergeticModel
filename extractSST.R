# ========================================== #
# SST  ============ #
# ========================================== #
# function to extract daily SST values between two given dates

extractDailySST <- function(file_sst, date_start, date_end){
  require(ncdf4)
  
  file=nc_open(file_sst) 
  sst <- ncvar_get(file, varid = "sst") 
  lat <- ncvar_get(file, varid = "latitude")
  lon <- ncvar_get(file, varid = "longitude")
  dts <- ncvar_get(file, varid = "time")
  # extract only date, cut off time
  dts_date = substr(as.POSIXlt(structure(dts,class=c('POSIXc','POSIXct'))), 1,10) 
  rm(file_sst)
  
  dts1 <- which(dts_date == date_start)
  dts2 <- which(dts_date == date_end)
  
  sst_by_date <- sst[,,dts1:dts2] 
  return(sst_by_date)
}



