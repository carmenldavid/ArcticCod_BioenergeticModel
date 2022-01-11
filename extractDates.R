# ============ #
# SST Dates    #
# ============ #
# function to extract a string of dates in the format of SST files

extractDates <- function(file_sst, date_start, date_end){
  require(ncdf4)
  
  file = nc_open(file_sst)
  dts <- ncvar_get(file, varid = "time")
  # extract only date, cut off time
  dts_date = substr(as.POSIXlt(structure(dts,class=c('POSIXc','POSIXct'))), 1,10) 
  rm(file_sst)
  
  dts1 <- which(dts_date == date_start)
  dts2 <- which(dts_date == date_end)
  
  dates <- dts_date[dts1:dts2] 
  return(dates)
}