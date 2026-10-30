
# Load libraries ----------------------------------------------------------
library(pacman)
p_load(terra, tidyverse, glue, stringr, fs, RColorBrewer, ggspatial, sf, rnaturalearth, rnaturalearthdata)

rm(list = ls())
sf_use_s2(FALSE)


# Load data ---------------------------------------------------------------

# Parameters --------------------------------------------------------------
root <- '../tif/anom_yearly'
fles <- dir_ls(root)

## Template mask 
tmpl <- terra::rast('../tif/baseline/raw_mex/thetao_6stats_baseline_2010_depthsurf.tif')

## Mexico 
mex0 <- geodata::gadm(country = 'MEX', level = 1, path = './tmpr')
land <- geodata::world(resolution = 1, level = 0, path = './tmpr')
plot(mex0)

## Template adjust
tmpl <- terra::mask(tmpl, land, inverse = T, touches = F)

# Function to interpolate: Bilinear + IDW  --------------------------------
anom <- terra::rast(fles[1])
tmpl <- tmpl
mxds <- 50

intl.bilf <- function(anom, tmpl, mxds = 50){
  
  ## To start 
  cat('To start the blinear + IDW fill: ', names(anom)[1], ' (', terra::nlyr(anom), ' layers) \n', sep = '')
  
  ## Bilinear in the open ocean 
  bil <- terra::resample(anom, tmpl, method = 'bilinear')
  plot(bil[[1]])
  
  ## IDW (nearest, power) 
  out <- map(1:nlyr(anom), function(i){
    cat('Layer : ', i, '\t')
    pts <- terra::as.data.frame(anom[[i]], xy = T, na.rm = T)
    pts <- as.matrix(pts)
    idw <- terra::interpIDW(tmpl, pts, radius = mxds / 111, power = 2, maxPoints = 8)
    rst <- terra::cover(bil[[i]], idw) 
  }) |> 
    rast() |> 
    mask(tmpl)
  
  plot(out[[1]])
  
  out <- terra::crop(out, anom[[1]], snap = 'out')
  
  ## To add the names and the time
  names(out) <- names(anom)
  if(terra::has.time(anom)) terra::time(out, tstep = terra::timeInfo(anom)$step) <- terra::time(anom)
  
  ## Finish 
  cat('Done!\n')
  return(out)
  
}

# Function to interpolate anomalies with TPS (Spline) ---------------------
anom <- terra::rast(fles[1]); tmpl <- tmpl; basn <- NULL; knts <- 800; mxds <- 50

intp.tps  <- function(anom, tmpl, basn = NULL, knts = 800, mxds = 50){
  
  ## To start 
  cat('To start the TPS: ', names(anom)[1], ' (', terra::nlyr(anom), ' layers)\n', sep = '')
  cos <- cos(mean(as.vector(terra::ext(tmpl))[3:4]) * pi / 180)   # Degrees [Coseno]
  
  ## Target cells: ocean cells of the template
  trg <- terra::as.data.frame(tmpl, xy = T, cells = T, na.rm = T) %>% dplyr::select(cell, x, y)
  
  
  
}



