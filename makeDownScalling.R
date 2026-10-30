
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
  })
  
  
  
}




