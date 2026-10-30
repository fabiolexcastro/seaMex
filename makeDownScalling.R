
# Load libraries ----------------------------------------------------------
library(pacman)
p_load(terra, tidyverse, glue, stringr, fs, RColorBrewer, ggspatial, sf, rnaturalearth, rnaturalearthdata)

rm(list = ls())
sf_use_s2(FALSE)


# Load data ---------------------------------------------------------------

# Parameters --------------------------------------------------------------
root <- './tif/anom'
fles <- dir_ls(root)

## Template mask 
tmpl <- terra::rast('./tif/baseline/raw_mex/thetao_6stats_baseline_2010_depthsurf.tif')



