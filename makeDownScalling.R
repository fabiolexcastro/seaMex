
# Load libraries ----------------------------------------------------------
library(pacman)
p_load(terra, tidyverse, glue, stringr, fs, RColorBrewer, ggspatial, sf, rnaturalearth, rnaturalearthdata)

rm(list = ls())
sf_use_s2(FALSE)


# Load data ---------------------------------------------------------------


