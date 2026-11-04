
# Load libraries ----------------------------------------------------------
library(pacman)
p_load(terra, tidyverse, glue, stringr, mgcv, fs, RColorBrewer, ggspatial, sf, rnaturalearth, rnaturalearthdata)

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
  dout <- './tif/anom_005/tps'; dir_create(dout)
  cat('To start the TPS: ', names(anom)[1], ' (', terra::nlyr(anom), ' layers)\n', sep = '')
  cos <- cos(mean(as.vector(terra::ext(tmpl))[3:4]) * pi / 180)   # Degrees [Coseno]
  
  ## Target cells: ocean cells of the template
  trg <- terra::as.data.frame(tmpl, xy = T, cells = T, na.rm = T) %>% dplyr::select(cell, x, y)
  
  ## Target cells / ocean cell
  cvr <- terra::resample(anom[[1]], tmpl, method = 'near')
  dst <- terra::distance(cvr)
  trg <- trg %>% mutate(dist = terra::extract(dst, cell)[, 1] / 1000) %>% filter(dist <= mxds)
  
  ## Simple conditional
  if(is.null(basn)){
    trg$bsn <- 0
    bsc     <- NULL
  }
  
  ## Layer by layer
  out <- map(1:nlyr(anom), function(i){
    
    ## Output file of the layers: it it already exists, akip it
    fle <- glue('{dout}/{names(anom)[i]}_tps005.tif')
    if(file.exists(fle)){
      cat('Layer ', i, ' of ', terra::nlyr(anom), ': already written -> skip\n', sep = '')
      return(terra::rast(fle))
    }
    
    cat('Layer ', i, ' of ', terra::nlyr(anom), '\n', sep = '')
    lyr <- anom[[i]]
    pnt <- terra::as.data.frame(lyr, xy = T, cells = T, na.rm = T) %>% setNames(c('cell', 'x', 'y', 'vlue'))
    nrow(pnt)
    
    ## A simple plot
    plot(anom[[1]] * 0 + 1)
    plot(mex0, add = T, border = 'grey80')
    points(pnt$x, pnt$y, pch = 16, col = 'red', cex = 0.1)
    
    ##
    pnt$bsn <- if(is.null(bsc)) 0 else terra::extract(bsc, pnt$cell)[, 1] %>% replace_na(0)
    
    ## TPS
    prd <- map_dfr(sort(unique(trg$bsn)), function(b){
      
      cat('Processing....!\t')
      pnb <- pnt |> filter(bsn == b) |> mutate(xs = x * cos)
      tgb <- trg |> filter(bsn == b) |> mutate(xs = x * cos)
      
      ## Conditional by points
      if(nrow(pnb) < 30){
        nrow(pnb)
        cat('Basin ', b, ': only ', nrow(pnb), ' coarse cells -> NA\n', sep = '')
        return(tgb %>% mutate(vlue = NA_real_))
      }
      
      kkk <- min(knts, floor(nrow(pnb) / 3))
      set.seed(1)
      
      ## Generalized additive models [GAM]
      mdl <- gam(vlue ~ s(xs, y, bs = 'tp', k = kkk), data = pnb, method = 'REML')
      edf <- sum(mdl$edf) - 1
      tgb$vlue <- predict(mdl, newdata = tgb) %>% as.numeric()
      return(tgb)
      
    })
    
    ##
    rst <- terra::rast(tmpl)
    rst[prd$cell] <- prd$vlue
    
    ## To write the raster
    terra::writeRaster(x = rst, filename = fle, overwrite = T)
    return(rst)
    
  })
  
  ## Same names and time as the input
  out <- terra::rast(out)
  names(out) <- names(anom)
  if(terra::has.time(anom)) terra::time(out, tstep = terra::timeInfo(anom)$step) <- terra::time(anom)
  
  cat('Done!\n')
  return(out)
  
}  

## Function to test the methods (perfect-model test) --------------------------
test.intp <- function(trth, fact = 3, basn = NULL, knts = 800){
  
  cat('To start!\n')
  # trth <- terra::rast(fles[1])[[1]]; fact <- 3; basn <- basn; knts <- 800
  frc <- terra::aggregate(!is.na(trth), fact = fact, fun = 'mean')
  crs <- terra::aggregate(trth, fact = fact, fun = 'mean', na.rm = T)
  crs <- terra::ifel(frc >= 0.5, crs, NA)
  
  ## Back to the fine grid with each method
  pts <- terra::as.data.frame(crs, xy = T, na.rm = T) %>% as.matrix()
  rsl <- list(
    near     = terra::resample(crs, trth, method = 'near'),
    bilinear = terra::resample(crs, trth, method = 'bilinear'),
    idw      = terra::interpIDW(trth, pts, radius = 3 * terra::res(crs)[1], power = 2, maxPoints = 8),
    tps      = intp.tps(anom = crs, tmpl = trth, basn = NULL, knts = knts),
    tps.basn = if(!is.null(basn)) intp.tps(anom = crs, tmpl = trth, basn = basn, knts = knts)
  ) %>% compact()
  
  
}


###




#