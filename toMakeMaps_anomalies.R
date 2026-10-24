# Load libraries ----------------------------------------------------------
library(pacman)
p_load(terra, tidyverse, glue, stringr, fs, RColorBrewer, ggspatial, sf, rnaturalearth, rnaturalearthdata)

rm(list = ls())
sf_use_s2(FALSE)

# Parameters --------------------------------------------------------------
root <- './tif/anom_yearly'
dout <- './png/anom'
refp <- 2006:2020   ## The same reference period used in calc.anom
qntl <- 0.98        ## Colour limits: quantile of |anomaly| over all the months and years; beyond it the colour saturates
mlbl <- c('Enero', 'Febrero', 'Marzo', 'Abril', 'Mayo', 'Junio', 'Julio', 'Agosto', 'Septiembre', 'Octubre', 'Noviembre', 'Diciembre')

## Palette by variable: the last colour is the positive anomaly
## Red = stress direction (warming, deoxygenation, acidification)
pals <- tribble(
  ~vari,    ~labl,                 ~pall,  ~revs, ~unts,
  'thetao', 'temperatura del mar', 'RdBu', TRUE,  '°C',
  'so',     'salinidad',           'PuOr', FALSE, 'PSU',
  'chl',    'clorofila-a',         'BrBG', FALSE, 'mg/m³',
  'o2',     'oxígeno disuelto',    'RdBu', FALSE, 'mmol/m³',
  'ph',     'pH',                  'RdBu', FALSE, 'pH'
)

# Load data ---------------------------------------------------------------
fles <- dir_ls(root, regexp = '\\.tif$') |> as.character()

## Only the monthly outputs of this reference period
fles <- tibble(fles = fles) |>
  mutate(
    name = basename(fles) |> path_ext_remove(),
    vari = str_extract(name, '^[^_]+'),
    tstp = str_extract(name, '(month|year)(?=_ref)'),
    rela = str_detect(name, '_anompct_')
  ) |>
  filter(str_detect(name, glue('_ref{min(refp)}-{max(refp)}$')))

## Variables with only yearly anomalies cannot be mapped by month
noms <- setdiff(fles$vari, filter(fles, tstp == 'month')$vari)
if(length(noms) > 0) warning(glue('Sin anomalías mensuales (corre calc.anom con tstp = "month"): {paste(noms, collapse = ", ")}'))
fles <- filter(fles, tstp == 'month') |> distinct(vari, .keep_all = TRUE)

stopifnot('No hay anomalías mensuales con ese período de referencia en root' = nrow(fles) > 0)
print(select(fles, name, vari, tstp, rela))

## Land (Natural Earth 1:50m)
land <- ne_countries(scale = 'medium', returnclass = 'sf') |> st_geometry()

# To get the year and month of each layer ---------------------------------
get.ym <- function(rst){
  
  ## From the layer names written by calc.anom (thetao_2000_01_anom)
  nms  <- names(rst)
  year <- str_extract(nms, '(19|20)\\d{2}') |> as.numeric()
  mnth <- str_match(nms, '(19|20)\\d{2}_(\\d{2})(?!\\d)')[, 3] |> as.numeric()
  
  ## If the names don't have them, from the time
  if((anyNA(year) | anyNA(mnth)) & has.time(rst)){
    tinf <- time(rst)
    year <- as.numeric(format(tinf, '%Y'))
    mnth <- as.numeric(format(tinf, '%m'))
  }
  
  dtes <- tibble(year = year, mnth = mnth)
  stopifnot(
    'El raster debe ser mensual (año y mes en cada capa)' = !anyNA(dtes$year) & !anyNA(dtes$mnth),
    'Fechas repetidas (año-mes) en las capas' = anyDuplicated(dtes) == 0
  )
  return(dtes)
  
}

# To get the colour limits ------------------------------------------------
get.lims <- function(rst, qntl = 0.98){
  
  ## Symmetric around 0 and from all the layers, so every month and year share the scale
  vals <- values(rst, mat = FALSE)
  lims <- quantile(abs(vals), qntl, na.rm = TRUE) |> unname() |> signif(2)
  return(c(-lims, lims))
  
}

# To make the map of one year (12 months) ---------------------------------
make.map <- function(rst, year, vari, rela = FALSE, lims = NULL, dout = './png/anom'){
  
  ## 
  # rst <- rast(fles$fles[1])
  # year <- 2006
  # vari <- 'chl'
  # lims <- NULL
  # rela <- FALSE
  # dout<- './png/anom'
  
  ## To start
  dtes <- get.ym(rst)
  stopifnot('El año no está en el raster' = year %in% dtes$year)
  cat('Mapping: ', vari, ' ', year, '\n', sep = '')
  
  ## Palette, units and label of the variable
  ppal <- filter(pals, vari == .env$vari)
  if(nrow(ppal) == 0) ppal <- tibble(vari = vari, labl = vari, pall = 'RdBu', revs = TRUE, unts = '')
  colr <- brewer.pal(11, ppal$pall)
  if(ppal$revs) colr <- rev(colr)
  unts <- if(rela) '%' else ppal$unts
  
  ## Colour limits: if not given, from all the months and years of the raster (not only this year)
  if(is.null(lims)) lims <- get.lims(rst)
  
  ## Layers of the year to table (only sea cells)
  lyrs <- which(dtes$year == year)
  lyer <- rst[[lyrs]]
  names(lyer) <- glue('m{sprintf("%02d", dtes$mnth[lyrs])}')
  tble <- as.data.frame(lyer, xy = TRUE, na.rm = NA) |>
    pivot_longer(-c(x, y), names_to = 'mnth', values_to = 'valu') |>
    drop_na(valu) |>
    mutate(mnth = factor(mlbl[as.numeric(str_extract(mnth, '\\d{2}'))], levels = mlbl))
  
  ## Missing months (an incomplete year) stay as empty panels
  miss <- setdiff(1:12, dtes$mnth[lyrs])
  capt <- glue('Referencia {min(refp)}–{max(refp)} · escala común a todos los meses y años, saturada en ±{lims[2]} {unts}')
  capt
  if(length(miss) > 0) capt <- glue('{capt} · sin datos: {paste(mlbl[miss], collapse = ", ")}')
  
  ## Land cropped to the raster extent
  extn <- ext(rst) |> as.vector()
  bbox <- st_bbox(c(xmin = extn[[1]], ymin = extn[[3]], xmax = extn[[2]], ymax = extn[[4]]), crs = st_crs(land))
  lndc <- st_crop(land, bbox) |> suppressMessages()
  
  ## Scale bar and north arrow only in the last panel
  last <- tibble(mnth = factor('Diciembre', levels = mlbl))
  
  ## Map: 12 panels, 4 rows x 3 columns
  gmap <- ggplot() +
    geom_raster(data = tble, aes(x = x, y = y, fill = valu)) +
    geom_sf(data = lndc, fill = 'grey88', col = 'grey45', linewidth = 0.15) +
    scale_fill_gradientn(colours = colr, limits = lims, oob = scales::squish, name = unts) +
    facet_wrap(~ mnth, nrow = 4, ncol = 3, drop = FALSE) +
    coord_sf(xlim = extn[1:2], ylim = extn[3:4], expand = FALSE) +
    annotation_scale(data = last, location = 'bl', width_hint = 0.3, text_cex = 0.5) +
    annotation_north_arrow(data = last, location = 'tr', which_north = 'true', height = unit(0.6, 'cm'), width = unit(0.6, 'cm'), style = north_arrow_fancy_orienteering) +
    labs(title = glue('Anomalía mensual de {ppal$labl}, {year}'), caption = capt, x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 6),
      strip.background = element_blank(),
      strip.text = element_text(face = 'bold'),
      legend.key.height = unit(2, 'cm'),
      plot.caption = element_text(hjust = 0, colour = 'grey30')
    )
  
  ## To write the map (one folder per variable)
  sufx <- if(rela) 'anompct' else 'anom'
  dir_create(glue('{dout}/{vari}'))
  ggsave(glue('{dout}/{vari}/{vari}_{sufx}_{year}_ref{min(refp)}-{max(refp)}.png'), gmap, width = 8.5, height = 10, dpi = 300)
  
  return(gmap)
  
}

# To apply ----------------------------------------------------------------
maps <- map(seq_len(nrow(fles)), function(i){
  
  ## Monthly anomalies of the variable
  fle <- slice(fles, i)
  cat('To start the process: ', fle$name, '\n')
  rst <- rast(fle$fles)
  stopifnot('El raster debe estar en lon/lat' = is.lonlat(rst))
  
  ## A single colour scale for all the months and years of the variable
  lims <- get.lims(rst, qntl = qntl)
  
  
})

