
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(reticulate, tidyverse, glue, fs)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Python ------------------------------------------------------------------
## Once, in the cmd with the same Python you use for copernicusmarine:
##   python -m pip install intake-esgf
py_config()
py_config()
py_exe()   

## Once: install intake-esgf in THAT Python
system2(py_exe(), c('-m', 'pip', 'install', '--upgrade', 'pip', 'setuptools'))
system2(py_exe(), c('-m', 'pip', 'install', 'intake-esgf'))

## Check and import
py_module_available('intake_esgf')   # must be TRUE

reticulate::py_require('intake-esgf')
esgf <- reticulate::import('intake_esgf')
py_config()

esgf <- reticulate::import('intake_esgf')

##

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(reticulate, tidyverse, glue, fs)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Python ------------------------------------------------------------------
## reticulate >= 1.41: declare the package; reticulate builds the environment with uv (no pip)
reticulate::py_require('intake-esgf')
esgf <- reticulate::import('intake_esgf')
py_config()

## Indices: the STAC ones have CMIP6 publishing errors (warning of intake-esgf) -> off
## The Globus ESGF2-US index stays on; CEDA and DKRZ (Solr) are added for European models
esgf$conf$set(indices = reticulate::dict('discovery.east.esgf.io' = F, 'discovery.west.esgf.io' = F,
                                         'esgf.ceda.ac.uk' = T, 'esgf-data.dkrz.de' = T))

# Parameters --------------------------------------------------------------
## Surface versions (2D, small) of thetao, so, ph, o2, chl
vars <- c('tos', 'sos', 'phos', 'o2os', 'chlos')
expr <- c('historical', 'ssp370')

# Function to search ------------------------------------------------------
srch.cmip <- function(vari, expr = c('historical', 'ssp370')){
  
  # vari <- 'chlos'; expr <- c('historical', 'ssp370')
  
  ## To start
  cat('To search: ', vari, '\n')
  ctl <- esgf$ESGFCatalog()
  
  ## Search: project = 'CMIP6' is mandatory (the default project of intake-esgf is CMIP7)
  rsl <- tryCatch({
    ctl$search(project = 'CMIP6', experiment_id = expr, table_id = 'Omon', variable_id = vari, quiet = T)
    ### pandas -> R through plain Python lists: py_get_attr() takes the raw pandas object, so
    ### the result does not depend on how each reticulate / pandas version converts a DataFrame
    col <- c('source_id', 'experiment_id', 'member_id', 'variable_id', 'grid_label')
    dfr <- reticulate::py_get_attr(ctl, 'df')$to_dict('list')
    map(col, ~ as.character(unlist(dfr[[.x]]))) %>% setNames(col) %>% as_tibble()
  }, error = function(e){
    cat('No results or error: ', conditionMessage(e), '\n')
    tibble()
  })
  
  cat('Datasets found: ', nrow(rsl), '\n')
  return(rsl)
  
}

# To apply the function ---------------------------------------------------
srch <- map_dfr(vars, srch.cmip, expr = expr)
stopifnot('Ninguna búsqueda devolvió resultados: revisa la conexión o los nombres' = nrow(srch) > 0)

## Per model and variable: the first member with historical AND ssp370 (gr preferred over gn)
avlb <- srch %>%
  distinct(source_id, variable_id, member_id, grid_label, experiment_id) %>%
  group_by(source_id, variable_id, member_id, grid_label) %>%
  filter(all(expr %in% experiment_id)) %>%
  ungroup() %>%
  distinct(source_id, variable_id, member_id, grid_label) %>%
  mutate(rnum = str_extract(member_id, '(?<=^r)\\d+') %>% as.numeric(),
         fnum = str_extract(member_id, '(?<=f)\\d+$') %>% as.numeric()) %>%
  arrange(source_id, variable_id, rnum, fnum, desc(grid_label == 'gr')) %>%
  group_by(source_id, variable_id) %>%
  slice(1) %>%
  ungroup()

## Table: model x variable -> member (grid)
tble <- avlb %>%
  mutate(vlue = as.character(glue('{member_id} ({grid_label})'))) %>%
  dplyr::select(source_id, variable_id, vlue) %>%
  pivot_wider(names_from = variable_id, values_from = vlue) %>%
  dplyr::select(source_id, any_of(vars)) %>%
  mutate(nvar = rowSums(!is.na(across(any_of(vars))))) %>%
  arrange(desc(nvar), source_id)

print(tble, n = Inf, width = Inf)

## To write the table
dir_create('./tbl')
write_csv(tble, './tbl/cmip6_omon_ssp370_availability.csv')
