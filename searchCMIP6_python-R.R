

install.packages(c('pacman', 'reticulate', 'terra', 'tidyverse', 'glue', 'fs'))
packageVersion('reticulate')   # >= 1.41

library(reticulate)

# "C:/scripts/.venv_copernicus/Scripts/python.exe"

Sys.setenv(RETICULATE_PYTHON = 'managed')   # force the uv environment of py_require()

reticulate::py_require('intake-esgf')
esgf <- reticulate::import('intake_esgf')
py_config()        
