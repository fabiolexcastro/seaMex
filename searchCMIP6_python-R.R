

install.packages(c('pacman', 'reticulate', 'terra', 'tidyverse', 'glue', 'fs'))
packageVersion('reticulate')   # >= 1.41

# "C:/scripts/.venv_copernicus/Scripts/python.exe"

esgf <- reticulate::import('intake_esgf')  
py_config()          