# Libraries --------------------------------------------------------------------

# Paths ----
# local project location
path <- "C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/"

# high performance computing cluster location
HPCpath <- "M:/Kelly/massvacc/"

# Packages ----
library(malariasimulation) # devtools::install_github('mrc-ide/malariasimulation', force=T)
library(cali)              # devtools::install_github('https://github.com/mrc-ide/cali')
library(netz)              # devtools::install_github('https://github.com/mrc-ide/netz')
library(tidyverse)
library(data.table)

# HPC
library(didehpc)

# to edit HPC username and password
# usethis::edit_r_environ

# If you are having trouble installing libraries on the HPC,
# manually install the dependencies, so that the pkgdepends does not get confused:

# ctx <- context::context_save(path = paste0(HPCpath, "contexts"), packages = c("dplyr", "malariasimulation","cali"))
# share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")
# config <- didehpc::didehpc_config(credentials = list(
#   username = Sys.getenv("DIDE_USERNAME"),
#   password = Sys.getenv("DIDE_PASSWORD")),
#   shares = share)

# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "later")
# obj$install_packages("mrc-ide/cali")
# obj$install_packages("mrc-ide/individual")
# obj$install_packages("mrc-ide/malariaEquilibrium")
# obj$install_packages("mrc-ide/malariasimulation")

# plotting
library(patchwork)
library(LaCroixColoR)      # devtools::install_github('johannesbjork/LaCroixColoR')
library(ggplot2)

# color palette
# lacroix_palette(type = "paired")

