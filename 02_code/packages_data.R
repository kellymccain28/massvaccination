# Loading packages and data and initial parameters 

library(malariasimulation)  
library(ggplot2)
library(mgcv)
library(malariaEquilibrium)
library(reshape2)
library(DiagrammeR)
library(cowplot)
library(dplyr)
library(tidyr)
library(patchwork)
library(foresite)
library(cali) #remotes::install_github("mrc-ide/cali")
library(data.table)
library(didehpc)
library(LaCroixColoR)# devtools::install_github('johannesbjork/LaCroixColoR')
library(purrr)
library(flextable)
library(gtsummary)
library(janitor)
library(stringr)
library(postie) #remotes::install_github("mrc-ide/postie")
library(ggforce)
library(forcats)
library(RColorBrewer)


knitr::opts_knit$set(message=FALSE, warning=FALSE)

`young children` <- 'lightblue'
`under 5s` <- 'maroon'
`school-aged` <- 'skyblue4'
`all children` <- 'palegreen4'
everyone <- 'tomato3'
        
col_SV_models <- list(`young children`, `under 5s`, `school-aged`, `all children`,everyone)
        
set.seed(1234)

col_template <- c('#8A7967', '#21677E', '#822F5A', '#B5D334', '#6A3B77', '#D07232', '#FFDB00')

# Paths ----
# local project location
path <- "C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Mass vaccination/massvaccination/"

# high performance computing cluster location
HPCpath <- "M:/Kelly/massvacc/"


# Setting the configuration
# share <- didehpc::path_mapping("malaria", "M:", "//fi--didenas1/malaria", "M:")
# share <- didehpc::path_mapping("PriorityPathogens", "P:", "//projects.dide.ic.ac.uk/PriorityPathogens", "P:")
# config <- didehpc::didehpc_config(credentials = list(
#   username = Sys.getenv("DIDE_USERNAME"),
#   password = Sys.getenv("DIDE_PASSWORD")),
#   shares = share)
# 
# drat:::add("mrc-ide")
# root<- paste0(HPCpath, "contexts/")
# 
# packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation', 'cali')
# src <- conan::conan_sources("github::mrc-ide/malariasimulation")
# 
# # save a context--------------------------------------------------------------
# ctx <- context::context_save(path = root, 
#                              packages = packages, 
#                              package_sources = src)
# 
# ## Manually install the dependencies, so that the pkgdepends does not
# ## get confused:--------------------------------------------------------------
# obj <- didehpc::queue_didehpc(ctx, config = config, provision = "later")
# obj$install_packages("mrc-ide/cali")
# obj$install_packages("mrc-ide/individual")
# obj$install_packages("mrc-ide/malariaEquilibrium")
# obj$install_packages("mrc-ide/malariasimulation")

