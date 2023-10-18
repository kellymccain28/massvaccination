# Loading packages and data and initial parameters 

library(malariasimulation)  #devtools::install_github("mrc-ide/malariasimulation@epic/hybrid_pev")
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
library(ggtext)


knitr::opts_knit$set(message=FALSE, warning=FALSE)

col515 <- c('#0E6251', '#148F77', '#1ABC9C', '#76D7C4')
col59 <- c('#4A235A', '#6C3483', '#A569BD', '#D2B4DE')
col5m15y <- c('#154360', '#1F618D', '#2980B9', '#7FB3D5')
col5m3y <- c('#7E5109', '#B9770E', '#F39C12', '#F8C471')
col5m5y <- c('#78281F','#B03A2E', '#E74C3C', '#F1948A')
colsAB <- c('#186A3B', '#239B56', '#58D68D', '#82E0AA')
colsSVhybrid <- c('#283747', '#85929E')
colors <- c(col515, col59, col5m15y, col5m3y, col5m5y, colsAB, colsSVhybrid)

CUcols <- c('#B03A2E','#6C3483','#1F618D','#148F77','#B9770E','#239B56','#283747', 'tan','#85929E')

        
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

