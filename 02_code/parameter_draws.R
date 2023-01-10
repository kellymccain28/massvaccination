# Set up draws of malariasimulation parameters, pulled from MalariaLaunchR
# 50 draws are used for uncertainty runs

source("./02_code/data_libraries.R")


# get draws --------------------------------------------------------------------
# pull all .txt files from MalariaLaunchR folder and combine
# https://github.com/mrc-ide/MalariaLaunchR/tree/master/inst/model_files/parameters
files <- list.files(path = paste0(path, "01_data/parameter_draws/"), pattern = "*.txt", full.names = TRUE)
dat_list <- lapply(files, function (x) read.delim(x, header = F, col.names = c('var', 'value')))

dat <- rbindlist(dat_list, fill = TRUE, idcol = "file")

d <- dat |> group_by(file) |> pivot_wider(names_from = 'var', values_from = 'value')
head(d); str(d)

# save to local machine and HPC access drive
saveRDS(d, paste0(path, "03_output/parameter_draws.rds"))
# saveRDS(d, paste0(HPCpath, "parameter_draws.rds"))
