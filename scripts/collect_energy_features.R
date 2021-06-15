library(foreach)
library(doParallel)

collect_energy = function(path, cores){
  new_path = paste0(path, "/pMEM")
  files = list.files(new_path, full.names = T)
  registerDoParallel(cores)
  features = foreach(f=files, .combine = "rbind", .packages = "tidyverse")%dopar%{
    source("scripts/ELA.R")
    load(f)
    name = strsplit(f, "/")[[1]] %>% .[grep("RData",.)] %>% strsplit("_") %>% .[[1]]
    info = data.frame(band = name[1], ID = name[2],
                      diagnosis = name[3], condition = name[4],
                      epoch = strsplit(name[5], "\\.")[[1]][1])
    params = results$parameters
    row = get_state_energy(params)$energy %>% t()
    cbind(info, row)
  }
  stopImplicitCluster()
  return(features)
}
