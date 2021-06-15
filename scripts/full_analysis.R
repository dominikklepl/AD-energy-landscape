run = function(name, window, down, base){
  source("scripts/pMEM.R")
  source("scripts/ELA.R")
  library(tictoc)

  root = paste0(base, name)
  fold_pMEM = paste0(root, "/pMEM")
  fold_pMEM_feat = paste0(fold_pMEM, "_features")
  dir.create(base)
  dir.create(root)
  dir.create(fold_pMEM)
  dir.create(fold_pMEM_feat)
  
  downsample = function(x, freq){
    len = length(x)/2000
    
    n = floor(len*freq)
    return(eegkit::eegresample(x, n))
  }
  
  run_folder = function(wave_band, selection, size, down){
    input = paste0("data/clean_", wave_band)
    files = list.files(input, full.names = T)
    files = files[c(-29,-30)]
    
    pMEM_results = matrix(NA, length(files), 7)
    colnames(pMEM_results) = c("ID", "diagnosis", "condition", "epoch", "band", "r", "R2") #29
    for(i in 1:length(files)){
      load(files[i])
      data = data[,selection]
      
      #downsample
      data = apply(data, 2, downsample, freq=down) %>% as.data.frame()
      
      results = run_pMEM(data, size)
      
      info = strsplit(strsplit(files[i],"/")[[1]][3],"_")[[1]]
      ID = info[1]
      diagnosis = info[2]
      condition = info[3]
      epoch = strsplit(info[4],"\\.")[[1]][1]
      
      #save to data/pMEM_3_500Hz
      filename = paste0(wave_band,"_",strsplit(files[i],"/")[[1]][3])
      save_to = paste0(fold_pMEM,"/",filename)
      save(results, file = save_to)
      
      results_row = cbind(ID, diagnosis, condition, epoch, wave_band, results$metrics$r, results$metrics$R2)
      
      pMEM_results[i,] = results_row
    }
    
    pMEM_results = as.data.frame(pMEM_results)
    pMEM_results = na.omit(pMEM_results)
    
    return(pMEM_results)
  }
  
  
  selections = readr::read_csv2("channel_selection/final_selection.csv")
  bands = c("alpha", "beta", "delta", "full", "gamma", "theta")
  
  registerDoParallel(6)
  result = foreach (b = bands, .combine = "rbind", .packages = "tidyverse") %dopar% {
    source("scripts/pMEM.R")
    selection = selections %>% dplyr::filter(band == b)
    selection = as.character(selection[1,1:10]) %>% na.omit()
    run_folder(wave_band = b, selection = selection, size=window, down = down)
  }
  stopImplicitCluster()
  readr::write_csv(result, paste0(root,"/pMEM_metrics.csv"))
  
  #collect J and h parameters
  files = list.files(fold_pMEM, full.names = T)
  bands = c("alpha", "beta", "delta", "full", "gamma", "theta")
  
  registerDoParallel(6)
  for(b in bands){
    Jh_features = foreach(f=files[grep(b,files)], .combine = "rbind", .packages = "tidyverse") %dopar% {
      load(f)
      J = results$parameters$J
      J[lower.tri(J)] = 0
      diag(J) = results$parameters$h + 1e-50
      J = reshape2::melt(J) %>% filter(value != 0) %>% unite("name",Var1:Var2)
      res = J$value
      names(res) = J$name
      res = as.data.frame(t(res))
      
      res$band = strsplit(f, "/")[[1]][4] %>% strsplit("_") %>% .[[1]] %>% .[1]
      res$ID = strsplit(f, "/")[[1]][4] %>% strsplit("_") %>% .[[1]] %>% .[2]
      res$diagnosis = strsplit(f, "/")[[1]][4] %>% strsplit("_") %>% .[[1]] %>% .[3]
      res$condition = strsplit(f, "/")[[1]][4] %>% strsplit("_") %>% .[[1]] %>% .[4]
      res$epoch = strsplit(f, "/")[[1]][4] %>% strsplit("_") %>% .[[1]] %>% .[5] %>% strsplit("\\.") %>% .[[1]] %>% .[1]
      res
    }
    
    save_to = paste0(fold_pMEM_feat,"/", b, ".csv")
    readr::write_csv(Jh_features, save_to)
    
    print(b)
  }
  stopImplicitCluster()
  
  
  #statistic pMEM features
  files = list.files(fold_pMEM_feat, full.names = T)
  
  shannon = function(x){
    dens = density(x)$y
    return(entropy::entropy.empirical(dens))
  }
  
  features = foreach(f = files, .combine = "rbind") %do% {
    band = strsplit(f, "/")[[1]][4] %>% strsplit("\\.") %>% .[[1]] %>% .[1]
    df = read.csv(f)
    
    #feature engineering
    J = df[,1:55]
    J = J[,c(-1,-3,-6,-10,-15,-21,-28,-36,-45,-55)]
    h = df[,c(1,3,6,10,15,21,28,36,45,55)]
    
    h$mean = apply(h[,1:10], 1, mean)
    h$sd = apply(h[,1:10], 1, sd)
    h$IQR = apply(h[,1:10], 1, IQR)
    h$entropy = apply(h[,1:10],1, shannon)
    colnames(h)[11:14] = paste0(colnames(h)[11:14],"_h")
    
    J$mean = apply(J[,1:45], 1, mean)
    J$sd = apply(J[,1:45], 1, sd)
    J$IQR = apply(J[,1:45], 1, IQR)
    J$entropy = apply(J[,1:45],1, shannon)
    colnames(J)[46:49] = paste0(colnames(J)[46:49],"_J")
    
    features = cbind(df[,56:60], h[,11:14], J[,46:49])
    features$band = band
    features
  }
  
  readr::write_csv(features,paste0(fold_pMEM_feat, "/statistic_features.csv"))
  
  #ML on connectivity
  library(caret)
  library(foreach)
  library(doParallel)
  
  files = list.files(fold_pMEM_feat, full.names = T)
  files = files[-grep("statistic",files)]
  
  fit_ctrl = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          verboseIter = F,
                          classProbs = T,
                          summaryFunction = twoClassSummary,
                          savePredictions = "final",
                          allowParallel = T)
  
  registerDoParallel(10)
  
  results = foreach(f = files, .combine = "rbind") %do% {
    band = strsplit(f, "/")[[1]][4] %>% strsplit("\\.") %>% .[[1]] %>% .[1]
    df = read.csv(f)
    EC = df %>% filter(condition=="EC") %>% select(-band, -condition, -epoch, -ID)
    EO = df %>% filter(condition=="EO") %>% select(-band, -condition, -epoch, -ID)
    
    model = train(diagnosis ~ .,
                  data = EC,
                  trControl = fit_ctrl,
                  method = 'svmRadial',
                  preProcess = c("scale", "center", "pca"),
                  metric = c("ROC"))
    
    result_EC = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EC$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EC$condition = "EC"
    
    model = train(diagnosis ~ .,
                  data = EO,
                  trControl = fit_ctrl,
                  method = 'svmRadial',
                  preProcess = c("scale", "center", "pca"),
                  metric = c("ROC"))
    
    result_EO = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EO$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EO$condition = "EO"
    
    result = rbind(result_EC, result_EO)
    result$band = band
    result
  }
  readr::write_csv(results, paste0(root,"/ML_connectivity.csv"))
  stopImplicitCluster()
  
  #ML on stats of connectivity
  data = read.csv(paste0(fold_pMEM_feat,"/statistic_features.csv"))
  bands = unique(data$band)
  
  grid = expand.grid(C = c(0.25,0.5,1),
                     sigma = seq(0,1,length.out = 20))
  
  fit_ctrl = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          verboseIter = F,
                          classProbs = T,
                          summaryFunction = twoClassSummary,
                          savePredictions = "final",
                          allowParallel = T)
  
  registerDoParallel(10)
  results = foreach(b = bands, .combine = "rbind") %do% {
    df = data %>% filter(band==b)
    
    EC = df %>% filter(condition=="EC") %>% select(-band, -condition, -epoch, -ID)
    EO = df %>% filter(condition=="EO") %>% select(-band, -condition, -epoch, -ID)
    
    model = train(diagnosis ~ .,
                  data = EC,
                  trControl = fit_ctrl,
                  method = 'svmRadial',
                  preProcess = c("scale", "center","pca"),
                  metric = "ROC",
                  tuneGrid = grid)
    
    result_EC = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EC$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EC$condition = "EC"
    
    model = train(diagnosis ~ .,
                  data = EO,
                  trControl = fit_ctrl,
                  method = 'svmRadial',
                  preProcess = c("scale", "center","pca"),
                  metric = c("ROC"))
    
    result_EO = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EO$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EO$condition = "EO"
    
    result = rbind(result_EC, result_EO)
    result$band = b
    result
  }
  stopImplicitCluster()
  readr::write_csv(results, paste0(root,"/ML_connectivity_stats.csv"))
  
  #extract energy features
  source("scripts/ELA.R")
  library(foreach)
  library(doParallel)
  
  analyze_energy_landscape = function(file){
    source("scripts/ELA.R")
    load(file)
    
    fname = strsplit(file, "/")[[1]][4]
    band = strsplit(fname, "_")[[1]][1]
    ID = strsplit(fname, "_")[[1]][2]
    diagnosis = strsplit(fname, "_")[[1]][3]
    condition = strsplit(fname, "_")[[1]][4]
    epoch = strsplit(strsplit(fname, "_")[[1]][5], "\\.")[[1]][1]
    
    info = data.frame(band, ID, diagnosis, condition, epoch)
    
    
    #data.frame for storing results
    features = data.frame(minima = NA, E_diff = NA, basin_sd = NA, duration = NA)
    
    states_energy = get_state_energy(results$parameters)
    adj = get_adjacency(states_energy$state)
    
    #local minima related features
    minima = find_minima(adj, states_energy)
    features$minima = nrow(minima)
    
    features$E_diff = mean_energy_difference(minima)
    if(is.na(features$E_diff)){features$E_diff=0}
    
    global_min = minima$state[minima$energy==min(minima$energy)]
    
    #basin size
    basin_df = basin_size_estimator(adj, states_energy)
    #basin_plot = plot_basin(basin_df)
    features$basin_sd = basin_df$membership %>% 
      group_by(basin) %>%
      summarise(count = n()) %>%
      pull() %>%
      sd()
    
    if(is.na(features$basin_sd)){features$basin_sd= nrow(basin_df$membership)}
    
    #simulate transitions
    duration = NULL
    while(is_empty(duration)){
      x = simulate_transitions(2e4, adj, states_energy)
      duration = data.frame(state = x) %>%
        left_join(basin_df$membership, by="state") %>%
        group_by(basin) %>%
        summarise(duration = n()) %>%
        mutate(duration = duration/500)
      duration = duration$duration[duration$basin==global_min]
    }
    features$duration = duration
    
    features = cbind(info, features)
    return(features)
  }
  
  
  files = list.files(fold_pMEM, full.names = T)
  
  registerDoParallel(10)
  tic()
  results = foreach(f = files, .combine = "rbind", .packages = c("tidyverse")) %dopar% {
    feat = try(analyze_energy_landscape(f))
    
    if(is.character(feat)){
      rep(NA,(9))
    } else {feat}
  }
  toc()
  stopImplicitCluster()
  
  nas = which(is.na(results$band))
  
  for (n in nas){
    results[n,] = analyze_energy_landscape(files[n])
  }
  
  readr::write_csv(results, paste0(root,"/energy_features.csv"))
  
  #test for differences
  fit_model = function(data){
    model = lmer(value ~ diagnosis + (1|ID), data)
    summary = as.data.frame(summary(model)$coefficients)
    
    result = data.frame(band=data$band[1],
                        condition=data$condition[1],
                        feature=data$feature[1],
                        intercept=summary$Estimate[1],
                        error_i=summary$`Std. Error`[1],
                        B=summary$Estimate[2],
                        error_B=summary$`Std. Error`[2],
                        t=summary$`t value`[2],
                        p=summary$`Pr(>|t|)`[2]
    )
    return(result)
  }
  
  library(dplyr)
  library(lmerTest)
  library(foreach)
  
  df = read.csv(paste0(root,"/energy_features.csv"))
  df = df %>% reshape2::melt(id=1:5)
  colnames(df)[6] = "feature"
  
  bands = unique(df$band)
  conditions = unique(df$condition)
  features = unique(df$feature)
  
  results = data.frame()
  for(b in bands){
    for(c in conditions){
      for(f in features){
        res = df %>% filter(band==b, condition==c, feature==f) %>% fit_model()
        results = rbind(results, res)
      }
    }
  }
  
  readr::write_csv(results, paste0(root,"/lineardiff_full.csv"))
  
  #get only significant diffs
  signif = results %>% filter(p<=0.05)
  readr::write_csv(signif, paste0(root,"/lineardiff_significant.csv"))
  
  #ML on energy
  library(caret)
  library(foreach)
  library(doParallel)
  
  data = read.csv(paste0(root,"/energy_features.csv"))
  data$diagnosis = as.factor(data$diagnosis)
  bands = unique(data$band)
  
  
  fit_ctrl = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          verboseIter = F,
                          classProbs = T,
                          summaryFunction = twoClassSummary,
                          savePredictions = "final",
                          allowParallel = T)
  
  registerDoParallel(10)
  
  method = 'svmPoly'
  results = foreach(b = bands, .combine = "rbind") %do% {
    df = data %>% filter(band==b)
    
    EC = df %>% filter(condition=="EC") %>% select(-band, -condition, -epoch, -ID)
    EO = df %>% filter(condition=="EO") %>% select(-band, -condition, -epoch, -ID)
    
    model = train(diagnosis ~ .,
                  data = EC,
                  trControl = fit_ctrl,
                  method = method,
                  preProcess = c("scale", "center"),
                  metric = "ROC")
    
    result_EC = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EC$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EC$condition = "EC"
    
    model = train(diagnosis ~ .,
                  data = EO,
                  trControl = fit_ctrl,
                  method = method,
                  preProcess = c("scale", "center"),
                  metric = c("ROC"))
    
    result_EO = model$resample %>%
      summarise(meanROC = mean(ROC),sdROC = sd(ROC),
                meanSens = mean(Sens), sdSens = sd(Sens),
                meanSpec = mean(Spec), sdSpec = sd(Spec))
    result_EO$Accuracy = confusionMatrix(model$pred$pred, model$pred$obs)$overall[1]
    result_EO$condition = "EO"
    
    result = rbind(result_EC, result_EO)
    result$band = b
    result
  }
  stopImplicitCluster()
  readr::write_csv(results,paste0(root,"/ML_energy.csv"))
}
