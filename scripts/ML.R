library(caret)
library(dplyr)
library(kernlab)

train_models = function(features, data, algorithm, cv){
  fit_ctrl = trainControl(method = cv,
                          number = 10,
                          repeats = 10)
  results = matrix(NA, length(features), 3)
  for (i in 1:length(features)){
    df = data %>% select(diagnosis, features[i])
    model = train(diagnosis ~ .,
                  data = df,
                  trControl = fit_ctrl,
                  method = algorithm)
    results[i,] = c(features[i], model$results$Accuracy, model$results$AccuracySD)
  }
  colnames(results) = c("feature", "accuracy", "sd")
  results = as.data.frame(results)
  return(results)
}

train_file = function(file, algorithm, cv){
  split = strsplit(strsplit(file,"/")[[1]][3],"_")[[1]]
  selection = split[1]
  band = strsplit(split[2], "\\.")[[1]][1]
  
  data = read.csv(file)
  EC = data %>% filter(condition == "EC")
  EC = EC %>% select(-condition,-ID)
  EO = data %>% filter(condition == "EO")
  EO = EO %>% select(-condition,-ID)
  
  result_EC = train_models(features, EC, algorithm, cv)
  result_EO = train_models(features, EO, algorithm, cv)
  
  result_EC$condition = "EC"
  result_EO$condition = "EO"
  
  result = rbind(result_EC, result_EO)
  result$selection = selection
  result$band = band
  
  return(result)
}
