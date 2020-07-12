library(tictoc) #measuring runtime
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gtools)

binarize = function(data){
  return(apply(data, 1, function(x) ifelse(x >= mean(x), 1, -1)) %>% t())
}

pMEM_PL = function(data, lr_rate = 0.08, iter_max = 1e5, stopping = 1e-9) {
  len = ncol(data)
  N = nrow(data)
  mean_emp = rowMeans(data)
  cor_emp = (data%*%t(data))/len
  
  h = matrix(0, nrow = N, ncol = 1)
  J = matrix(0, nrow = N, ncol = N)
  ones = matrix(1, nrow = 1, ncol = len)
  
  for (i in 1:iter_max) {
    likelihood_h = -mean_emp + rowMeans(tanh(J%*%data + h%*%ones))
    h = h - lr_rate*likelihood_h
    
    likelihood_J = -cor_emp + 0.5*data %*% t((tanh(J%*%data + h%*%ones)))/len+0.5*t((data%*%t((tanh(J%*%data + h%*%ones)))))/len
    likelihood_J = likelihood_J - diag(diag(likelihood_J))
    J = J - lr_rate*likelihood_J
    
    if (sqrt(norm(likelihood_J, type = "F")^2 + norm(likelihood_h, type="2")^2)/N/(N+1) < stopping) {break
    }
  }
  cat("\nIterations:",i)
  return(list(h=h,J=J))
}

get_states = function(N, type = c("1","0")) {
  if (type == "1") {return(t(gtools::permutations(2,N,c(1,-1), repeats.allowed = T)))}
  if (type == "0") {return(t(gtools::permutations(2,N,c(1,0), repeats.allowed = T)))}
}

get_probability = function(parameter_list) {
  h = parameter_list$h
  J = parameter_list$J
  
  N_ROI = nrow(h)
  states = get_states(N_ROI, type = "1")
  states_len = ncol(states)
  states_merged = states %>% t() %>% as.data.frame() %>% unite("state", sep = "")
  
  Z1 = -diag(t(0.5*J%*%states) %*% states)
  Z2 = rowSums(t((h %*% matrix(1,1,states_len)) * states))
  
  Z = sum(exp(-(Z1-Z2)))
  probs = exp(-(Z1-Z2))/Z
  probs = probs[order((states_merged$state))]
  return(probs)
}

get_empirical_prob = function(data){
  unique = t(data) %>%
    as.data.frame() %>%
    unite("state", sep = "") %>%
    group_by(state) %>%
    summarise(count = n())
  
  states_emp = unique$state
  states = get_states(nrow(data),"1") %>% t() %>% as.data.frame() %>% unite("state", sep = "") 
  
  add_states = data.frame(state=states[which(!(states$state %in% states_emp)),])
  add_states$count = rep(1,nrow(add_states))
  
  unique_new = rbind(unique,add_states)
  unique_new = unique_new[order(unique_new$state),]
  
  probs_emp = unique_new$count/sum(unique_new$count)
  
  return(probs_emp)
}

get_energy = function(parameter_list) {
  h = parameter_list$h
  J = parameter_list$J
  
  N_ROI = nrow(h)
  states = get_states(N_ROI, type = "1")
  states_len = ncol(states)
  
  Z1 = -diag(t(0.5*J%*%states) %*% states)
  Z2 = rowSums(t((h %*% matrix(1,1,states_len)) * states))
  
  return((Z1-Z2))
}

get_R2 = function(empirical, predicted){
  loglog = data.frame(empirical = empirical, predicted = predicted)
  
  m1 = lm(empirical ~ predicted, loglog)
  return(summary(m1)$r.squared)
}

pMEM_eval = function(data,parameters){
  len = ncol(data)
  N = nrow(data)
  
  #get all states
  states = get_states(N, "1")
  N_states = 2^N
  
  #Compute empirical probability and its entropy
  prob_emp = as.matrix(get_empirical_prob(data))
  
  Ent_emp = sum(-prob_emp * log2(prob_emp))
  
  #Compute first-order model entropy (only h)
  active = rowMeans(data==1)
  inactive = 1-active
  
  active_mat = matrix(active, N, N_states)
  inactive_mat = matrix(inactive, N, N_states)
  
  prob_1 = ((states+1)/2) * active_mat + ((1-states)/2) * inactive_mat
  prob_1 = as.matrix(Rfast::colprods(prob_1))
  
  Ent_1 = sum(-prob_1 * log2(prob_1))
  
  #Compute entropy of fitted pMEM
  pMEM_prob = get_probability(parameters)
  Ent_2 = sum(-pMEM_prob * log2(pMEM_prob))
  
  #Accuracy index - r
  r = round((Ent_1-Ent_2)/(Ent_1-Ent_emp),3)
  
  
  #Accuracy with KL divergence
  D1_mat = prob_emp * log2(prob_emp/prob_1)
  D1 = sum(D1_mat)
  
  D2_mat = prob_emp * log2(prob_emp/pMEM_prob)
  D2 = sum(D2_mat)
  
  d = round((D1-D2)/D1,3)
  
  #R2
  R2 = get_R2(prob_emp,pMEM_prob)
  
  result = list(r = d, R2 = R2)
  
  return(result)
}

run_pMEM = function(data){
  df = binarize(data)
  params = pMEM_PL(df)
  eval = pMEM_eval(df,params)
  
  result = list(parameters = params,
                metrics = eval)
  return(result)
}

run_folder = function(wave_band, selection, selection_name){
  input = paste0("data/clean_", wave_band)
  files = list.files(input, full.names = T)
  files = files[c(-29,-30)]
  
  pMEM_results = matrix(NA, length(files), 5)
  colnames(pMEM_results) = c("ID", "diagnosis", "condition", "r", "R2") #29
  for (i in 1:length(files)){
    load(files[i])
    data = data[,selection]
    data = t(data)
    results = run_pMEM(data)
    save_to = paste0("data/pMEM_", selection_name, "/",wave_band, "/",strsplit(files[i],"/")[[1]][3])
    save(results, file = save_to)
    
    info = strsplit(strsplit(files[i],"/")[[1]][3],"_")[[1]]
    ID = info[1]
    diagnosis = info[2]
    condition = strsplit(info[3],"\\.")[[1]][1]
    
    results_row = cbind(ID, diagnosis, condition, results$metrics$r, results$metrics$R2)
    
    pMEM_results[i,] = results_row
  }
  
  pMEM_results = as.data.frame(pMEM_results)
  pMEM_results = na.omit(pMEM_results)
  
  output = paste0("results/pMEM_metrics/",selection_name, "_", wave_band, ".csv")
  readr::write_csv(pMEM_results, output)
}
