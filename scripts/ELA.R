library(tidyverse)

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

get_states = function(N, type = c("1","0")) {
  if (type == "1") {return(t(gtools::permutations(2,N,c(1,-1), repeats.allowed = T)))}
  if (type == "0") {return(t(gtools::permutations(2,N,c(1,0), repeats.allowed = T)))}
}

get_state_energy = function(pMEM){
  energy = get_energy(pMEM)
  states = get_states(length(pMEM$h),"0") %>% t() %>%
    as.data.frame() %>%
    unite("state", sep = "")
  
  states_energy = data.frame(state = states$state,
                             energy = energy)
  return(states_energy)
}

get_adjacency = function(states_list){
  len = length(states_list)
  adj_df = gtools::permutations(len,2, states_list, repeats.allowed = T) %>%
    as.data.frame() %>%
    mutate(link = ifelse(stringdist::stringdist(V1,V2, method = "hamming")==1,1,0))
  adj_mat = matrix(adj_df$link, nrow = len,ncol = len)
  colnames(adj_mat) = states_list
  rownames(adj_mat) = states_list
  return(adj_mat)
}

is_minimum = function(state,energy=states_energy,A=adj){
  one = A[,state]==1
  neighbours = rownames(A)[one]
  
  state_E = energy$energy[energy$state==state]
  rest_E = energy$energy[energy$state %in% neighbours]
  
  return(ifelse(min(rest_E) > state_E, T,F))
}

find_minima = function(adj, states_energy){
  minima = sapply(colnames(adj),is_minimum)
  minima = colnames(adj)[minima]
  minima = data.frame(state = minima, energy = states_energy$energy[states_energy$state %in% minima])
  
  return(minima)
}