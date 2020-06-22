library(tidyverse)
library(tidygraph)
library(ggraph)

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

basin_size_estimator = function(adjancency, states_energy){
  basin_membership = rep(NA, ncol(adjancency))
  positions = data.frame(from=NA,to =NA)
  
  iter = 1e8 #random large number
  
  for (i in 1:ncol(adjancency)){
    start = colnames(adjancency)[i]
    for (j in 1:iter){
      start_E = states_energy$energy[states_energy$state==start]
      neighbours = rownames(adjancency)[adjancency[,start]==1]
      neighbours_E =  states_energy$energy[states_energy$state %in% neighbours]
      
      states = c(start,neighbours)
      E = c(start_E,neighbours_E)
      new = states[order(E)[1]]
      positions[nrow(positions)+1,] = c(start,new)
      
      if (start==new){break}
      else {start=new}
    }
    basin_membership[i] = start
  }
  
  basin_membership = data.frame(state = colnames(adjancency),
                                basin = basin_membership)
  result = list(membership = basin_membership,positions = positions)
  return(result)
}

BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

plot_basin = function(basin_df) {
  pos_dec = basin_df$positions[-1,] %>% distinct()
  mem_dec = basin_df$membership
  
  pos_dec$from = sapply(pos_dec$from,BinToDec)
  pos_dec$to = sapply(pos_dec$to,BinToDec)
  
  mem_dec$state = sapply(mem_dec$state,BinToDec)
  mem_dec$basin = sapply(mem_dec$basin,BinToDec)
  colnames(mem_dec)[1] = "from"
  
  basin_network = pos_dec %>%
    full_join(mem_dec, by = "from")
  basin_network$basin = as.factor(basin_network$basin)
  
  tidy = as_tbl_graph(basin_network) 
  
  plot = tidy %>%
    activate(nodes) %>%
    mutate(color = .E()$basin) %>%
    activate(edges) %>%
    ggraph(layout="stress")+
    geom_edge_link()+
    geom_node_point(aes(colour=color),size=1)+
    scale_color_tableau()+
    guides(color=guide_legend(title = "Local minimum"))+
    theme_few()
  return(plot)
}