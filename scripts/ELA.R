library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggdendro)
library(ggplot2)

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

find_energy_barrier = function(minima_df,states_energy_df, adj_matrix) {
  #vector with local minima
  min_list = minima_df$state
  
  #get all pairs of minima
  pairs = gtools::permutations(length(min_list), 2, min_list, repeats.allowed = T)
  barier_df = as.data.frame(pairs)
  barier_df$saddle = NA
  barier_df$barrier = NA
  
  #find shortest path between the pairs
  #selected 1
  for (i in 1:nrow(pairs)){
    states_df = states_energy_df
    A = adj_matrix
    pick = pairs[i,]
    network = igraph::graph_from_adjacency_matrix(A, "undirected")
    path = igraph::shortest_paths(network, from = pick[1], to = pick[2])$vpath[[1]]$name
    #find highest energy on the path
    maxE = max(states_df$energy[states_df$state %in% path])
    if (pick[1]==pick[2]){
      barier_df$saddle[i] = 0
      barier_df$barrier[i] = 0
    }
    else{
      while(length(path)!=0){
        maxE_last = maxE
        max_state = states_df$state[states_df$energy==maxE]
        states_df = states_df[states_df$energy<maxE,]
        A = get_adjacency(states_df$state)
        network = igraph::graph_from_adjacency_matrix(A, "undirected")
        path = try(igraph::shortest_paths(network, from = pick[1], to = pick[2])$vpath[[1]]$name,silent = T)
        maxE = try(max(states_df$energy[states_df$state %in% path]),silent = T)
      }
      barier_df$saddle[i] = max_state
      barier_df$barrier[i] = maxE_last
    }
  }
  colnames(barier_df)[1:2] = c("min1","min2")
  barier_df = barier_df %>% filter(barrier!=0)
  return(barier_df)
}

plot_disconnectivity = function(barrier_df) {
  min_list = unique(barrier_df$min1)
  min_n = length(min_list)
  #bariers as matrix
  barrier_mat = matrix(barrier_df$barrier, nrow = min_n,ncol=min_n)
  dist = abs(barrier_mat*lower.tri(barrier_mat, diag = T))
  dist = as.dist(dist)
  hc = hclust(dist)
  
  #make custom dendogram
  disc_graph = ggdendrogram(hc)+
    theme_few()+
    labs(x = NULL,
         y = "Energy barrier")
  
  #create matrix showing the activation patterns
  patterns = as.numeric(unlist(strsplit(min_list,"")))
  activations = matrix(patterns,nrow = length(patterns)/length(min_list), ncol = length(min_list)) %>% reshape2::melt() %>% mutate(Var1 = as.factor(Var1),Var2 = as.factor(Var2),value=as.factor(value))
  pattern_plot = ggplot(data = activations, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color="white")+
    theme_minimal()+
    coord_fixed()+
    labs(y = "Minima")+
    scale_fill_manual(values=c("#F2F2F2","black"))+
    guides(fill=F)
  
  library(patchwork)
  return(disc_graph/pattern_plot)
}

simulate_transitions = function(n_iter, adjacency_mat, states_df){
  #number of transitions
  n_iter = n_iter+2000
  
  #vector for storing visited states
  visited = rep(0,n_iter)
  
  #random start
  start = sample(1:nrow(adjacency_mat),1)
  current = colnames(adjacency_mat)[start]
  
  for (i in 1:n_iter){
    #store current position
    visited[i] = current
    
    E_current = states_df$energy[states_df$state==current]
    
    #get neighbours
    neighbours = rownames(adjacency_mat)[adjacency_mat[,current]==1]
    
    #select one of the neighbours with P=1/N
    proposal = sample(neighbours,1)
    #get E
    E_proposal = states_df$energy[states_df$state==proposal]
    
    if (E_proposal>E_current){
      p_move = exp(E_current-E_proposal)
      random = runif(1)
      current = ifelse(p_move>random, proposal, current)
    } else current = proposal
  }
  visited = visited[-1:-2000]
  return(visited)
}