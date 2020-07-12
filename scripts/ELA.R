library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggdendro)
library(ggplot2)
library(patchwork)

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

is_minimum = function(state,energy,A){
  one = A[,state]==1
  neighbours = rownames(A)[one]
  
  state_E = energy$energy[energy$state==state]
  rest_E = energy$energy[energy$state %in% neighbours]
  
  return(ifelse(min(rest_E) > state_E, T,F))
}

find_minima = function(adjacency, states_energy){
  minima = sapply(colnames(adjacency),is_minimum, energy=states_energy,A=adjacency)
  minima = colnames(adjacency)[minima]
  minima = data.frame(state = minima, energy = states_energy$energy[states_energy$state %in% minima])
  
  return(minima)
}

mean_energy_difference = function(minima_df){
  global = min(minima_df$energy)
  mean_diff = minima_df %>%
    mutate(diff = energy - min(energy)) %>%
    filter(diff != 0) %>%
    summarise(mean = mean(diff)) %>%
    as.numeric()
  return(mean_diff)
}

basin_size_estimator = function(adjacency, states_energy){
  basin_membership = rep(NA, ncol(adjacency))
  positions = data.frame(from=NA,to =NA)
  
  iter = 1e8 #random large number
  
  for (i in 1:ncol(adjacency)){
    start = colnames(adjacency)[i]
    for (j in 1:iter){
      start_E = states_energy$energy[states_energy$state==start]
      neighbours = rownames(adjacency)[adjacency[,start]==1]
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
  
  basin_membership = data.frame(state = colnames(adjacency),
                                basin = basin_membership)
  result = list(membership = basin_membership,positions = positions)
  return(result)
}

BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

plot_basin = function(basin_df, states_df) {
  states_df$state = sapply(states_df$state, BinToDec)
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
  
  tidy = as_tbl_graph(basin_network, directed = T) 
  
 tidy = tidy %>%
    activate(nodes) %>%
    mutate(energy = states_df$energy[states_df$state %in% name],
           basin = .E()$basin) %>%
    group_by(basin) %>%
    mutate(minimum = ifelse(energy==min(energy), T, F)) %>%
    ungroup() %>%
    activate(edges)
 
 plot = ggraph(tidy)+
   geom_node_voronoi(aes(alpha=energy), fill="navyblue")+
   theme_map()
      
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
  barier_df$maxE = NA
  
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
      
      #compute energy barrier
      E_min1 = states_df$energy[states_df$state==pick[1]]
      E_min2 = states_df$energy[states_df$state==pick[2]]
      diff1 = maxE_last - E_min1
      diff2 = maxE_last - E_min2
      barier_df$barrier[i] = min(diff1, diff2)
      barier_df$maxE[i] = maxE_last
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
    labs(x = "Channel activity",y = "Minima")+
    scale_fill_manual(values=c("#F2F2F2","black"))+
    guides(fill=F)
  
  return(disc_graph/pattern_plot)
}

construct_LMTS = function(barrier_df, states_energy){
  united = barrier_df %>% unite("pair", min1:min2)
  mins = unique(barrier_df$min1)
  find = gtools::combinations(length(mins), 2, mins) %>% as.data.frame() %>% unite("pair", 1:2)
  idx = united$pair %in% find$pair
  
  LMTS_construct = barrier_df[idx,]
  LMTS = matrix(NA, 0, 5)
  saddles = unique(LMTS_construct$saddle) %>% sapply(BinToDec)
  
  #also compute path length 
  path_len = rep(NA,nrow(LMTS_construct))
  
  for (i in 1:nrow(LMTS_construct)){
    row = LMTS_construct[i,]
    barrier = row$maxE
    weight = exp(-row$barrier)
    
    #get reduced network
    states_reduced = states_energy[states_energy$energy<=barrier,]
    A_reduced = get_adjacency(states_reduced$state)
    network_reduced = igraph::graph_from_adjacency_matrix(A_reduced, "undirected")
    transfer = igraph::shortest_paths(network_reduced, from = row$min1, to = row$min2) $vpath[[1]]$name
    
    sub = igraph::induced_subgraph(network_reduced, transfer)
    edgelist = igraph::as_edgelist(sub)
    edgelist = cbind(edgelist,
                     rep(row$min1, nrow(edgelist)), 
                     rep(row$min2, nrow(edgelist)),
                     rep(weight, nrow(edgelist)))
    
    LMTS = rbind(LMTS, edgelist)
    path_len[i] = length(transfer)
  }
  
  LMTS[,1] = sapply(LMTS[,1], BinToDec)
  LMTS[,2] = sapply(LMTS[,2], BinToDec)
  LMTS[,3] = sapply(LMTS[,3], BinToDec)
  LMTS[,4] = sapply(LMTS[,4], BinToDec)
  
  LMTS = as.data.frame(LMTS) %>% unite("Path", 3:4, sep = "-")
  colnames(LMTS)[c(1:2,4)] = c("from","to","weight")
  
  #extract features from LMTS
  n_nodes = length(unique(c(unique(LMTS$from),unique(LMTS$to))))
  n_nodes = n_nodes/length(mins)
  
  #effective path length - path_len - hamming distance
  #stringdist::stringdist(V1,V2, method = "hamming")
  hamming = LMTS_construct %>%
    select(min1, min2) %>%
    mutate(hamming = stringdist::stringdist(min1,min2, method = "hamming")) %>%
    pull()
  
  effective_path = path_len
  effective_path_mean = mean(effective_path)
  
  tidy = as_tbl_graph(LMTS, directed = T) 
  
  mins = sapply(mins,BinToDec)
  tidy = tidy %>% activate(nodes) %>%
    mutate(Node_type = as.factor(case_when(name %in% mins ~ "Minimum",
                                 name %in% saddles ~ "Saddle",
                                 T ~ "1")))
  
  result = list(graph = tidy,
                features = data.frame(nodes = n_nodes,
                                      effective_path = effective_path_mean))
  
  return(result)
}

plot_LMTS = function(tidy_graph){
  plot = tidy_graph %>%
    ggraph()+
    geom_edge_fan(aes(edge_color=Path))+
    geom_node_point(aes(color=Node_type),size=5)+
    scale_color_manual(breaks = c("Minimum", "Saddle"),
                       values = c("grey", "blue", "red"))+
    guides(color=guide_legend(title = "Node type"),
           edge_color = F)+
    theme_few()+
    labs(x = NULL,
         y = NULL)
  return(plot)
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

analyze_energy_landscape = function(file){
  load(file)
  
  fname = strsplit(file, "/")[[1]][4]
  ID = strsplit(fname, "_")[[1]][1]
  diagnosis = strsplit(fname, "_")[[1]][2]
  condition = strsplit(strsplit(fname, "_")[[1]][3],"\\.")[[1]][1]
  
  info = data.frame(ID, diagnosis, condition)
  
  
  #data.frame for storing results
  features = data.frame(minima = NA, E_diff = NA, basin_sd = NA, LMTS_nodes = NA,
                        LMTS_effective = NA, duration = NA)
  
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
  
  #disconnetivity analysis
  barrier_df = find_energy_barrier(minima,states_energy,adj)
  #disconnectivity_graph = plot_disconnectivity(barrier_df)
  
  #LMTS
  if(features$minima>1){
    LMTS = construct_LMTS(barrier_df, states_energy)
    #LMTS_network = plot_LMTS(LMTS$graph)
    features[4:5] = as.numeric(LMTS$features)
  } else {features[4:5] = c(1,0)}
  
  #simulate transitions
  x = simulate_transitions(2e4, adj, states_energy)
  duration = data.frame(state = x) %>%
    left_join(basin_df$membership, by="state") %>%
    group_by(basin) %>%
    summarise(duration = n()) %>%
    mutate(duration = duration/2000)
  features$duration = duration$duration[duration$basin==global_min]
  
  features = cbind(info, features)
  
  #plots = list(basins = basin_plot, disconnectivity_graph = disconnectivity_graph, LMTS = LMTS_network)
  #results = list(features = features,
  #               plots = plots)
  #return(results)
  return(features)
}
