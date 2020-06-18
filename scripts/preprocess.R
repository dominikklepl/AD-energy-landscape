library(ggplot2)
library(ggthemes)
library(wesanderson)
library(eegkit)
library(fastICA)
library(tictoc)

get_info = function(filename) {
  slash_split = strsplit(filename,"/")
  diagnosis = strsplit(slash_split[[1]][3], " ")[[1]][1]
  condition = strsplit(slash_split[[1]][3], " ")[[1]][2]
  ID = strsplit(slash_split[[1]][4], " ")[[1]][1]
  return(cbind(ID,diagnosis,condition))
}

load_data = function(path) {
  data = read.delim(path, header=T)
  
  if (nrow(data)%%2!=0){
    data = data[-nrow(data),]
  }
  
  #attach participant info as attribute
  attr(data,"info") = get_info(path)
  
  #fix channel names
  colnames(data) = c("Time","F8_F4","F7_F3","F4_C4","F3_C3","F4_FZ","F3_FZ","FZ_CZ","T4_C4","T3_C3","C4_CZ","C3_CZ","CZ_PZ","C4_P4","C3_P3","T4_T6","T3_T5","P4_PZ","P3_PZ","T6_O2","T5_O1","P4_O2","P3_O1","O1_O2")
  
  return(data)
}

#filters
filter_notch = function(x, FS = 2000, lower=49, upper=50) {
  len = length(x)
  
  transform = fft(x)/len
  N = FS/len
  freq = seq(0,N*(len/2-1),N)
  
  trans = data.frame(f=transform,freq=c(freq,rev(freq)))
  trans$f[trans$freq>lower&trans$freq<upper] = 0
  
  cleaned = Re(fft(trans$f, inverse = T))
  return(cleaned)
}

filter_bandpass = function(x, FS = 2000, lower=NULL, upper=NULL) {
  len = length(x)
  
  transform = fft(x)/len
  N = FS/len
  freq = seq(0,N*(len/2-1),N)
  
  trans = data.frame(f=transform,freq=c(freq,rev(freq)))
  
  if(is.null(upper)==F){trans$f[trans$freq>upper] = 0}
  if(is.null(lower)==F){trans$f[trans$freq<lower] = 0}
  
  cleaned = Re(fft(trans$f, inverse = T))
  return(cleaned)
}

amplitude_envelope = function(x){
  s  =  hht::HilbertTransform(x)
  return(hht::HilbertEnvelope(s))
}

#Plotting
geom_head = function(){
  #head shape
  rad = 12.5
  head = data.frame(xx = rad * cos(seq(0, 2 * pi, length.out = 360)),
                    yy = rad * sin(seq(0, 2 * pi, length.out = 360)))
  
  #nose
  nose1 = data.frame(xx=c(head$xx[81], 0),
                     yy=c(head$yy[81], rad * 1.175))
  nose2 = data.frame(xx=c(-head$xx[81], 0),
                     yy=c(head$yy[81], rad * 1.175))
  
  #ears
  ear1 = data.frame(xx = (0.5 * cos(seq(0, 2 * pi, l = 360)))-13,
                    yy = (2.5 * sin(seq(0, 2 * pi, l = 360))))
  ear2 = data.frame(xx = (0.5 * cos(seq(0, 2 * pi, l = 360)))+13,
                    yy = (2.5 * sin(seq(0, 2 * pi, l = 360))))
  pal = wesanderson::wes_palette("Zissou1", 100, type = "continuous")
  plot = ggplot(head,aes(xx,yy))+
    geom_path()+
    coord_cartesian(xlim=(c(-13, 13)),
                    ylim = c(-12.5, 14))+
    geom_path(data=nose1, aes(xx,yy))+
    geom_path(data=nose2, aes(xx,yy))+
    geom_path(data=ear1, aes(xx,yy))+
    geom_path(data=ear2, aes(xx,yy))+
    theme_void()+
    labs(x=NULL,
         y=NULL)+
    scale_color_gradientn(colours = pal)
  return(plot)
}

get_channel_loc = function(){
  #channel positions
  data(eegcoord)
  
  #get first part of the bipolar channels
  names = c("F8_F4","F7_F3","F4_C4","F3_C3","F4_FZ","F3_FZ","FZ_CZ","T4_C4","T3_C3","C4_CZ","C3_CZ","CZ_PZ","C4_P4","C3_P3","T4_T6","T3_T5","P4_PZ","P3_PZ","T6_O2","T5_O1","P4_O2","P3_O1","O1_O2")
  channel_positions = matrix(0,length(names),2)
  colnames(channel_positions) = c("x","y")
  rownames(channel_positions) = names
  
  for (i in 1:length(names)){
    pair_name = strsplit(names[i],"_")[[1]]
    pair =  eegcoord[pair_name,4:5]
    
    #compute middle
    middle = c(min(pair$xproj)+abs(((pair$xproj[1]-pair$xproj[2])/2)),
               min(pair$yproj)+abs(((pair$yproj[1]-pair$yproj[2])/2)))
    channel_positions[i,] = middle
  }
  
  #manually fill missing channels
  #T4-C4
  channel_positions["T4_C4",] = c(eegcoord["C4",4]+abs((eegcoord["C4",4]-eegcoord["T8",4])/2),
                                  eegcoord["T8",5]+abs((eegcoord["C4",5]-eegcoord["T8",5])/2))
  
  #T3-C3
  channel_positions["T3_C3",] = c(eegcoord["T7",4]+abs((eegcoord["C3",4]-eegcoord["T7",4])/2),
                                  eegcoord["T7",5]+abs((eegcoord["C3",5]-eegcoord["T7",5])/2))
  
  #T4-T6
  channel_positions["T4_T6",] = c(eegcoord["P8",4]+abs((eegcoord["P8",4]-eegcoord["T8",4])/2),
                                  eegcoord["P8",5]+abs((eegcoord["P8",5]-eegcoord["T8",5])/2))
  
  #T3-T5
  channel_positions["T3_T5",] = c(eegcoord["T7",4]+abs((eegcoord["P7",4]-eegcoord["T7",4])/2),
                                  eegcoord["P7",5]+abs((eegcoord["P7",5]-eegcoord["T7",5])/2))
  
  #T6-O2
  channel_positions["T6_O2",] = c(eegcoord["O2",4]+abs((eegcoord["P8",4]-eegcoord["O2",4])/2),
                                  eegcoord["O2",5]+abs((eegcoord["P8",5]-eegcoord["O2",5])/2))
  
  #T5-O1
  channel_positions["T5_O1",] = c(eegcoord["P7",4]+abs((eegcoord["P7",4]-eegcoord["O1",4])/2),
                                  eegcoord["O1",5]+abs((eegcoord["P7",5]-eegcoord["O1",5])/2))
  return(as.data.frame(channel_positions))
}

plot_ica = function(data,column){
  geom_head()+
    geom_point(data=data, aes_string("x","y",color=column),size=5)+
    theme(legend.title=element_blank())
}

run_full_ica = function(data,ncomp = 17, locations){
  ica = fastICA(data[,-1], ncomp)
  
  #compute variance explained
  vars = rowSums(ica$A^2)
  vars = (vars * nrow(ica$X))/sum(ica$X^2)
  
  mixing_ordered = ica$A[order(vars,decreasing = T),]
  
  icadata = cbind(locations,t(mixing_ordered))
  colnames(icadata)[3:(ncomp+2)] = paste0("Comp",1:ncomp)
  
  ica_plots = lapply(colnames(icadata)[3:(ncomp+2)],plot_ica,data=icadata)
  
  plot_grid = ggpubr::ggarrange(plotlist = ica_plots, labels = colnames(icadata)[3:(ncomp+2)], ncol = 2, nrow = round(ncomp/2,0)+1)
  
  return(plot_grid)
}

run_half_ica = function(data,ncomp = 8, locations, group){
  data = data[,group]
  locations = locations[group,]
  
  ica = fastICA(data, ncomp)
  
  #compute variance explained
  vars = rowSums(ica$A^2)
  vars = (vars * nrow(ica$X))/sum(ica$X^2)
  
  mixing_ordered = ica$A[order(vars,decreasing = T),]
  
  icadata = cbind(locations,t(mixing_ordered))
  colnames(icadata)[3:(ncomp+2)] = paste0("Comp",1:ncomp)
  
  ica_plots = lapply(colnames(icadata)[3:(ncomp+2)],plot_ica,data=icadata)
  
  plot_grid = ggpubr::ggarrange(plotlist = ica_plots, labels = colnames(icadata)[3:(ncomp+2)], ncol = 2, nrow = round(ncomp/2,0)+1)
  
  return(plot_grid)
}


preprocess = function(path){
  data = load_data(path)
  
  #filter - notch and highpass (100 Hz)
  data[,-1] = apply(data[,-1], 2, filter_notch)
  data[,-1] = apply(data[,-1], 2, filter_bandpass, upper=50)
  
  data[,-1] = apply(data[,-1],2, amplitude_envelope)
  
  return(data)
}