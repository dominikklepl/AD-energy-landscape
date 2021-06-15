library(tidyverse)
library(ggthemes)
library(patchwork)
library(lmerTest)

merge_results = function(pMEM_files, selection){
  result = data.frame()
  for (i in 1:length(pMEM_files)){
    if (grepl(selection, pMEM_files[i])){
      one_res = read.csv(pMEM_files[i])
      one_res$band = strsplit(strsplit(pMEM_files[i],"_")[[1]][3], "\\.")[[1]][1]
      result = rbind(result, one_res)
    }
  }
  return(result)
}

evaluate_pMEM = function(result_df, selection){
  r_test = t.test(r ~ diagnosis, result_df)
  R2_test = t.test(R2 ~ diagnosis, result_df)
  
  r_test = paste0("t(",round(r_test$parameter,2),") = ", round(r_test$statistic,2), ", p-value = ", round(r_test$p.value,2))
  R2_test = paste0("t(",round(R2_test$parameter,2),") = ", round(R2_test$statistic,2), ", p-value = ", round(R2_test$p.value,2))
  
  test_result = as.data.frame(rbind(r_test, R2_test))
  write_to = paste0("results/pMEM_evaluation/", selection, ".csv")
  readr::write_csv(test_result, write_to)
  
  summary_df = result_df %>%
    group_by(diagnosis, condition) %>%
    summarise(mean_r = mean(r),
              sd_r = sd(r),
              mean_R2 = mean(R2),
              sd_R2 = sd(R2))
  
  r_plot = ggplot(summary_df, aes(x=diagnosis, y = mean_r,fill=diagnosis))+
    geom_bar(stat="identity")+
    geom_errorbar(aes(ymin = mean_r - sd_r, ymax = mean_r + sd_r), width = 0.4)+
    theme_few()+
    scale_fill_tableau()+
    facet_wrap(~condition)+
    labs(x = NULL,
         y = bquote(r[D]))
  
  R2_plot = ggplot(summary_df, aes(x=diagnosis, y = mean_R2,fill=diagnosis))+
    geom_bar(stat="identity")+
    geom_errorbar(aes(ymin = mean_R2 - sd_R2, ymax = mean_R2 + sd_R2), width = 0.4)+
    theme_few()+
    scale_fill_tableau()+
    facet_wrap(~condition)+
    guides(fill=F)+
    labs(x = NULL,
         y = bquote(R^2))
  
  pMEM_plot = r_plot + R2_plot + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
  
  pMEM_plot
  save_to = paste0("results/pMEM_evaluation/", selection, ".png")
  ggsave(save_to, width = 7, height = 5)
}

test_effect = function(files, feature){
  df_EC = matrix(NA, length(files), 9)
  colnames(df_EC) = c("selection","wave","condition","intercept","error_i","ß","error_ß", "t", "p-value")
  df_EC = as.data.frame(df_EC)
  df_EO = matrix(NA, length(files), 9)
  colnames(df_EO) = c("selection","wave","condition","intercept","error_i","ß","error_ß", "t", "p-value")
  df_EO = as.data.frame(df_EO)
  
  df_EC$condition = "EC"
  df_EO$condition = "EO"
  
  for(i in 1:length(files)){
    file = files[i]
    data = read.csv(file) %>%
      select("diagnosis","condition","ID",all_of(feature)) %>%
      rename("value" = feature)
    
    split = strsplit(strsplit(strsplit(file, "/")[[1]][3], "\\.")[[1]][1],"_")[[1]]
    
    df_EC$selection[i] = split[1]
    df_EC$wave[i] = split[2]
    df_EO$selection[i] = split[1]
    df_EO$wave[i] = split[2]
    
    EC = data %>% filter(condition=="EC")
    EO = data %>% filter(condition=="EO")
    
    #run t.tests - in fact mixed effects glms
    EC_lm = lmer(value ~ diagnosis + (1|ID), EC)
    EO_lm = lmer(value ~ diagnosis + (1|ID), EO)
    
    EC_summary = as.data.frame(summary(EC_lm)$coefficients)
    EO_summary = as.data.frame(summary(EO_lm)$coefficients)
    
    df_EC$intercept[i] = EC_summary$Estimate[1]
    df_EC$error_i[i] = EC_summary$`Std. Error`[1]
    df_EC$ß[i] = EC_summary$Estimate[2]
    df_EC$error_ß[i] = EC_summary$`Std. Error`[2]
    df_EC$t[i] = EC_summary$`t value`[2]
    df_EC$`p-value`[i] = EC_summary$`Pr(>|t|)`[2]
    
    df_EO$intercept[i] = EO_summary$Estimate[1]
    df_EO$error_i[i] = EO_summary$`Std. Error`[1]
    df_EO$ß[i] = EO_summary$Estimate[2]
    df_EO$error_ß[i] = EO_summary$`Std. Error`[2]
    df_EO$t[i] = EO_summary$`t value`[2]
    df_EO$`p-value`[i] = EO_summary$`Pr(>|t|)`[2]
  }
  
  result = rbind(df_EC, df_EO) %>% arrange(selection, wave)
  return(result)
}

plot_difference_data = function(results, feature){
  data = results[[feature]]
  
  result = foreach(i=1:nrow(data), .combine = "rbind") %do% {
    read.csv(paste0("results/ELA/", data$selection[i], "_", data$wave[i], ".csv")) %>%
      filter(condition == data$condition[i]) %>%
    mutate(label = data$wave[i]) %>%
      select(diagnosis, condition, label, feature)
  }
  colnames(result)[4] = "value"
  result$diagnosis = ifelse(result$diagnosis=="AD", "Alzheimer's", "Healthy")
  result$condition = ifelse(result$condition=="EC", "Eyes Closed", "Eyes Open")
  
  plot = ggplot(result, aes(x = label, y = value, fill = diagnosis))+
    stat_summary(geom = "bar", fun = mean, position = "dodge2") +
    stat_summary(geom = "errorbar", fun.data = mean_cl_normal, position = position_dodge2(width = 0.5, padding = 0.5))+
    theme_few()+
    scale_fill_manual(values = c("#BA0C2F","#00629B"), name = "Diagnosis")+
    labs(x = NULL,
         y = feature)+
    facet_wrap(~condition, scales = "free_x")+
    theme(axis.text = element_text(size=8),
          axis.title.y = element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=9),
          strip.text = element_text(size=8))
  
  return(plot)
}
