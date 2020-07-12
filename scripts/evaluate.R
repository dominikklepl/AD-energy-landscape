library(tidyverse)
library(ggthemes)
library(patchwork)

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
  df_EC = matrix(NA, length(files), 7)
  colnames(df_EC) = c("selection","wave","condition","ß", "error", "t", "p-value")
  df_EC = as.data.frame(df_EC)
  df_EO = matrix(NA, length(files), 7)
  colnames(df_EO) = c("selection","wave","condition","ß", "error", "t", "p-value")
  df_EO = as.data.frame(df_EO)
  
  df_EC$condition = "EC"
  df_EO$condition = "EO"
  
  for(i in 1:length(files)){
    file = files[i]
    data = read.csv(file) %>%
      select("diagnosis","condition",all_of(feature)) %>%
      rename("value" = feature)
    
    split = strsplit(strsplit(strsplit(file, "/")[[1]][3], "\\.")[[1]][1],"_")[[1]]
    
    df_EC$selection[i] = split[1]
    df_EC$wave[i] = split[2]
    df_EO$selection[i] = split[1]
    df_EO$wave[i] = split[2]
    
    EC = data %>% filter(condition=="EC")
    EO = data %>% filter(condition=="EO")
    
    #run t.tests
    EC_lm = as.data.frame(summary(lm(value ~ diagnosis, EC))$coefficients)
    EO_lm = as.data.frame(summary(lm(value ~ diagnosis, EO))$coefficients)
    
    df_EC$ß[i] = EC_lm$Estimate[2] %>% round(3)
    df_EC$error[i] = EC_lm$`Std. Error`[2] %>% round(3)
    df_EC$t[i] = EC_lm$`t value`[2]  %>% round(3)
    df_EC$`p-value`[i] = EC_lm$`Pr(>|t|)`[2]  %>% round(3)
    
    df_EO$ß[i] = EO_lm$Estimate[2] %>% round(3)
    df_EO$error[i] = EO_lm$`Std. Error`[2] %>% round(3)
    df_EO$t[i] = EO_lm$`t value`[2]  %>% round(3)
    df_EO$`p-value`[i] = EO_lm$`Pr(>|t|)`[2]  %>% round(3)
  }
  
  result = rbind(df_EC, df_EO) %>% arrange(selection, wave)
  return(result)
}

plot_difference = function(results, feature){
  data = signif_results[[feature]]
  all_diffs = data.frame()
  for (j in 1:nrow(data)){
    load_this = paste0("results/ELA/", data$selection[j],"_",data$wave[j], ".csv")
    diff = read.csv(load_this) %>%
      group_by(condition, diagnosis) %>%
      select("condition", "diagnosis", feature) %>%
      rename(value = feature) %>%
      summarise(mean = mean(value),
                sd = sd(value)/1.98) %>%
      ungroup() %>%
      filter(condition==data$condition[j])
    selection = ifelse(data$selection[j]=="R", "ML", data$selection[j])
    diff$label = paste0(selection,"-",data$wave[j])
    all_diffs = rbind(all_diffs, diff)
  }
  all_diffs$condition = ifelse(all_diffs$condition=="EC", "Eyes closed", "Eyes open")
  all_diffs$diagnosis = ifelse(all_diffs$diagnosis=="AD", "Alzheimer's", "Healthy")
  
  plot = ggplot(all_diffs, aes(x = label, y = mean, fill = diagnosis))+
    geom_bar(stat = "identity", position="dodge2")+
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), position = position_dodge2(width = 0.5, padding = 0.5))+
    theme_few()+
    scale_fill_tableau(name = "Diagnosis")+
    labs(x = NULL,
         y = feature)+
    facet_grid(~condition, scales = "free_x", space = "free_x")
  
  return(plot)
}
