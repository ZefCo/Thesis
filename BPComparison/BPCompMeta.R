rm(list=ls())
library(ggplot2)

score_path <- "D:\\Coding\\Thesis\\Data_Files\\BPComp\\Fusion_Scores_min100.csv"

score_data <- read.csv(score_path, header = T)

score_data$Classification <- factor(score_data$Classification)

global_list <- c("H3_Seq_H5_Intron", "H3_Seq_H5_Exon", "H3_Seq_T5_Seq", "H3_Seq_T3_Intron", "H3_Seq_T3_Exon", "H5_Intron_H5_Exon", "H5_Intron_T5_Seq", "H5_Intron_T3_Intron", "H5_Intron_T3_Exon", "H5_Exon_T5_Seq", "H5_Exon_T3_Intron", "H5_Exon_T3_Exon", "T5_Seq_T3_Intron", "T5_Seq_T3_Exon", "T3_Intron_T3_Exon")
selected_list <- c("H3_Seq_T3_Intron", "H3_Seq_T3_Exon", "H5_Intron_H5_Exon", "H5_Exon_T5_Seq")

global_scores <- data.frame()
for (i in 1:length(global_list)) {
  local_score <- score_data[[global_list[[i]]]]
  local_factor <- factor(rep(global_list[[i]], length(local_score)))
  
  local_frame <- data.frame(scores = local_score, location = local_factor)
  
  global_scores <- rbind(global_scores, local_frame)

  
}


selected_scores <- data.frame()
for (i in 1:length(selected_list)) {
  selected_score <- score_data[[selected_list[[i]]]]
  selected_factor <- factor(rep(selected_list[[i]], length(selected_score)))
  
  selected_frame <- data.frame(scores = selected_score, location = selected_factor)
  
  selected_scores <- rbind(selected_scores, selected_frame)
}

randomGlobal <- data.frame(x = global_list, y = rep(0.75, length(global_list)))
randomSelected <- data.frame(x = selected_list, y = rep(0.75, length(selected_list)))


ggglobal <- ggplot(global_scores, aes(x=location, y=scores, fill=location))  + 
  geom_boxplot() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + theme_grey() 
ggselected <- ggplot(selected_scores, aes(x=location, y=scores, fill=location)) + 
  geom_boxplot() + theme(axis.title.x = element_blank()) + theme_grey()

#ggglobal <- ggglobal + ggplot(randomGlobal, aes(x = x, y = y))

ggglobal
ggselected
