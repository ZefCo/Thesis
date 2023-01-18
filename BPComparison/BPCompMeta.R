rm(list=ls())
library(ggplot2)


score_path <- "D:\\Coding\\Thesis\\Data_Files\\BPComp\\Exon2ExonJunctionScores_v3.csv"

score_data <- read.csv(score_path, header = T)

score_data$Classification <- factor(score_data$Classification)
score_data$Ctype <- factor(score_data$Ctype)
# print(levels(score_data$Ctype))
score_data <- score_data[score_data$Ctype == "BRCA",]

global_list <- c("H3_Seq_H5_Intron", "H3_Seq_H5_Exon", "H3_Seq_T5_Seq", "H3_Seq_T3_Intron", "H3_Seq_T3_Exon", "H5_Intron_H5_Exon", "H5_Intron_T5_Seq", "H5_Intron_T3_Intron", "H5_Intron_T3_Exon", "H5_Exon_T5_Seq", "H5_Exon_T3_Intron", "H5_Exon_T3_Exon", "T5_Seq_T3_Intron", "T5_Seq_T3_Exon", "T3_Intron_T3_Exon")
selected_list <- c("H3_Seq_T3_Intron", "H3_Seq_T3_Exon", "H5_Intron_T5_Seq", "H5_Exon_T5_Seq")

global_scores <- data.frame()
for (i in 1:length(global_list)) {
  local_score <- score_data[[global_list[[i]]]]
  #local_norm <- local_score / sum(local_score)
  local_factor <- factor(rep(global_list[[i]], length(local_score)))
  
  local_frame <- data.frame(scores = local_score, location = local_factor) #, norm = local_norm)
  
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

#global_scores <- global_scores / global_scores

hcolors <- c("H3_Seq_T3_Intron" = "red", "H3_Seq_T3_Exon" = "green", "H5_Intron_T5_Seq" = "blue", "H5_Exon_T5_Seq" = "purple")

ggglobal <- ggplot(global_scores, aes(x=location, y=scores, fill=location))  + 
  geom_boxplot() + theme_grey() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + labs(title = "All Comparisons")
print(ggglobal)
ggselected <- ggplot(selected_scores, aes(x=location, y=scores, fill=location)) + 
  geom_boxplot() + theme_grey() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + labs(title = "Selected Comparisons")
print(ggselected)
gghisto <- ggplot(selected_scores, aes(x = scores, fill = location)) + geom_histogram(position = "dodge", bins = 50) #+ scale_fill_manual(values = hcolors)
#ggglobal <- ggglobal + ggplot(randomGlobal, aes(x = x, y = y))
print(gghisto)

distanceplot <- ggplot(score_data, aes(x = Classification, color = Classification, fill = Classification)) + geom_bar() + labs(title = "Global Distance Classification")
print(distanceplot)

# breast_data <- score_data[score_data$Ctype == "BRCA"]

# ggglobal
# ggselected
# gghisto
