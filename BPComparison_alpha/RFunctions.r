library(NGCHM)
library(psych)
library(openxlsx)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gdata)

rm(list=ls())



# Very boring: just imports the data.
importData <- function(file_path) {

  if (endsWith(file_path, ".csv")) {
    imported_data <- read.csv(file_path, header=T, sep=',')
  } else if (endsWith(file_path, ".xlsx")) {
    imported_data <- read.xlsx(file_path)

  }

  return(imported_data)
  
}


# 
findMineLength <- function(input_data, target_col = "SeqLen", min_length = 20) {

  return_data <- input_data[input_data[target_col] >= min_length, ]

  return(return_data)
}

