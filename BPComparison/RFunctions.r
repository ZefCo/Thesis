library(gdata)


import_data <- function(datapath) {
    return_data <- read.csv(datapath, header = T)

    return()
}



bar_and_whiskers <-function(inputData) {
    header_list <- c("H3_Seq_H5_Intron", "H3_Seq_H5_Exon", "H3_Seq_T5_Seq", "H3_Seq_T3_Intron", "H3_Seq_T3_Exon", "H5_Intron_H5_Exon", "H5_Intron_T5_Seq", "H5_Intron_T3_Intron", "H5_Intron_T3_Exon", "H5_Exon_T5_Seq", "H5_Exon_T3_Intron", "H5_Exon_T3_Exon", "T5_Seq_T3_Intron", "T5_Seq_T3_Exon", "T3_Intron_T3_Exon")

    # ggplot(inputData, aes(x=))
}
