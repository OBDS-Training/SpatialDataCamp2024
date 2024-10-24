library(tidyverse)

file <- "/ceph/project/simmons_hts/aantanav/_r_projects/CCB_SUMMER_SCHOOL_2024/0_DATASETS/XENIUM_COLON/transcripts.csv.gz"
transcripts <- read_csv(file)
transcripts$cell_id2 <- as.numeric(as.factor(transcripts$cell_id))
transcripts$cell_id2 [ transcripts$cell_id == "UNASSIGNED"] <- 0
transcripts$cell_id <- transcripts$cell_id2
write.csv(transcripts, file="/ceph/project/simmons_hts/aantanav/_r_projects/CCB_SUMMER_SCHOOL_2024/2_BAYSOR_REFORMAT/COLON_XENIUM/detected_transcripts_for_baysor.csv", row.names = FALSE)
  




  