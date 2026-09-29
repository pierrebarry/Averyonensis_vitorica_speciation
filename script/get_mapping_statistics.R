e =list.files("/shared/projects/origamis/wgs/mapping",full.names = T)

file = c()
for (i in e){
  if (grepl("flagstat",i)){
    file = c(file,i)
  }
}

library(tidyverse)

mapping_stats = data.frame(SAMPLE = as.character(),
                           TOTAL_READS = as.numeric(),
                           SECONDARY = as.numeric(),
                           SUPPLEMENTARY= as.numeric(),
                           DUPLICATES = as.numeric(),
                           PERCENTAGE_DUPLICATES = as.numeric(),
                           MAPPED = as.numeric(),
                           PERCENTAGE_MAPPED = as.numeric(),
                           PAIRED_IN_SEQUENCING = as.numeric(),
                           PROPERLY_PAIRED = as.numeric(),
                           PERCENTAGE_PROPERLY_PAIRED = as.numeric(),
                           SINGLETON = as.numeric(),
                           PERCENTAGE_SINGLETON = as.numeric(),
                           MATE_DIFFERENT_CHROMOSOME = as.numeric()
)

for (j in file){
  tmp = read.delim(j,header=F) 
  
  
  tmp_data = data.frame(SAMPLE = as.character(strsplit(i,".txt")[[1]][1]),
                        TOTAL_READS = as.numeric(strsplit(tmp[1,]," ")[[1]][1]),
                        SECONDARY = as.numeric(strsplit(tmp[2,]," ")[[1]][1]),
                        SUPPLEMENTARY= as.numeric(strsplit(tmp[3,]," ")[[1]][1]),
                        DUPLICATES = as.numeric(strsplit(tmp[4,]," ")[[1]][1]),
                        PERCENTAGE_DUPLICATES = round(as.numeric(strsplit(tmp[4,]," ")[[1]][1])*100/as.numeric(strsplit(tmp[1,]," ")[[1]][1]),2),
                        MAPPED = as.numeric(strsplit(tmp[5,]," ")[[1]][1]),
                        PERCENTAGE_MAPPED = as.numeric(str_extract(strsplit(strsplit(tmp[5,]," ")[[1]][5],":")[[1]],"\\d+\\.*\\d*")),
                        PAIRED_IN_SEQUENCING = as.numeric(strsplit(tmp[6,]," ")[[1]][1]),
                        PROPERLY_PAIRED = as.numeric(strsplit(tmp[9,]," ")[[1]][1]),
                        PERCENTAGE_PROPERLY_PAIRED = as.numeric(str_extract(strsplit(strsplit(tmp[9,]," ")[[1]][6],":")[[1]],"\\d+\\.*\\d*")),
                        SINGLETON = as.numeric(strsplit(tmp[11,]," ")[[1]][1]),
                        PERCENTAGE_SINGLETON = as.numeric(str_extract(strsplit(strsplit(tmp[11,]," ")[[1]][5],":")[[1]],"\\d+\\.*\\d*")),
                        MATE_DIFFERENT_CHROMOSOME = as.numeric(strsplit(tmp[12,]," ")[[1]][1]))
  
  mapping_stats = rbind(mapping_stats,
                        tmp_data)
  
}

write.csv(dd,file="/shared/projects/origamis/wgs/mapping/mapping.csv",row.names = F,quote=F)

