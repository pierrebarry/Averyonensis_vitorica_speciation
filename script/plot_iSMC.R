library(ggplot2)
library(stringi)
ll = list.files("/shared/projects/origamis/wgs/output/iSMC/", full.names = T)

data = data.frame(SAMPLE = as.character(),
                  CHROM = as.character(),
                  WINDOW = as.character(),
                  START = as.numeric(),
                  RHO = as.numeric())

for (i in ll){
  
  if (stri_sub(i, -4)==".csv"){
    tmp = read.csv(i)
    
    data = rbind(data,
                 data.frame(SAMPLE = c(strsplit(strsplit(i,"//")[[1]][2],"_")[[1]][2]),
                            CHROM = c(strsplit(strsplit(i,"//")[[1]][2],"_")[[1]][1]),
                            WINDOW = c(strsplit(strsplit(strsplit(i,"//")[[1]][2],"_")[[1]][4],".csv")[[1]][1]),
                            START = tmp$POS,
                            RHO = tmp$RHO
                 ))
  }
  
  
}

data$CHROM <- factor(data$CHROM)



p<-data[data$WINDOW=="1000000kb",] %>%
  mutate(chrom_color_group = case_when(as.numeric(CHROM) %% 2 != 0 ~ "even",
                                       CHROM == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(CHROM,labels=seq(1,length(levels(factor(CHROM)))))) %>%
  ggplot(aes(x=START,y=RHO,group=SAMPLE)) +
  geom_step(aes(col=chrom_color_group),alpha=0.5,lwd=1) +
  theme_bw() +
  facet_grid(WINDOW~chromosome,scales="free_x",switch="x",space = "free_x") +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_color_manual(values = c("grey50", "black"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression(rho)) +
  xlab("Chromosome")

pdf("/shared/projects/origamis/GBS/figures/rho_1Mb.pdf",width=20,height=2.5)
print(p)
dev.off()

p<-data[data$WINDOW=="5000000kb",] %>%
  mutate(chrom_color_group = case_when(as.numeric(CHROM) %% 2 != 0 ~ "even",
                                       CHROM == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(CHROM,labels=seq(1,length(levels(factor(CHROM)))))) %>%
  ggplot(aes(x=START,y=RHO,group=SAMPLE)) +
  geom_step(aes(col=chrom_color_group),alpha=0.5,lwd=1) +
  theme_bw() +
  facet_grid(WINDOW~chromosome,scales="free_x",switch="x",space = "free_x") +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_color_manual(values = c("grey50", "black"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression(rho)) +
  xlab("Chromosome")

pdf("/shared/projects/origamis/GBS/figures/rho_5Mb.pdf",width=20,height=2.5)
print(p)
dev.off()

p<-data[data$WINDOW=="10000000kb",] %>%
  mutate(chrom_color_group = case_when(as.numeric(CHROM) %% 2 != 0 ~ "even",
                                       CHROM == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(CHROM,labels=seq(1,length(levels(factor(CHROM)))))) %>%
  ggplot(aes(x=START,y=RHO,group=SAMPLE)) +
  geom_step(aes(col=chrom_color_group),alpha=0.5,lwd=1) +
  theme_bw() +
  facet_grid(WINDOW~chromosome,scales="free_x",switch="x",space = "free_x") +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_color_manual(values = c("grey50", "black"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab(expression(rho)) +
  xlab("Chromosome")

pdf("/shared/projects/origamis/GBS/figures/rho_10Mb.pdf",width=20,height=2.5)
print(p)
dev.off()
