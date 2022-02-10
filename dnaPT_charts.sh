library(tidyverse)
library(ggrepel)
library(cowplot)
rca<-read.table("simulans_dnaPipeTE_0.15_1_t20/reads_per_component_and_annotation", fill = T, na.strings=c("","NA"))
rca<-separate(rca, V6, c("subclass", "superfamily"), sep = "/",fill = "right") 
rca$superfamily[is.na(rca$superfamily)]<-"Unknown"
rca$subclass[is.na(rca$subclass)]<-"Unknown"
# print("consolidate SF column...")
rca$superfamily<-paste(rca$subclass, rca$superfamily, sep = "/")
# print("renames NA in SF columns...")
rca$superfamily[rca$superfamily == "Unknown/Unknown"]<-"Unknown"
names(rca)<-c("reads", "bp", "contig", "q_length", "target_name", "t_subc", "t_supf", "hit_cov")
reads.c<-as.numeric(system("grep -c '>' melano_dnaPipeTE_0.15_1_t20/renamed.blasting_reads.fasta", intern = T))
reads.b<-as.numeric(system("tail -n 2 melano_dnaPipeTE_0.15_1_t20/Counts.txt | head -n 1 | cut -f 2", intern = T))
threshold=0.001
rca<-rca[order(rca$bp, decreasing = T),]
leftover<-c(
  sum(rca[rca$bp < threshold*reads.b/100,]$reads),
  sum(rca[rca$bp < threshold*reads.b/100,]$bp),
  rep(NA, 3),
  paste("repeats_under_", threshold, "%", sep = ""),
  paste("repeats_under_", threshold, "%", sep = ""),
  NA)
rca_t<-rca[rca$bp >= threshold*reads.b/100,]
rca_t<-rbind(rca_t, leftover)
rca_t$perc_g<-as.numeric(rca_t$bp)/reads.b*100
rca_t$x2<-cumsum(as.numeric(rca_t$perc_g))
rca_t$x1<-as.numeric(rca_t$x2)-as.numeric(rca_t$perc_g)

# pick the colors
cols<-read.table("../dnaPT_utils/colors.land", sep = "\t")
cols[length(cols$V1)+1,]<-c(paste("repeats_under_", threshold, "%", sep = ""), "grey10")
col.bars<-rep("", length((levels(as.factor(rca_t$t_subc)))))
for(i in 1:length(levels(as.factor(rca_t$t_subc)))){
  col.bars[i]<-cols$V2[grep(pattern = paste("^", levels(as.factor(rca_t$t_subc))[i], "$", sep = ""), x = cols$V1)]
}

barplot<-ggplot(rca_t)+
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = perc_g, fill = t_subc))+
  scale_fill_manual(values = col.bars, name = "repeat type")+
  xlab("cumulative genome %")+
  ylab("genome %")+
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank())

############ PIE CHARTS #############################

rca_pie<-as.data.frame(tapply(rca_t$perc_g, rca_t$t_subc, sum))
rca_pie$subc<-rownames(rca_pie)
names(rca_pie)<-c("perc_g", "subc")
rca_pie$perc_r<-rca_pie$perc_g/sum(rca_pie$perc_g)*100


###### plot in % TE ################################
df2 <- rca_pie %>% 
  mutate(csum = rev(cumsum(rev(perc_r))), 
         pos = perc_r/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc_r/2, pos))
df2$subc<-factor(df2$subc, levels = df2$subc, labels = paste(paste0(round(df2$perc_r,2), "% ",df2$subc)))
pie_rel<-ggplot(df2, aes(x = "", y = perc_r, fill = subc))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y")+
  scale_fill_manual(values = col.bars, name = "repeat type (relative %)")+
  # geom_label_repel(data = df2,
  #                  aes(y = pos, label = paste0(round(perc_r,2), "%")),
  #                  size = 4.5, nudge_x = 1.5, show.legend = FALSE, fill = "white")+ #, col = col.bars)+
  theme_void()
  
###### plot in % genome ############################
nr<-100-sum(rca_pie$perc_g)
rca_pie_g<-rbind(rca_pie[,1:2], c(nr, "non-repetitive"))
col.bars_g<-c(col.bars, "grey90")
rca_pie_g$perc_g<-as.numeric(rca_pie_g$perc_g)
dfg <- rca_pie_g %>% 
  mutate(csum = rev(cumsum(rev(perc_g))), 
         pos = perc_g/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc_g/2, pos))
dfg$subc<-factor(dfg$subc, levels = dfg$subc, labels = paste(paste0(round(dfg$perc_g,2), "% ",dfg$subc)))
pie_g<-ggplot(dfg, aes(x = "", y = perc_g, fill = subc))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y")+
  scale_fill_manual(values = col.bars_g, name = "repeat type (genome %)")+
  #geom_label_repel(data = dfg,
  #                 aes(y = pos, label = paste0(round(perc_g,2), "%")),
  #                 size = 4.5, nudge_x = 3,min.segment.length = 15, show.legend = FALSE, fill = "white")+ #, col = col.bars)+
  theme_void()


top<-plot_grid(pie_rel, pie_g, ncol = 2)
plot_grid(top, barplot, ncol = 1)
