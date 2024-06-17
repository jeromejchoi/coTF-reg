
library(readxl)
library(dplyr)

peak_gene <- read_excel("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Data/Panos/adg3754_Tables_S1_to_S14 3.xlsx", sheet = 4)
peak_gene <- data.frame(peak_gene[-c(1,2),])
colnames(peak_gene) <- c("correlation","gene","peak","zscore","pvalue","pvalue.fdr")
dim(peak_gene) #7291 edges
length(unique(peak_gene$gene)) # 3082 unique TGs

diff_peaks <- read_excel("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Data/Panos/adg3754_Tables_S1_to_S14 3.xlsx", sheet = 3)
diff_peaks <- data.frame(diff_peaks[-c(1,2),])
colnames(diff_peaks) <- c("p_val",	"avg_log2FC",	"pct.1",	"pct.2",	"p_val_adj",	"celltype",	"peak")

diff_peaks_oligo <- diff_peaks[diff_peaks$celltype=="Oligodendrocytes",]
dim(diff_peaks_oligo) # 9167 peaks

oligo_peak_gene <- peak_gene[peak_gene$peak %in% diff_peaks_oligo$peak,]
dim(oligo_peak_gene) # 930 edges
length(unique(oligo_peak_gene$gene)) # 606 TGs

intersect(peak_gene$peak, diff_peaks_oligo$peak)

count_group <- oligo_peak_gene %>% 
  group_by(gene) %>% 
  summarise(n = n()) %>%
  arrange(desc(n))
count_group

count_group1_TGs <- count_group[count_group$n >1,]
dim(count_group1_TGs) # 192 genes have more than one peak

Group1 <- oligo_peak_gene[oligo_peak_gene$gene %in% count_group1_TGs$gene,]
dim(Group1)
head(Group1)

library(tidyr)
Group1_df <- Group1 %>%
  separate(peak, into=c("chr","start","end"), sep="-")
head(Group1_df)

Group1_df$length <- as.numeric(Group1_df$end) - as.numeric(Group1_df$start)
summary(Group1_df$length)
length(unique(Group1_df$gene))
hist(as.numeric(Group1_df$correlation))

write.csv(Group1_df, "/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Group1_df.csv", row.names = F)

Group2 <- oligo_peak_gene[!oligo_peak_gene$gene %in% count_group1_TGs$gene,]
dim(Group2)
head(Group2)
length(unique(Group2$gene))
hist(as.numeric(Group2$correlation))

Group2_df <- Group2

write.csv(Group2_df, "/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Group2.csv",row.names = F)

"OLIG2" %in% count_group_list$gene


# Group 1: Peak_gene & oligo-diff peaks & TGs have more than 1 peak.
# Group 2: Peak_gene & oligo-diff peaks & TGs have only one peak.
# Group 3: Oligo-diff peaks & TGs have more than 1 peak.
# Group 4: Oligo-diff peaks & TGs have only one peak.

# Group 5: All peaks overlap with Bing Ren peaks
# Group 6: Oligo-diff peaks overlap with Bing Ren peak

Group2_df <- read.csv("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Group2.csv")
dim(Group2_df)

library(dplyr);library(tidyr)
Group2_df_check <- Group2_df %>%
  separate(peak, into=c("chr","start","end"), sep ="-")
head(Group2_df_check)

Group2_df_check$end <- as.numeric(Group2_df_check$end)
Group2_df_check$start <- as.numeric(Group2_df_check$start)

summary(Group2_df_check$end - Group2_df_check$start)
hist(Group2_df_check$end - Group2_df_check$start)


Group1_df <- read.csv("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Group1_df.csv")
dim(Group1_df)

library(dplyr);library(tidyr)
Group1_df_check <- Group1_df

Group1_df_check$end <- as.numeric(Group1_df_check$end)
Group1_df_check$start <- as.numeric(Group1_df_check$start)

summary(Group1_df_check$end - Group1_df_check$start)
hist(Group1_df_check$end - Group1_df_check$start)
