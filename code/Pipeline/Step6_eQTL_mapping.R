library(dplyr);library(tidyr);library(IRanges);library(gtools)

loadRData <- function(file_name){
  load(file_name)
  get(ls()[ls() != "file_name"])
}

#  Percentage of eQTL SNPs that are located outside or inside of their paired gene region (enhmoter)
# When SNPs were associated with expression of multiple genes, they were counted multiple times.

# SNP position
SNP_posi <- read.table("snp_pos.txt.gz", 
                       header = T)
dim(SNP_posi)
head(SNP_posi)

SNP_posi <- SNP_posi %>%
  separate(SNP_id_hg38, into=c("chr_SNP","position"),sep = ":")

# oligo-diff data
diff_peaks <- read_excel("adg3754_Tables_S1_to_S14.xlsx", sheet = 3)
diff_peaks <- data.frame(diff_peaks[-c(1,2),])
colnames(diff_peaks) <- c("p_val",	"avg_log2FC",	"pct.1",	"pct.2",	"p_val_adj",	"celltype",	"peak")

diff_peaks_oligo <- diff_peaks[diff_peaks$celltype=="Oligodendrocytes",]
dim(diff_peaks_oligo)
length(unique(diff_peaks_oligo$peak))

diff_peaks_oligo <- diff_peaks_oligo %>%
  filter(!duplicated(peak)) %>%
  separate(peak, into=c("chr","start","end"), sep = "-")
dim(diff_peaks_oligo)

# GRN
All_TFBS <- read.csv("sub_median_1_DEG_final.csv")

All_TFBS <- All_TFBS %>%
  filter(!duplicated(peak)) %>%
  separate(peak,into=c("chr","coord"),sep=":") %>%
  separate(coord,into=c("start","end"),sep="-")

All_TFBS$start <- as.numeric(All_TFBS$start)
All_TFBS$end <- as.numeric(All_TFBS$end)

# chr (function)
total_SNP <- list()
signi_SNP <- list()

oligo_diff_SNP <- list()

oligo_all_SNP <- list()
oligo_all_eGene <- list()
oligo_all_TG <- list()
check <- list()

chromatin_inter_SNP <- list()

for (i in 1:22){

  print(i)
  
  chr_df <- read.table(paste0("Oligodendrocytes_chr",i))
  colnames(chr_df) <- c("eGene","SNP_id","distance","p","beta")
  
  total_SNP[[i]] <- length(unique(chr_df$SNP_id))
  
  chr_df_0.05 <- chr_df %>%
    filter(p<0.05)
  
  signi_SNP[[i]] <- length(unique(chr_df_0.05$SNP_id))
  
  SNP_posi_chr <- SNP_posi %>%
    filter(chr_SNP == paste0("chr",i))
  SNP_posi_chr <- SNP_posi_chr[,c("SNP","chr_SNP","position")]
  colnames(SNP_posi_chr)[1] <- "SNP_id"
  
  chr_map <- merge(chr_df_0.05[,c(1,2,3)],SNP_posi_chr, by="SNP_id")
  
  chr_map$position <- as.numeric(chr_map$position)
  isnps <- with(chr_map, IRanges(position, width=1, names=SNP_id))
  
  # For oligo-diff
  diff_chr <- diff_peaks_oligo %>%
    filter(chr == paste0("chr",i))
  
  diff_chr$start <- as.numeric(diff_chr$start)
  diff_chr$end <- as.numeric(diff_chr$end)
  diff_chr$id <- paste0(diff_chr$chr,":",diff_chr$start,"_",diff_chr$end)
  diff_peak <- with(diff_chr, IRanges(start, end),names=id)
  
  olaps_diff <- findOverlaps(isnps, diff_peak)
  chr_eQTLs_diff <- cbind(chr_map[queryHits(olaps_diff),], diff_chr[subjectHits(olaps_diff),])
  
  chr_eQTLs_diff$comp_ID <- paste0(chr_eQTLs_diff$SNP_id,"&",chr_eQTLs_diff$eGene)
  
  chr_eQTLs_diff <- chr_eQTLs_diff%>%
    filter(!duplicated(comp_ID))
  dim(chr_eQTLs_diff)
  
  chr_eQTLs_diff_df <- chr_eQTLs_diff %>%
    separate(eGene, into=c("eGene","junk"),sep="_")
  
  chr_eQTLs_diff_df <- chr_eQTLs_diff_df[,c(1,2,5,6,10,11,12,13,14,15)]

  oligo_diff_SNP[[i]] <- length(unique(chr_eQTLs_diff_df$SNP_id))
  
  write.csv(chr_eQTLs_diff_df, file=paste0("diff_chr",i,".csv"),row.names = F)
 
  # For peak-gene & oligo-diff
  All_TFBS0 <- All_TFBS %>%
    filter(chr == paste0("chr",i))

  All_TFBS0$start <- as.numeric(All_TFBS0$start)
  All_TFBS0$end <- as.numeric(All_TFBS0$end)
  All_peak <- with(All_TFBS0, IRanges(start, end))
  
  olaps_fidd <- findOverlaps(isnps, All_peak)
  chr_eQTLs <- cbind(chr_map[queryHits(olaps_fidd),], All_TFBS0[subjectHits(olaps_fidd),])
  
  chr_eQTLs$comp_ID <- paste0(chr_eQTLs$SNP_id,"&",chr_eQTLs$eGene)
  
  chr_eQTLs <- chr_eQTLs%>%
    filter(!duplicated(comp_ID))
  dim(chr_eQTLs)
  
  chr_eQTLs_df <- chr_eQTLs %>%
    separate(eGene, into=c("eGene","junk"),sep="_")
  
  chr_eQTLs_all_summary <- chr_eQTLs_df[,c(1,2,5,6,10,11,12,13,14,15,16,17,19)]
  chr_eQTLs_all_summary$check <- ifelse(chr_eQTLs_all_summary$eGene == chr_eQTLs_all_summary$TG,"same","different")
  
  oligo_all_SNP[[i]] <- length(unique(chr_eQTLs_df$SNP_id))
  oligo_all_eGene[[i]] <- length(unique(chr_eQTLs_df$eGene))
  oligo_all_TG[[i]] <- length(unique(chr_eQTLs_df$TG))
  check[[i]] <- table(chr_eQTLs_all_summary$check)[2]
  
  write.csv(chr_eQTLs_all_summary, file=paste0("ALL_chr",i,".csv"),row.names = F)
  
}

summary <- rbind.data.frame(unlist(total_SNP), unlist(signi_SNP),
                            unlist(oligo_diff_SNP), 
                            unlist(oligo_all_SNP),unlist(oligo_all_eGene), unlist(oligo_all_TG), unlist(check))

summary <- t(summary)
dim(summary)
head(summary)

colnames(summary) <- c("total_SNPs","signi_SNPs","oligo_specific_SNPs","oligo_peak_gene_SNPs",
                       "eGenes","TGs","check","chromatin_interaciton_SNPs")

rownames(summary) <- paste0(rep("chr",22),1:22)

write.csv(summary, "eQTL_summary_0415124.csv",row.names = T)


## Combine files
oligo_specific_eQTLs <- list.files(path="Mapping_oligo_specific/", pattern=".csv", all.files=TRUE, 
                                        full.names=TRUE)
oligo_specific_eQTLs <- mixedsort(oligo_specific_eQTLs)

oligo_peak_gene_eQTLs <- list.files(path="Mapping_oligo_peak_gene/", pattern=".csv", all.files=TRUE, 
                                   full.names=TRUE)
oligo_peak_gene_eQTLs <- mixedsort(oligo_peak_gene_eQTLs)

# Combine the summaries
oligo_specific_eQTLs_total = data.frame()
for (i in 1:22){
  print(paste0("i=",i))
  df <- assign(paste0("mat_",i), read.csv(oligo_specific_eQTLs[i],header = T))
  oligo_specific_eQTLs_total <- rbind(oligo_specific_eQTLs_total,df)
}

oligo_peak_gene_eQTLs_total = data.frame()
for (i in 1:22){
  print(paste0("i=",i))
  df <- assign(paste0("mat_",i), read.csv(oligo_peak_gene_eQTLs[i],header = T))
  oligo_peak_gene_eQTLs_total <- rbind(oligo_peak_gene_eQTLs_total,df)
}


# oligo_peak_gene_eQTLs_total - check TG==eGene
oligo_peak_gene_eQTLs_total$peak <- paste0(oligo_peak_gene_eQTLs_total$chr,":",
                                           oligo_peak_gene_eQTLs_total$start,"-",
                                           oligo_peak_gene_eQTLs_total$end)

All_TFBS <- read.csv("sub_median_1_DEG_final.csv")

full_eQTL_mapping <- merge(oligo_peak_gene_eQTLs_total[,c(1,2,3,4,15)],All_TFBS[,-c(1:3)], by="peak")
head(full_eQTL_mapping)
dim(full_eQTL_mapping)

full_eQTL_mapping$id <- paste0(full_eQTL_mapping$peak,"&",full_eQTL_mapping$SNP_id,"&",
                               full_eQTL_mapping$eGene,"&",full_eQTL_mapping$chr_SNP,"&",full_eQTL_mapping$position,"&",
                               full_eQTL_mapping$TG,"&", full_eQTL_mapping$TF1,"_",full_eQTL_mapping$TF2)

full_eQTL_mapping_df <- full_eQTL_mapping %>%
  filter(!duplicated(id))
dim(full_eQTL_mapping_df)

full_eQTL_mapping_df$check <- ifelse(full_eQTL_mapping_df$eGene == full_eQTL_mapping_df$TG,
                                     "same","different")

table(full_eQTL_mapping_df$check)

oligo_peak_gene_eQTLs_total_yes <- full_eQTL_mapping_df[full_eQTL_mapping_df$check == "same",]
dim(oligo_peak_gene_eQTLs_total_yes)
head(oligo_peak_gene_eQTLs_total_yes)
length(unique(oligo_peak_gene_eQTLs_total_yes$SNP_id))
length(unique(oligo_peak_gene_eQTLs_total_yes$TG))
eQTL_TGs <- unique(oligo_peak_gene_eQTLs_total_yes$TG)
write.csv(eQTL_TGs, "eQTL_TGs.csv",row.names = F)
length(unique(oligo_peak_gene_eQTLs_total_yes$peak))

oligo_peak_gene_eQTLs_total_yes_Key <- oligo_peak_gene_eQTLs_total_yes[oligo_peak_gene_eQTLs_total_yes$TF1 %in% key_TF_list1 |
                                                                         oligo_peak_gene_eQTLs_total_yes$TF2 %in% key_TF_list1,]
dim(oligo_peak_gene_eQTLs_total_yes_Key)
head(oligo_peak_gene_eQTLs_total_yes_Key)
length(unique(oligo_peak_gene_eQTLs_total_yes_Key$SNP_id))
length(unique(oligo_peak_gene_eQTLs_total_yes_Key$TG))
eQTL_key_TGs <- unique(oligo_peak_gene_eQTLs_total_yes_Key$TG)
write.csv(eQTL_key_TGs, "eQTL_key_TGs.csv",row.names = F)
length(unique(oligo_peak_gene_eQTLs_total_yes_Key$peak))


# ref TGs
overlaps_genes <- read.csv("overlaps_genes.csv")

unique(oligo_peak_gene_eQTLs_total_yes_Key$TG)[unique(oligo_peak_gene_eQTLs_total_yes_Key$TG) %in% overlaps_genes$x]

length(unique(oligo_peak_gene_eQTLs_total_yes_Key$TG)) # 16%

unique(oligo_peak_gene_eQTLs_total_yes_Key$TF_id)


oligo_peak_gene_eQTLs_total_yes_Key_refTG <- oligo_peak_gene_eQTLs_total_yes_Key[oligo_peak_gene_eQTLs_total_yes_Key$TG %in% overlaps_genes$x,]
dim(oligo_peak_gene_eQTLs_total_yes_Key_refTG)

unique(oligo_peak_gene_eQTLs_total_yes_Key_refTG$TG)

full_effect_df_filter <-read.csv("interaction_summary_median_1_CoV_0.5.csv")

oligo_peak_gene_eQTLs_Key_rank <- merge(oligo_peak_gene_eQTLs_total_yes_Key, full_effect_df_filter, by="TF_id")
hist(oligo_peak_gene_eQTLs_Key_rank$percentile)


## non key
oligo_peak_gene_eQTLs_total_yes_non_Key <- oligo_peak_gene_eQTLs_total_yes[!(oligo_peak_gene_eQTLs_total_yes$TF1 %in% key_TF_list1 |
                                                                         oligo_peak_gene_eQTLs_total_yes$TF2 %in% key_TF_list1),]

oligo_peak_gene_eQTLs_non_Key_rank <- merge(oligo_peak_gene_eQTLs_total_yes_non_Key, full_effect_df_filter, by="TF_id")
hist(oligo_peak_gene_eQTLs_non_Key_refTG_rank$percentile)

## t test
t.test(oligo_peak_gene_eQTLs_Key_rank$percentile,
       oligo_peak_gene_eQTLs_non_Key_rank$percentile)

oligo_peak_gene_eQTLs_Key_rank$type <- "key"
oligo_peak_gene_eQTLs_non_Key_rank$type <- "non_key"


combine_eQTL_key_nonkey <- rbind.data.frame(oligo_peak_gene_eQTLs_Key_rank,oligo_peak_gene_eQTLs_non_Key_rank)

library(plyr)
mu <- ddply(combine_eQTL_key_nonkey, "type", summarise, grp.mean=median(percentile))

ggplot(combine_eQTL_key_nonkey, aes(x=percentile, color=type, fill=type)) +
  geom_histogram(aes(y=..density..), position="identity",,alpha=0.5)+
  #geom_density(alpha=0.6)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),
             linetype="dashed")+
  scale_color_manual(values=c("#56B4E9","purple"),labels=c('Cooperative TF pairs (Key)', 'co-binding TF pairs (non-Key)'))+
  scale_fill_manual(values=c("#56B4E9","purple"),labels=c('Cooperative TF pairs (Key)', 'co-binding TF pairs (non-Key)'))+
  #labs(title="Cooperativity of co-binding TFs",x="Interaction scores (percentile)", y = "Density")+
  theme(axis.text.y=element_text(size=20), axis.text.x=element_text(size=20),
        axis.title=element_text(size=20), plot.title = element_text(size=20),
        legend.text=element_text(size=12), legend.title =element_blank(),
        legend.direction = "horizontal",legend.position="top")+
  xlab("Percentile (Interaction score)") + ylab("density")


full_yes_df <- read.csv("/median1_eQTLs_full_yes.csv")
length(unique(full_yes_df$TG))

full_yes_key_df <- read.csv("edian1_eQTLs_full_yes_keyTF.csv")
length(unique(full_yes_key_df$TG))


