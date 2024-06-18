library(dplyr)

key_sub <- read.csv("key_sub.csv")
head(key_sub)
key_sub$TF_id <- paste0(key_sub$TF1,"_",key_sub$TF2)

summary(key_sub$rank)

key_sub_high <- key_sub %>%
  filter(rank < 15)
dim(key_sub_high)

All_TG_TF1_TF2_df <- read.csv("All_TG_TF1_TF2_df_final.csv")
head(All_TG_TF1_TF2_df)
All_TG_TF1_TF2_df$TF_id <- paste0(All_TG_TF1_TF2_df$TF1,"_",All_TG_TF1_TF2_df$TF2)

key_sub_df_high <- All_TG_TF1_TF2_df[All_TG_TF1_TF2_df$TF_id %in% key_sub_high$TF_id,]
dim(key_sub_df_high)
head(key_sub_df_high)

write.csv(key_sub_df_high[,c(1:3)], "key_sub_df_high.csv")


# Load full GRN
GRN_full <- read.csv("/sub_median_1_DEG_final.csv")
# Load Interactions
key_CoV_df <- read.csv("/Interaction_summary/key_CoV.csv")
non_key_CoV_df <- read.csv("non_key_CoV.csv")

# Sub GRNs
key_CoV_df <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean))

key_CoV_NKX6.2 <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="NKX6.2" | TF2 == "NKX6.2")
dim(key_CoV_NKX6.2)

key_CoV_SOX10 <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="SOX10" | TF2 == "SOX10")
dim(key_CoV_SOX10)

key_CoV_MYRF <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="MYRF" | TF2 == "MYRF")
dim(key_CoV_MYRF)

key_CoV_OLIG1 <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="OLIG1" | TF2 == "OLIG1")
dim(key_CoV_OLIG1)

key_CoV_OLIG2 <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="OLIG2" | TF2 == "OLIG2")
dim(key_CoV_OLIG2)

key_CoV_TCF7L2 <- key_CoV_df %>%
  arrange(desc(Interaction_score_mean)) %>%
  filter(TF1=="TCF7L2" | TF2 == "TCF7L2")
dim(key_CoV_TCF7L2)

intersect(intersect(key_CoV_NKX6.2$TG,key_CoV_SOX10$TG),key_CoV_MYRF$TG)
intersect(key_CoV_SOX10$TG,key_CoV_OLIG2$TG)

GRN_SOX10_OLIG2 <- GRN_full %>%
  filter(TF1 %in% c("SOX10","OLIG2") | TF2 %in% c("SOX10","OLIG2"))
dim(GRN_SOX10_OLIG2)

GNR_top10 <- GRN_full[GRN_full$TF_id %in% key_CoV_df$TF_id[1:10],]
dim(GNR_top10)
head(GNR_top10)

GNR_top10_Cyto <- data.frame(GNR_top10[,c("TG","TF1","TF2")])
from <- as.character(rbind(GNR_top10_Cyto$TF1,GNR_top10_Cyto$TF2))
to <- rep(GNR_top10_Cyto$TG,2)

GNR_top10_Cyto_df <- data.frame(from,to)
GNR_top10_Cyto_df$weight <- 1
GNR_top10_Cyto_df$type <- "hyperlink"

write.csv(GNR_top10_Cyto_df,"GNR_top10_Cyto_df.csv",row.names = F)

# top 9 for each key TF
GNR_top_each <- rbind(key_CoV_NKX6.2[1:9,],key_CoV_SOX10[1:9,],key_CoV_MYRF[1:9,],
                      key_CoV_OLIG1[1:9,],key_CoV_OLIG2[1:9,],key_CoV_TCF7L2[1:9,])
GNR_top_each_df <-GRN_full[GRN_full$TF_id %in% GNR_top_each$TF_id,]
dim(GNR_top_each_df)
head(GNR_top_each_df)
unique(c(GNR_top_each_df$TF1,GNR_top_each_df$TF2))
unique(GNR_top_each_df$TG)
unique(GNR_top_each_df$TG[duplicated(GNR_top_each_df$TG)])

GNR_each_Cyto <- data.frame(GNR_top_each_df[,c("TG","TF1","TF2")])
from <- as.character(rbind(GNR_each_Cyto$TF1,GNR_each_Cyto$TF2))
to <- rep(GNR_each_Cyto$TG,2)

GNR_each_Cyto_df <- data.frame(from,to)
GNR_each_Cyto_df$weight <- 1
GNR_each_Cyto_df$type <- "hyperlink"

write.csv(GNR_each_Cyto_df,"/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Manuscript/GRNs/GRN_Cyto/GNR_each_Cyto_df.csv",row.names = F)
GNR_each_Cyto_df <- read.csv("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Manuscript/GRNs/GRN_Cyto/GNR_each_Cyto_df.csv")


## TG overlaps
SOX10_TCF12 <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="SOX10_TCF12",]$TG)
TCF7L2_RBPJ <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="TCF7L2_RBPJ",]$TG)
RORA_OLIG2 <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="RORA_OLIG2",]$TG)
FOXP1_NKX6.2 <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="FOXP1_NKX6.2",]$TG)
FOXP1_OLIG1 <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="FOXP1_OLIG1",]$TG)
SMAD2_MYRF <- c(GNR_top_each_df[GNR_top_each_df$TF_id=="SMAD2_MYRF",]$TG)

co_TF_list <- list(SOX10_TCF12,TCF7L2_RBPJ,RORA_OLIG2,FOXP1_NKX6.2,FOXP1_OLIG1,
                   SMAD2_MYRF)

intersect(SOX10_TCF12,TCF7L2_RBPJ)
#  [1] "AMOTL2"  "AMPD3"   "BOK"     "CECR2"   "EVI5L"   "FAM178B" "FGFR2"   "MATN2"  
# [9] "OLIG2"   "SLC45A3" "SYNGR2"  "ADA"     "INF2"    "OLIG1"   "SOX13"   "TMEM230"
intersect(SOX10_TCF12,RORA_OLIG2)
# "ANO4"  "CALD1"
intersect(SOX10_TCF12,FOXP1_OLIG1)
# [1] "FA2H"  "CALD1"
intersect(SOX10_TCF12,SMAD2_MYRF)
# [1] "JAKMIP3" "POU2AF1" "SEMA6A" 
intersect(SOX10_TCF12,FOXP1_NKX6.2)
# [1] "CPM"     "OLIG2"   "PLLP"    "ZNF536"  "CHD7"    "OLIG1"   "ST3GAL5"


overlaps_genes <- read.csv("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Final_pipline/Step6/overlaps_genes.csv")

intersect(overlaps_genes$x,intersect(SOX10_TCF12,TCF7L2_RBPJ))
#[1] "SOX13"
intersect(overlaps_genes$x,intersect(SOX10_TCF12,RORA_OLIG2))
intersect(overlaps_genes$x,intersect(SOX10_TCF12,FOXP1_OLIG1))
#[1] "FA2H"
intersect(overlaps_genes$x,intersect(SOX10_TCF12,SMAD2_MYRF))
intersect(overlaps_genes$x,intersect(SOX10_TCF12,FOXP1_NKX6.2))
#[1] "PLLP"    "ST3GAL5"

intersect(overlaps_genes$x,intersect(SMAD2_MYRF,TCF7L2_RBPJ))
intersect(overlaps_genes$x,intersect(SMAD2_MYRF,RORA_OLIG2))
intersect(overlaps_genes$x,intersect(SMAD2_MYRF,FOXP1_OLIG1))
intersect(overlaps_genes$x,intersect(SMAD2_MYRF,FOXP1_NKX6.2))

intersect(overlaps_genes$x,intersect(TCF7L2_RBPJ,RORA_OLIG2))
intersect(overlaps_genes$x,intersect(TCF7L2_RBPJ,FOXP1_OLIG1))
intersect(overlaps_genes$x,intersect(TCF7L2_RBPJ,FOXP1_NKX6.2))

intersect(overlaps_genes$x,intersect(FOXP1_OLIG1,RORA_OLIG2))
intersect(overlaps_genes$x,intersect(FOXP1_OLIG1,FOXP1_NKX6.2))

intersect(overlaps_genes$x,intersect(RORA_OLIG2,FOXP1_NKX6.2))

"CALD1" %in% overlaps_genes$x
"PPP1R16B" %in% overlaps_genes$x

library(tidyverse)
l <- mget(ls(pattern = '^set\\d'))

overlaps <- map(2:length(l),
    ~ combn(l, .x, \(x)
            list(reduce(x, intersect)) %>%
              set_names(str_c(names(
                x
              ), collapse = ' & ')),
            simplify = FALSE)) %>%
  unlist(FALSE) %>%
  unlist(FALSE) %>%
  c(.,
    map(seq_along(l), ~reduce(l[-.x], setdiff,.init = l[[.x]])) %>% 
      set_names(names(l))
  ) %>% 
  keep(~ length(.x) > 0)


set1[overlaps$set1 %in% intersect(SOX10_TGs,SOX10_ref_genes)]
set6[overlaps$set6 %in% intersect(MYRF_TGs,MYRF_ref_brain_genes)]
set6[overlaps$set6 %in% intersect(MYRF_TGs,MYRF_ref_GE_genes)]
set5[overlaps$set5 %in% intersect(OLIG1_TGs,OLIG1_ref_GE_genes)]
set3[overlaps$set3 %in% intersect(OLIG2_TGs,OLIG2_ref_GE_genes)]


## Non-key oligo
Interaction_non_keyTFs <- read.csv("/Users/jeromechoi/Documents/jerome/Documents/WISC/BMI/Daifeng Wang/Oligo Project/Pipeline/chromatin_interaction/Final_pipline/Step2/Interaction_summary/Interaction_non_keyTFs.csv")
dim(Interaction_non_keyTFs)

Interaction_non_keyTFs <- Interaction_non_keyTFs %>%
  arrange(desc(Interaction_score))

Interaction_non_keyTFs_top <- Interaction_non_keyTFs[1:55,]
Interaction_non_keyTFs_top



