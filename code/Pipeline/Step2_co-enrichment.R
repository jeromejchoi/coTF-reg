library(motifmatchr)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(readxl)

loadRData <- function(file_name){
  load(file_name)
  get(ls()[ls() != "file_name"])
}


load("JASPAR2022_030724.RData")

PFMatrixList <- JASPAR2022_030724
pwm_type <- "prob"
pwmlist <- TFBSTools::toPWM(PFMatrixList, type = pwm_type)
TF_names <- TFBSTools::name(pwmlist)
names(pwmlist) = TF_names

names(PFMatrixList) <- TF_names

TF_names_splited = sapply(TF_names,data.table::tstrsplit,'::|\\(var.2\\)|\\(var.3\\)')

# Load the peaks
load("All_TFBS_input.RData")

All_TFBS$start <- as.numeric(All_TFBS$start)
All_TFBS$end <- as.numeric(All_TFBS$end)

TFBSs <- GRanges(seqnames = All_TFBS$chr,
                 IRanges::IRanges(start=All_TFBS$start,
                                  end=All_TFBS$end))
                 #seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:24])

peak_gr <- GRanges(seqnames = All_TFBS$chr,
                 IRanges::IRanges(start=All_TFBS$start,
                                  end=All_TFBS$end))

# Get motif matches for example motifs in peaks
motif_ix <- matchMotifs(pwmlist, TFBSs,
                        genome = "hg38",
                        out = "positions")
motif_ix_df <- data.frame(motif_ix)

motif_ix_df_gr <- GRanges(seqnames = motif_ix_df$seqnames,
                          IRanges::IRanges(start=motif_ix_df$start,
                                           end=motif_ix_df$end))

motif_ix_overlaps <- findOverlaps(query = motif_ix_df_gr, subject = peak_gr, type = "any")
motif_ix_overlaps

motif_ix_overlaps_df <- data.frame(motif_ix_df[queryHits(motif_ix_overlaps),], All_TFBS[subjectHits(motif_ix_overlaps),])
dim(motif_ix_overlaps_df)
motif_ix_overlaps_df$id <- paste0(motif_ix_overlaps_df$seqnames,":",motif_ix_overlaps_df$start,"-",
                                  motif_ix_overlaps_df$end)

motif_ix_overlaps_df$peak <- paste0(motif_ix_overlaps_df$chr,":",motif_ix_overlaps_df$start.1,"-",motif_ix_overlaps_df$end.1)
length(unique(motif_ix_overlaps_df$peak))
#motif_ix_overlaps_df <- motif_ix_overlaps_df %>%
#  separate(group_name, into=c("TF","junk"),sep=" ")
Step1_overlaps_motif_peak <- motif_ix_overlaps_df
save(Step1_overlaps_motif_peak,file="Step2_overlaps_motif_peak.RData")
write.csv(Step1_overlaps_motif_peak[,-13],"Step2_overlaps_motif_peak.csv",row.names = F)


motif_df <- data.frame(motif_ix)
head(motif_df)
dim(motif_df)
length(unique(motif_df$group_name))


motif_df <- motif_df[,c(2,3,4,5)]


TF_pair <- data.frame(t(combn(unique(TF_names),2)))
colnames(TF_pair) <- c("TF1","TF2")

TF_pair$id <- paste0(TF_pair$TF1,"_",TF_pair$TF2)
length(unique(TF_pair$id))

TF_pair <- TF_pair %>%
  filter(!duplicated(id))

TFBS_over_list <- list()
peak_over_list <- list()
peak_list <- list()
hypertest <- list()

n = length(c(TF_pair$TF1,TF_pair$TF2)) * (length(c(TF_pair$TF1,TF_pair$TF2))-1) * 0.5
for (i in 1:n){
  print(i)
  # Exact TFBS overlaps
  TF1_motif_df <- motif_df[motif_df$group_name == TF_pair$TF1[i],]
  TF2_motif_df <- motif_df[motif_df$group_name == TF_pair$TF2[i],]
  
  TF1_df_G <- GRanges(seqnames = TF1_motif_df$seqnames,
                      IRanges::IRanges(start=TF1_motif_df$start,
                                       end=TF1_motif_df$end))
  TF2_df_G <- GRanges(seqnames = TF2_motif_df$seqnames,
                      IRanges::IRanges(start=TF2_motif_df$start,
                                       end=TF2_motif_df$end))
  
  TFBS_overlaps <- findOverlaps(query = TF1_df_G, subject = TF2_df_G, type = "any")
  TFBS_over_list[[i]] <- length(TFBS_overlaps)
  
  # Peak overlaps
  TF1_peak_df <- motif_ix_overlaps_df[motif_ix_overlaps_df$group_name == TF_pair$TF1[i],]
  TF1_peak_df$id <- paste0(TF1_peak_df$chr,":",TF1_peak_df$start.1,"-",TF1_peak_df$end.1)
  TF2_peak_df <- motif_ix_overlaps_df[motif_ix_overlaps_df$group_name == TF_pair$TF2[i],]
  TF2_peak_df$id <- paste0(TF2_peak_df$chr,":",TF2_peak_df$start.1,"-",TF2_peak_df$end.1)
  
  peak_over_list[[i]] <- length(unique(intersect(TF1_peak_df$id,TF2_peak_df$id)))
  peak_list[[i]] <- unique(intersect(TF1_peak_df$id,TF2_peak_df$id))
  
  # Hypergeometric test
  hyper <- phyper(length(unique(intersect(TF1_peak_df$id,TF2_peak_df$id)))-1,
                  length(unique(TF2_peak_df$id)),
                  length(unique(motif_ix_overlaps_df$peak)),
                  length(unique(TF1_peak_df$id)),
                  lower.tail=FALSE)
  
  hypertest[[i]] <- hyper
}
summary(unlist(TFBS_over_list))
hist(unlist(TFBS_over_list))
summary(unlist(peak_over_list))
hist(unlist(peak_over_list))
summary(unlist(hypertest))
hist(unlist(hypertest))

summary <- data.frame(unlist(TFBS_over_list),unlist(peak_over_list),unlist(hypertest),unlist(peak_list))
dim(summary)

rownames(summary) <- paste0(TF_pair$TF1,"_",TF_pair$TF2)
colnames(summary) <- c("num_TFBSs","num_peaks","p","peak_list")

## Add FDR adj
summary$p_fdr <- p.adjust(summary$p,method="fdr",n=length(summary$p))

summary <- summary[,c(1,2,4,3,5)]

write.csv(summary, "Step2_co-enrichment_summary.csv")
