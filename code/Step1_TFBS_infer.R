library(TFBSTools);library(JASPAR2022)

loadRData <- function(file_name){
  load(file_name)
  get(ls()[ls() != "file_name"])
}

TFBS_infer <- function(df, database = JASPAR2022::JASPAR2022, species_type = 9606, min_score = 0.7,
                          pwm_type = 'prob',num_cores = 2){
  
  opts <- list()
  opts[["species"]] <- species_type
  opts[["all_versions"]] <- TRUE
  PFMatrixList <-loadRData("/JASPAR2022_111624.RData")
  pwmlist <- TFBSTools::toPWM(PFMatrixList, type = pwm_type)
  TF_names <- TFBSTools::name(pwmlist)
  names(TF_names) = NULL
  
  TF_names_splited = sapply(TF_names,data.table::tstrsplit,'::|\\(var.2\\)|\\(var.3\\)')
  
  df$promoter_id <- paste(df$gene_chr,':',df$promoter_start,'-',df$promoter_end,sep = '')
  df = data.table::data.table(df)
  
  df_p <- df[,c('gene_chr','promoter_start','promoter_end','promoter_id')]
  
  suppressWarnings( G1 <- GenomicRanges::GRanges(seqnames = df_p$gene_chr,
                                                 IRanges::IRanges(start=df_p$promoter_start,
                                                                  end=df_p$promoter_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:24]))
  G1 <- GenomicRanges::trim(G1)
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  `%dopar%` <- foreach::`%dopar%`
  df_p$promoter_TF <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G1[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg38",
                                                                 out = "scores"
                                         )
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))
                                         
                                       }
  parallel::stopCluster(cl)
  
  
  df$promoter_TF <- df_p$promoter_TF[match(df_p$promoter_id,df$promoter_id)]

  df <- df[, c('gene','promoter_id',
               'promoter_TF')]
  colnames(df) <- c('gene','promoter',
                    'promoter_TF')
  return(df)
  
  
}

Group1 <- read.csv("Group1_df.csv")
Group1 <- Group1[,c(3,4,5,6)]
colnames(Group1) <- c("gene","gene_chr","promoter_start","promoter_end")

Group1$promoter_start <- as.numeric(Group1$promoter_start)
Group1$promoter_end <- as.numeric(Group1$promoter_end)

Group1_TFBS <- TFBS_infer(Group1)

save(Group1_TFBS, file="Group1_TFBS_011424.RData")

