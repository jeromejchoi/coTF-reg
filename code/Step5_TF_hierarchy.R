
library(ggplot2)

# Load GRN
GRN_full <- read.csv("sub_median_1_DEG_final.csv")
TFs <- unique(c(GRN_full$TF1,GRN_full$TF2))
TGs <- unique(c(GRN_full$TG))

unique(GRN_full$TG)[unique(GRN_full$TG) %in% unique(c(GRN_full$TF1,GRN_full$TF2))]


in_degree <- list()
for (i in 1:445){
  print(i)
  TF_sub <- GRN_full[GRN_full$TG == TGs[i], ]
  in_TF <- length(unique(c(TF_sub$TF1,TF_sub$TF2)))
  in_degree[[i]] <- in_TF
}
in_degree_df <- data.frame(TGs,unlist(in_degree))
colnames(in_degree_df) <- c("gene","in_degree")

out_degree <- list()
for (i in 1:206){
  print(i)
  TG_sub <- GRN_full[GRN_full$TF1 == TFs[i] | GRN_full$TF2 == TFs[i], ]
  out_TG <- length(unique(TG_sub$TG))
  
  out_degree[[i]] <- out_TG
}
out_degree_df <- data.frame(TFs,unlist(out_degree))
colnames(out_degree_df) <- c("gene","out_degree")

df=merge(out_degree_df,in_degree_df,by="gene")
dim(df)
head(df)

df$h_metric <- (df$out_degree-df$in_degree)/(df$out_degree+df$in_degree)
df$profile <- ifelse(df$h_metric >0.33,"top",ifelse(df$h_metric >-0.33,"middle","bottom"))
table(df$profile)

df$gene <- italic(df$gene)

# Plot in/out degrees
ggplot(df, aes(x=out_degree, y=in_degree, color=profile, size=h_metric))+ geom_point(aes(fill=profile),shape = 21,color = "black") + geom_text(label=df$gene, position = position_dodge(width = 1),
            vjust = 2, size = 4,colour="black", fontface=3) + 
  scale_fill_manual(values=c('#E0E0E0','#F5F076', 'orange'))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  xlab("Out-degree") + ylab("In-degree")

top_middle <- expand.grid(df$gene[df$profile=="top"], df$gene[df$profile=="middle"])
colnames(top_middle) <- c("from","to")
middle_bottom <- expand.grid(df$gene[df$profile=="middle"], df$gene[df$profile=="bottom"])
colnames(middle_bottom) <- c("from","to")

top_middle_bottom <- rbind.data.frame(top_middle,middle_bottom)
top_middle_bottom$weight <- 1
top_middle_bottom$type <- "hyperlink"

write.csv(top_middle_bottom, "top_middle_bottom.csv", row.names = F)

