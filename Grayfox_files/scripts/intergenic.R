cf_gf_af_fst <- read_tsv("~/Downloads/cf4_gf_af_window_phyloP.tsv") %>% 
  dplyr::rename("gf_chr"="seqnames","gf_st"="start","gf_end"="end") %>% 
  select(gf_chr,gf_st,gf_end,af_chr,af_st,af_end,outlier_af,cf4_chr,cf4_st,cf4_end,outlier_cf4) %>% 
  mutate(gf_chr = as.character(gsub("^chr", "", gf_chr)), af_chr = as.character(gsub("^chr", "", af_chr)),cf4_chr = as.character(gsub("^chr", "", cf4_chr))) 


gf_genic <- fread("~/Desktop/Foxes/DATA/grayFox.gtf.gz",col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
af_genic <- fread("~/Desktop/Foxes/DATA/arcticFox.ncbiRefSeq.gtf.gz",col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

gfchroms <- read_tsv("~/Desktop/Foxes/DATA/gf_chroms.txt", col_names = c("scaf","chrom")) 
afchroms <- read_tsv("~/Desktop/Foxes/DATA/arc_chroms.txt", col_names = c("chrom","name","scaff")) %>% select(chrom,scaff) 


gf_genic %>% left_join(gfchroms,by=c("seqname"="scaf")) %>% 
  relocate(chrom) %>% select(-seqname) %>% drop_na(chrom) %>% write_tsv("~/Desktop/Foxes/DATA/grayfox.chr.gtf",col_names = F,quote = "none")

af_genic %>% left_join(afchroms,by=c("seqname"="scaff")) %>% 
  relocate(chrom) %>% select(-seqname) %>% drop_na(chrom) %>% write_tsv("~/Desktop/Foxes/DATA/arcticfox.chr.gtf",col_names = F,quote = "none")

fread("~/Desktop/Foxes/DATA/canFam4.ncbiRefSeq.gtf.gz",col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))%>%
  filter(!grepl("^chr[^0-9]", seqname)) %>% dplyr::rename("chrom"="seqname")%>% write_tsv("~/Desktop/Foxes/DATA/canfam4.chr.gtf",col_names = F,quote = "none")


gfanf <- gf_genic %>% left_join(gfchroms,by=c("seqname"="scaf")) %>% relocate(chrom) %>% select(-seqname) %>%
  mutate(chrom = as.character(gsub("^chr", "", chrom)),start = as.numeric(start), end = as.numeric(end)) %>% drop_na(chrom) 

afanf <- af_genic %>% left_join(afchroms,by=c("seqname"="scaff")) %>% relocate(chrom) %>% select(-seqname) %>%
  mutate(chrom = as.character(gsub("^chr", "", chrom)),start = as.numeric(start), end = as.numeric(end))%>% drop_na(chrom) 

cfanf <- fread("~/Desktop/Foxes/DATA/canFam4.ncbiRefSeq.gtf.gz",col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))%>%
  filter(!grepl("^chr[^0-9]", seqname)) %>% dplyr::rename("chrom"="seqname")%>%
  mutate(chrom = as.character(gsub("^chr", "", chrom)),start = as.numeric(start), end = as.numeric(end))



# Canfam4 genic regions
cf_gr <- GRanges(seqnames = cfanf$chrom, 
  ranges = IRanges(start = cfanf$start, end = cfanf$end))

# Gray fox genic regions
gf_gr <- GRanges(seqnames = gfanf$chrom, 
  ranges = IRanges(start = gfanf$start, end = gfanf$end))

# Arctic fox genic regions
af_gr <- GRanges(seqnames = afanf$chrom, 
  ranges = IRanges(start = afanf$start, end = afanf$end))

cf_gf_af_fst_gr <- GRanges(
  seqnames = as.character(cf_gf_af_fst$gf_chr),  # Use Gray fox as the reference coordinate system
  ranges = IRanges(start = cf_gf_af_fst$gf_st, end = cf_gf_af_fst$gf_end))



compute_region_counts <- function(windows_gr, genic_gr) {
  # Find overlaps between FST windows and genic regions
  overlaps <- findOverlaps(windows_gr, genic_gr)
  
  # Create an empty results dataframe
  results <- data.frame(
    window_id = seq_along(windows_gr),  # Window index
    genic_count = 0,  # Initialize genic count
    intergenic_count = 0  # Initialize intergenic count
  )
  
  # If there are overlaps, count the number of unique genic regions per window
  if (length(overlaps) > 0) {
    genic_counts <- table(queryHits(overlaps))  # Count genic regions per FST window
    results$genic_count[as.numeric(names(genic_counts))] <- as.numeric(genic_counts)
  }
  
  # Intergenic count = 1 if no genic regions are assigned to the window
  results$intergenic_count <- ifelse(results$genic_count == 0, 1, 0)
  
  return(results)
}


cf_results <- compute_region_counts(cf_gf_af_fst_gr, cf_gr)  # Canfam4
gf_results <- compute_region_counts(cf_gf_af_fst_gr, gf_gr)  # Gray fox
af_results <- compute_region_counts(cf_gf_af_fst_gr, af_gr)  # Arctic fox



cf_gf_af_fst <- cf_gf_af_fst %>%
  mutate(
    genic_cf = cf_results$genic_count,
    intergenic_cf = cf_results$intergenic_count,
    genic_gf = gf_results$genic_count,
    intergenic_gf = gf_results$intergenic_count,
    genic_af = af_results$genic_count,
    intergenic_af = af_results$intergenic_count
  )

cf_gf_af_fst_long <- cf_gf_af_fst %>%
  pivot_longer(cols = c(genic_cf, intergenic_cf, genic_gf, intergenic_gf, genic_af, intergenic_af),
               names_to = c("region_type", "genome"),
               names_pattern = "(genic|intergenic)_(cf|gf|af)",
               values_to = "count") %>%
  mutate(
    genome = recode(genome, cf = "Canfam4", gf = "Gray fox", af = "Arctic fox"),
    region_type = recode(region_type, genic = "Genic", intergenic = "Intergenic"))

gf_af_fst_summary <- cf_gf_af_fst_long %>% 
  filter(!outlier_af=="match_rest") %>% 
  filter(!genome=="Canfam4") %>% 
  select(outlier_af,genome,region_type,count)%>% 
  mutate(outlier_af = recode(outlier_af, `af_outlier` = "heterospecific only", `gf_outlier` = "grayfox only",`match_outlier` = "both")) %>% 
  dplyr::rename("outlier"="outlier_af")

gf_cf4_fst_summary <- cf_gf_af_fst_long %>% 
  filter(!outlier_af=="match_rest") %>% 
  filter(!genome=="Arctic fox") %>% 
  select(outlier_cf4,genome,region_type,count)%>% 
  mutate(outlier_cf4 = recode(outlier_cf4, `cf4_outlier` = "heterospecific only", `gf_outlier` = "grayfox only",`match_outlier` = "both")) %>% 
  dplyr::rename("outlier"="outlier_cf4")


fst_region_summary <- rbind(gf_af_fst_summary, gf_cf4_fst_summary)

summary_table <- fst_region_summary %>% 
  filter(!outlier=="match_rest") %>% 
  group_by(region_type,genome,outlier) %>% 
  summarise(count=sum(count),window=n()) %>% mutate(prop=count/window)




ggplot(summary_table)+
  geom_col(aes(x=outlier,y=prop),position="dodge")+
  facet_wrap(~genome+region_type,scales="free",nrow=3)

ggplot(summary_table)+
  geom_col(aes(x=genome,y=prop,fill=outlier),position="dodge")+
  facet_wrap(~region_type,scales="free")+
  scale_fill_manual(values=c("gray20","orange1","darkorchid"))+
  theme_classic()+
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.key.size = unit(0.8,"line"),
        legend.text = element_text(size=10),
        legend.background = element_rect(fill="gray96"),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA,color="black",linewidth=0.5),
        axis.text=element_text(color="black",size=10))+
  ylab("Count per Window")+xlab("Annotation")



comp_gf_af <- cf_gf_af_fst_long %>% 
  filter(!outlier_af=="match_rest") %>% 
  filter(!genome=="Canfam4") %>% 
  group_by(region_type,genome,outlier_af) %>% 
  summarise(count=sum(count),window=n()) %>% mutate(prop=count/window)





