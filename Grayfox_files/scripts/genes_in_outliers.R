library(tidyverse)
library(data.table)
library(GenomicRanges)

cf_gf_af_fst <- read_tsv("~/Downloads/cf4_gf_af_window_phyloP.tsv") %>% 
  dplyr::rename("gf_chr"="seqnames","gf_st"="start","gf_end"="end") %>% 
  select(gf_chr,gf_st,gf_end,af_chr,af_st,af_end,outlier_af,cf4_chr,cf4_st,cf4_end,outlier_cf4) %>% 
  mutate(gf_chr = as.character(gsub("^chr", "", gf_chr)), af_chr = as.character(gsub("^chr", "", af_chr)),cf4_chr = as.character(gsub("^chr", "", cf4_chr))) 

af_gf <- cf_gf_af_fst %>% filter(outlier_af=="af_outlier") %>% 
  select(af_chr,af_st,af_end) %>% mutate(af_chr=paste0("chr",af_chr))

cf_gf <- cf_gf_af_fst %>% filter(outlier_cf4=="cf4_outlier") %>% 
  select(cf4_chr,cf4_st,cf4_end) %>% mutate(cf4_chr=paste0("chr",cf4_chr))


cfanf <- fread("~/Desktop/Foxes/DATA/canfam4.chr.gtf.gz",
               col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))



afanf <- fread("~/Desktop/Foxes/DATA/arcticfox.chr.gtf.gz",
               col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))



# Convert annotation data frame to GRanges object
annotation_gr <- GRanges(
  seqnames = afanf$seqname,
  ranges = IRanges(start = afanf$start, end = afanf$end),
  strand = afanf$strand,
  mcols = DataFrame(
    feature = afanf$feature,
    attribute = afanf$attribute))

# Convert window positions to GRanges object
windows_gr <- GRanges(
  seqnames = af_gf$af_chr,
  ranges = IRanges(start = af_gf$af_st, end = af_gf$af_end))

# Make seqlevels compatible by setting the same style
seqlevelsStyle(annotation_gr) <- seqlevelsStyle(windows_gr)

# Find overlaps between annotations and windows
overlaps <- findOverlaps(annotation_gr, windows_gr)

# Extract features that fall within windows
features_in_windows <- annotation_gr[queryHits(overlaps)]

# Filter for CDS features
cds_in_windows <- features_in_windows[features_in_windows$mcols.feature == "CDS"]

# Extract gene IDs from the attribute field
gene_ids <- gsub('.*gene_id ""([^"]*)"".*', "\\1", cds_in_windows$mcols.attribute)

genes_af<-as.data.frame(gene_ids) %>% filter(!str_detect(gene_ids,"LOC")) %>% distinct(gene_ids)

write_tsv(genes_af,"~/Downloads/af_outlier_genes.txt",col_names = F)

