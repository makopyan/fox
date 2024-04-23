library(data.table)
library(tidyverse)

files_dir <- "/project/jazlynmo_738/DataRepository/perChrom_recombination"

# Get list of files
east_file_names <- list.files(files_dir, pattern = "^gf_east_n14_\\d+\\.rmap$", full.names = TRUE)

west_file_names <- list.files(files_dir, pattern = "^gf_west_n12_\\d+\\.rmap$", full.names = TRUE)

# Read files and add a column with the file integer
east_data_list <- lapply(seq_along(east_file_names), function(i) {
  data <- fread(east_file_names[i],col.names = c("start","end","recomb"))
  data$chrom <- i
  return(data)
})

west_data_list <- lapply(seq_along(west_file_names), function(i) {
  data <- fread(west_file_names[i],col.names = c("start","end","recomb"))
  data$chrom <- i
  return(data)
})

# Combine data frames into a single data frame
reast <- bind_rows(east_data_list)

rwest <- bind_rows(west_data_list)

window_size <- 50000

east_rmap_df <- reast %>% 
  mutate(window_pos1 = floor(start / window_size) * window_size) %>%
  mutate(window_pos2 = floor(start / window_size) * 
           window_size + window_size) %>%
  group_by(chrom, window_pos1) %>%
  summarise(recomb_avg = mean(recomb, na.rm = TRUE)) %>%
  ungroup

west_rmap_df <- rwest %>% 
  mutate(window_pos1 = floor(start / window_size) * window_size) %>%
  mutate(window_pos2 = floor(start / window_size) * 
           window_size + window_size) %>%
  group_by(chrom, window_pos1) %>%
  summarise(recomb_avg = mean(recomb, na.rm = TRUE)) %>%
  ungroup

rmap_df <- bind_rows("east"=east_rmap_df,"west"=west_rmap_df,.id="pop")


#write_tsv(rmap_df,"rmap_50kb.txt")

ggplot(data=rmap_df, aes(x=window_pos1/1e6, y=recomb_avg*1e8, color=as.factor(chrom))) + 
  scale_color_manual(values = rep(c("gray60", "black"), 100)) +
  geom_smooth(method="loess",span=0.1, se=F, linewidth=0.7) +
  xlab("Position (mb)")+ylab("Recombination Rate (cM/Mb)")+
  facet_grid(cols = vars(chrom), rows=vars(pop),
             space = "free_x",  switch="x",
             scales = "free_x") +     
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background=element_blank(),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.text.x = element_blank(),
    axis.line.y = element_line(colour = "black", linewidth = 0.5),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = element_text(colour = "black", size = 12),
    legend.position="none",
    strip.text = element_text(colour = "black", size = 12),
    axis.title.y = element_text(colour = "black", size = 12),   
    panel.spacing = unit(0.0, "cm"),
    strip.background=element_rect(fill="gray90")) 

#ggsave("rmap_50kb.png", width=11, height = 4, units="in")

