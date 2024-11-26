Grayfox
================

-   [Reference genome](#reference-genome)
-   [Samples](#samples)
    -   [PCA](#pca)
    -   [Map (Fig1A)](#map-fig1a)
-   [Demographies](#demographies)
    -   [Masks for smc++](#masks-for-smc)
    -   [Input for smc++](#input-for-smc)
    -   [Run smc++ estimate and plot](#run-smc-estimate-and-plot)
    -   [Plot smc++ output](#plot-smc-output)
    -   [MSMC2](#msmc2)
    -   [Plot msmc2 output](#plot-msmc2-output)
    -   [Run Stairway Plot 2](#run-stairway-plot-2)
    -   [Plot MSMC2 Stairway Plot 2](#plot-msmc2-stairway-plot-2)

## Reference genome

> Armstrong et al. (2024). Chromosome-level assembly of the gray fox
> (Urocyon cinereoargenteus) confirms the basal loss of PRDM9 in
> Canidae, G3 Genes\|Genomes\|Genetics, 14:4jkae034,
> <https://doi.org/10.1093/g3journal/jkae034>

<img src="Grayfox_files/input/circos_plot.jpeg" width="400"/>

## Samples

Whole genome sequencing of gray fox

Data from
[ncbi.nlm.nih.gov/bioproject/PRJNA966176/](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA966176/)

> Preckler-Quisquater et al. (2023). Can demographic histories explain
> long-term isolation and recent pulses of asymmetric gene flow between
> highly divergent grey fox lineages? Molecular Ecology, 32, 5323–5337.
> <https://doi.org/10.1111/mec.17105>

Filter VCF for minor allele frequency = 0.05 and remove missing data

``` bash
vcftools --gzvcf grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz --maf 0.05 --max-missing 1 --recode --out grayfox_n41
```

### PCA

Conduct principal component analysis on resulting VCF

``` r
library("SNPRelate")

# Read in VCF with all 41 samples
vcf.fn <-"grayfox_n41.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)

pca <- ccm_pca$eigenval[1:20]
pca_perc <- pca/sum(pca)*100
pve <- data.frame(PC = 1:20, pve = pca/sum(pca)*100)

tab <- data.frame(sample.id = ccm_pca$sample.id,
                  EV1 = ccm_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = ccm_pca$eigenvect[,2],    # the second eigenvector
                  EV3 = ccm_pca$eigenvect[,3],
                  EV4 = ccm_pca$eigenvect[,4],
                  EV5 = ccm_pca$eigenvect[,5],
                  stringsAsFactors = FALSE)

write_tsv(tab,"input/pca_pc1-5.txt")
```

Plot PCA

``` r
gfpop <- read_tsv("Grayfox_files/input/mainland_WGSmetadata.txt") %>% select(SRA,Region,State)

pca_tab <- read_tsv("Grayfox_files/input/pca_pc1-5.txt") %>% left_join(gfpop,by=c("sample.id"="SRA")) %>% 
  mutate(id=as.numeric(str_remove(sample.id, "SRR24465"))) %>% 
  mutate(label=case_when(id %in% c(296,287,288,284,292,290) ~"west",
                         id %in% c(286,275,272,271,305,304) ~"east",
                         TRUE ~ "drop"))

ggplot()+
  geom_point(data=subset(pca_tab,label=="drop"),aes(x=EV2,y=EV1,color=label),size=2)+
  geom_point(data=subset(pca_tab,label!="drop"),aes(x=EV2,y=EV1,color=label),size=2)+
  theme_classic()+ 
  xlab("PC2") + ylab("PC1")+
  theme(legend.position=c(0.15,0.2),
        legend.title = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.spacing.y = unit(0, "cm"),
        legend.box.background = element_rect(color="black",fill=NA),
        panel.border = element_rect(color="black",fill=NA),
        axis.line = element_blank(),
        axis.text = element_text(color="black"))+
  scale_color_manual(values=c("red","gray","black"))
```

![](Grayfox_files/figure-gfm/plotPCA-1.png)<!-- -->

Map of all samples

``` r
usa <- map_data("state")


pop <- read_tsv("Grayfox_files/input/mainland_WGSmetadata.txt") %>% 
  dplyr::rename("long"=Lon,"lat"=Lat) %>% 
  filter(Dataset=="Sacks") %>%   
  mutate(SRA=as.numeric(str_remove(SRA, "SRR24465"))) %>% 
  mutate(label=case_when(SRA %in% c(296,287,288,284,292,290) ~"west",
                         SRA %in% c(286,275,272,271,305,304) ~"east",
                         TRUE ~ "drop"))


# Custom labeling function for longitude
lon_labels <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "W", "E"))}

# Custom labeling function for latitude
lat_labels <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "S", "N"))}


ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
               fill = "white", color="gray70", linewidth=0.2) + 
  geom_point(data=pop, aes(x=long, y = lat, color=label,size=`Depth(x)`),stroke=0.25) +
    #geom_text(data=subset(pop), aes(x=long, y = lat, color=State,label=SRA),size=3) +
  scale_color_manual(values=c("red","gray70","black"))+
  coord_fixed(xlim = c(-123, -78),  ylim = c(25, 40)) +
  scale_x_continuous(labels = lon_labels,position = "top") +
  scale_y_continuous(labels = lat_labels,position = "right") +
  xlab("") + ylab("") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.tag.position = c(0.04,0.98),
        plot.margin = margin(0.5,0.1,0,-0.5, "lines"),
        plot.tag = element_text(size=13,face="bold"),
        legend.position = "none",
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black"),
        panel.background=element_rect(color="white",fill="white",linewidth = 0.5))
```

![](Grayfox_files/figure-gfm/map-1.png)<!-- -->

### Map (Fig1A)

``` r
ggplot() +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), 
               fill = "white", color="gray70", linewidth=0.2) + 
  geom_point(data=subset(pop,label!="drop"), aes(x=long, y = lat, color=label),stroke=0.25,size=2) +
  coord_fixed(xlim = c(-123, -72.8),  ylim = c(25, 43.8)) +
  scale_color_manual(values = c("gray70", "black"))+
  geom_point(aes(x=-72.20392,y=44.14799),shape=18,size=3,color="orange1")+  
  geom_point(aes(x=-79.5,y=30),shape=18,size=3, color="orange1")+ 
  geom_text(aes(x=-74.5,y=30),label="Reference", size=3.5)+
  geom_text(aes(x=-74.5,y=28.5),label="genome", size=3.5)+
  scale_x_continuous(labels = lon_labels) +
  scale_y_continuous(labels = lat_labels) +
  xlab("") + ylab("") +
  labs(tag="A)")+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        legend.position = c(0.18,0.17),
        legend.direction = "horizontal",
        plot.tag.position = c(0.04,0.98),
        plot.margin = margin(0.5,0.1,0,-0.5, "lines"),
        plot.tag = element_text(size=13,face="bold"),
        legend.text = element_text(size=10),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.spacing.y = unit(0, "cm"),
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black"),
        panel.background=element_rect(color="white",fill="white",linewidth = 0.5))+
  guides(color = guide_legend(title.position = "top",title.hjust = 0.5,title="WGS"))
```

![](Grayfox_files/figure-gfm/Fig1A-1.png)<!-- -->

## Demographies

### Masks for smc++

Genic regions

<details>
<summary>
Code for generating genic bed files
</summary>

``` r
##Arctic Fox genic regions ±1kb

arcanno <- fread("arcticFox.ncbiRefSeq.gtf.gz")

arcch <- read_tsv("arc_chroms.txt", col_names = c("chrom","name","scaff")) %>% select(chrom,scaff) 

arcanf <- arcanno %>% rename("V1"="#!annotation-source NCBI Vulpes lagopus Annotation Release 100") %>% 
  filter(str_detect(V3, "CDS")) %>% mutate(len=V5-V4) %>% separate(V9, sep=";", into=c("gene","transcript")) %>%
  left_join(arcch,by=c("V1"="scaff")) 

arcanch <- arcanf  %>% group_by(gene) %>% slice(which.max(len)) 

arc1kb <- arcanch %>%  ungroup() %>% 
  mutate(newstart=V4-1000, newend=V5+1000) %>% select(chrom,newstart,newend) %>% arrange(chrom,newstart)

arc1kb %>% write_tsv("arcticfox.gtfgenic.bed", col_names = F)


##Gray Fox genic regions ±1kb

gfanno <- fread("grayFox.gtf.gz")

gfchroms <- read_tsv("gf_chroms.txt", col_names = c("scaf","chrom")) %>% 
  separate(scaf, remove=F, c(NA,NA,NA,NA,NA,"length"))

gfanf <- gfanno %>% left_join(gfchroms,by=c("V1"="scaf")) %>% 
  filter(str_detect(V3, "CDS")) %>% mutate(len=V5-V4) %>% separate(V9, sep=";", into=c("t","g")) %>%  
  mutate(g = str_remove_all(g, "\"")) %>% separate(g, sep="d", into=c("id","gid")) %>%
  mutate(gid=as.numeric(str_remove(gid,"g"))) %>% arrange(gid,len)

gfanchr <- gfanf %>% group_by(gid) %>% filter(len==max(len)) %>% drop_na(chrom) %>% ungroup()

gf1kb <- gfanchr %>% mutate(newstart=V4-1000, newend=V5+1000) %>% 
  mutate(start=case_when(newstart<0 ~ 1, TRUE ~ newstart)) %>%
  mutate(end=case_when(newend>as.numeric(length) ~ as.numeric(length), TRUE ~ newend)) %>%
  dplyr::select(chrom,start,end)

gf1kb %>% write_tsv("grayfox.gtfgenic.bed", col_names = F)


##CanFam4 genic regions ±1kb

cfanno <- fread("canFam4.ncbiRefSeq.gtf.gz")

cfch <- read_tsv("cf4_chroms.txt",col_names = F) 

cfanf <- cfanno %>% 
  filter(str_detect(V3, "CDS")) %>% mutate(len=V5-V4) %>% separate(V9, sep=";", into=c("gene","transcript")) %>% left_join(cfch,by=c("V1"="X1")) 

cfanch <- cfanf  %>% drop_na(X2) %>% group_by(gene) %>% slice(which.max(len)) 

cf1kb <- cfanch %>%  ungroup() %>% 
  mutate(newstart=V4-1000, newend=V5+1000) %>% select(V1,newstart,newend) %>% arrange(V1,newstart)

cf1kb %>% write_tsv("canfam4.gtfgenic.bed", col_names = F)
```

</details>

<br>

Merge bed files of genic regions with bed files of low mappability and
repeat-masked regions to generate final mask files

mask for grayfox: **grayfox.rep.map.genic.bed.gz**

mask for arcticfox: **arcticfox.rep.map.genic.bed.gz**

mask for canfam4: **canfam4.rep.map.genic.bed.gz**

### Input for smc++

Convert VCF to smc++ input

``` bash
sbatch Grayfox_files/scripts/smc_vcf_gf.sh

sbatch Grayfox_files/scripts/smc_vcf_af.sh

sbatch Grayfox_files/scripts/smc_vcf_cf4.sh
```

### Run smc++ estimate and plot

combine input files to generate a composite likelihood estimate for
fitting a population size history

do six times, once for each combination of population and reference
genome

``` bash
sbatch Grayfox_files/scripts/smc_est_plot.sh
```

### Plot smc++ output

![](Grayfox_files/figure-gfm/plotsmc-1.png)<!-- -->

### MSMC2

Filter gVCF, make masks, split files, and run MSMC2

``` bash
sbatch Grayfox_files/scripts/gvcf_filt.sh

## do for each reference genome, grayfox example
sbatch Grayfox_files/scripts/gvcf_split.sh

## run generate_multihetsep.py, grayfox example
sbatch Grayfox_files/scripts/msmc_inp.sh

## run MSMC2, grayfox east example
sbatch Grayfox_files/scripts/msmc_size_e.sh
```

### Plot msmc2 output

### Run Stairway Plot 2

<details>
<summary>
Blueprint file for running Stairway Plot
</summary>

``` bash
```

</details>

<br>

### Plot MSMC2 Stairway Plot 2

![](Grayfox_files/figure-gfm/plotdem-1.png)<!-- -->
