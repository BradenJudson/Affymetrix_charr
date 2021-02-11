#Braden Judson
#February 2020
#Formatting and filtering SNP chip data from cave charr samples
#Generalized and personalized version of @PhDMattyB's Script: https://github.com/PhDMattyB/ped_file_formating

setwd("C:/Users/Brade/OneDrive/Affymetrix_charr/Affymetrix_charr")


#### 1. LIBRARIES ####

library(tidyverse)
library(dplyr)
library(readr)

mytheme <- theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#### 2a. MAP FILE ANNOTATION ####


MAP <- read_tsv('cavecharr_genotypes.map',                                      # Read in .map file from affymetrix.
                col_names = TRUE)                                               # Retain column names.


affy_snp_id <- read_tsv('array_markers.tsv', col_names = TRUE) %>%              # Read in array marker file. 
                        rename(probeset_id = Name)                              # Rename probeset identifiers. 


affy_ID <- affy_snp_id[affy_snp_id$probeset_id %in% MAP$`Marker ID`, ] %>%      # Filter by IDs in map file.
  rename('Marker ID' = probeset_id) %>%                                         # Rename probeset IDs. 
  select(-c(Chromosome, Position))                                              # Remove unnecessary columns. 

vcf <- read_tsv("array_markers.NCBI.Salp.V2.vcf",                               # Read in variant call file.
                skip = 3, col_names = TRUE)                                     # Skip first 3 rows, retain column names.


vcf_loci <- vcf[vcf$ID %in% affy_ID$`Affx ID`,] %>%                             # Retain loci from affymetrix marker set.
  select("#CHROM", "ID", "POS") %>%                                             # Select chromosome and loci identifiers. 
  rename("Affx ID" = ID)                                                        # Rename ID to be more specific. 

geno <- merge(vcf_loci, affy_ID) %>%                                            # Merge info from vcf and affymetrix file.
  merge(MAP)  %>%                                                               # Merge info from map.
  rename("#Chromosome" = ChromNum)                                              # Rename. Important to keep # symbol so headers aren't incorporated into .ped or .map files. 

sample(unique(geno$`#Chromosome`), 20)                                          # Chromosomes are identified numerically.

# Salvelinus alpinus has 37 autosomal chromosomes. 
# AC6 and AC4 are split into two independent acrocentric chromosomes (i.e., serves as 39 linkage groups) (work done by Christine Ouellet).
# All other loci not placed on a LG scaffold are placed on AC40 pseudochromosome. 
geno = mutate(.data = geno,
              LG = as.factor(case_when(
                `#Chromosome` == '1' ~ 'AC01',
                `#Chromosome` == '2' ~ 'AC02',
                `#Chromosome` == '3' ~ 'AC03',
                `#Chromosome` == '4' ~ 'AC04p',
                `#Chromosome` == '5' ~ 'AC04q.1:29',
                `#Chromosome` == '6' ~ 'AC04q.2',
                `#Chromosome` == '7' ~ 'AC05',
                `#Chromosome` == '8' ~ 'AC06',
                `#Chromosome` == '9' ~ 'AC06',
                `#Chromosome` == '10' ~ 'AC07',
                `#Chromosome` == '11' ~ 'AC08',
                `#Chromosome` == '12' ~ 'AC09',
                `#Chromosome` == '13' ~ 'AC10',
                `#Chromosome` == '14' ~ 'AC11',
                `#Chromosome` == '15' ~ 'AC12',
                `#Chromosome` == '16' ~ 'AC13',
                `#Chromosome` == '17' ~ 'AC14',
                `#Chromosome` == '18' ~ 'AC15',
                `#Chromosome` == '19' ~ 'AC16',
                `#Chromosome` == '20' ~ 'AC17',
                `#Chromosome` == '21' ~ 'AC18',
                `#Chromosome` == '22' ~ 'AC19',
                `#Chromosome` == '23' ~ 'AC20',
                `#Chromosome` == '24' ~ 'AC21',
                `#Chromosome` == '25' ~ 'AC22',
                `#Chromosome` == '26' ~ 'AC23',
                `#Chromosome` == '27' ~ 'AC24',
                `#Chromosome` == '28' ~ 'AC25',
                `#Chromosome` == '29' ~ 'AC26',
                `#Chromosome` == '30' ~ 'AC27',
                `#Chromosome` == '31' ~ 'AC28',
                `#Chromosome` == '32' ~ 'AC30',
                `#Chromosome` == '33' ~ 'AC31',
                `#Chromosome` == '34' ~ 'AC32',
                `#Chromosome` == '35' ~ 'AC33',
                `#Chromosome` == '36' ~ 'AC34',
                `#Chromosome` == '37' ~ 'AC35',
                `#Chromosome` == '38' ~ 'AC36',
                `#Chromosome` == '39' ~ 'AC37',
                `#Chromosome` > '39' ~ 'AC40')))    # Pseudochromosome contig. 

geno$LG[is.na(geno$LG)] = 'AC40'                    # Place unanchored loci on contigs group.
length(unique(geno$LG))                             # 39 levels. 
class(geno$`Physical position`)                     # numeric.

# Split AC4.1q:29. Christine Ouellet's work suggests split ~ 40,000,000bp. 

AC04q.1_29 <- geno %>%                              # Isolate AC04.1q:29.
  filter(LG == 'AC04q.1:29') %>%  
  select(-c(LG)) %>%                                # Makes next step easier.
  mutate(LG = as.factor(case_when(                  # Split above/below 40,000,000bp. 
     POS > 40000000 ~ 'AC04q.1',                    # Above is AC04q.1
     POS < 40000000 ~ 'AC29')))                     # Below is AC29. 


genobind <- bind_rows(geno[geno$LG != "AC04q.1:29",],   # Combine to original. 
                  AC04q.1_29) %>% 
  rename(probeset_id = 'Marker ID')                     # Helpful downstream. 

dim(genobind)                                           # 65,235 rows by 24 columns. 
 

#### 1b. MAP FILE FORMATTING ####

#Read in the summary file from affymetrix and filter out the non-polymorphic snps.
polySNP = read_tsv('other_data.txt', col_names = T) %>% 
  filter(ConversionType %in% c('PolyHighResolution', 'NoMinorHom'))



#Join the map and the polymorphic snp summary file by the matching name so that the map file only retains polymorphic SNPs.
poly_snps_only = left_join(polySNP, genobind, by = 'probeset_id') %>%           # Retain polymorphic loci.
  select('LG', probeset_id,                                                     # Select most relevant columns.
         `Genetic distance`, `POS`) %>% 
  rename(`Position` = `POS`)                                                    # Rename the columns to make downstream processing easier.

poly_snps_only <- na.omit(poly_snps_only)                                       # Remove loci not anchored to a physical position in the genome (n = 183).
sum(is.na(poly_snps_only$`Position`))                                           # No missing data.

unique(poly_snps_only$LG)                                                       # 41 levels. 

PolyAC04_adj <- poly_snps_only[poly_snps_only$LG %in% 'AC04q.1',] %>%           # Correct position change from moving block to AC29. 
  mutate(posadj = (Position - 40000000)) %>% select(-c(Position)) %>% 
  rename("Position" = "posadj")

poly_snps_only <- rbind(poly_snps_only[poly_snps_only$LG != 'AC04q.1',],       # Bind corrected positions to original dataframe.
                        PolyAC04_adj)

write_tsv(x = poly_snps_only, col_names = FALSE, 'charrSNP.map')                # Write to harddrive. Can then use in PGDSPIDER, PLINK, etc. 

SNP_perLG <- poly_snps_only %>%                                                 # SNPs per linkage group. 
  group_by(LG) %>%                                                              # Average of 97 SNPs per LG.
  summarise(SNPs = length(LG))                                                  # 

sum(SNP_perLG[SNP_perLG$LG != 'AC40', 2]) / sum(SNP_perLG$SNPs) * 100           # 61.7% of SNPS anchored to significant linkage group.
sum(SNP_perLG[SNP_perLG$LG == 'AC40', 2]) / sum(SNP_perLG$SNPs) * 100           # 38.2% of SNPs on pseudochromosome contig. 


ggplot(poly_snps_only[poly_snps_only$LG != 'AC40',], aes(x = Position)) +
  geom_histogram() + ylab("Frequency") + theme_bw()                             # Right skewed. 

#### 2. PED FILE FORMATTING ####

dim(poly_snps_only)                                        # 3870 x 4. 

origmap <- read_tsv("cc_genotypes.map",                    # Original MAP file.
                    col_names = TRUE, skip = 3)            # skip 3 rows, retain col names.
origped <- read_tsv("genotype_ped_info.txt",
                    col_names = TRUE)
loci_head <- c("Sample", origmap$`Marker ID`)              # 66262 loci.


PED <- read_tsv("cave_charr_genotype.ped",                 # Read in. Takes a minute, pretty big file. Skip first row.
                col_names = loci_head, skip = 1)           # important that the order of header is the same as order in PED file.
rownames(PED) <- PED$Sample                                # set row names as 'Sample' column is discarded later.

dim(PED)                                                   # 19 x 66263 (ie., 19 individuals across >66k loci).

polyPED <- PED[colnames(PED) %in% poly_snps_only$probeset_id] # Retain polymorphic loci only.
rownames(polyPED) <- PED$Sample                               # Set rownames to sample IDs. 
polyPED <- rownames_to_column(polyPED) %>% 
  rename(Sample = rowname)

colnames(origped)[7] <- 'Sample'                           # Rename for easier joining.
ccPED <- left_join(origped, polyPED, by = 'Sample')        # left join by common column.
ccPED <- ccPED[,-7]                                        # Remove joining column - not needed.
write_tsv(x = ccPED,                                       # Write file as PED. 
          col_names = FALSE,                               
          'cavecharr_polySNP.ped')                         # Format acceptable for PLINK. 


##### Plot distribution of loci #####

assembly <- read.csv("LGstats.csv", header = TRUE)                              # Read in contig lengths, summary stats.

polySNPs <- merge(poly_snps_only[poly_snps_only$LG != 'AC40',],                 # Discard 'contigs' (i.e. AC40)
                  assembly[,c(5,11)], by = "LG")                                # add coding variable for each chromosome.

chroms <- ggplot() +                                                            # Plot chromosome sizes.
  geom_segment(data = assembly,                                                 
               aes(x = num, xend = num, y = 0, yend = as.numeric(Length)),      # Lengths for each linkage group, rounded ends.
                   lineend = "round", color = "grey", size = 3) +
  scale_x_continuous(breaks = 1:nrow(assem4),                                   # Properly labels LGs based on coding variable.
                   labels = c(assem4$LG)) + mytheme                             # custom theme.

plot(chroms)                                                                    # Plot LGs.

loci <- chroms +                                                                # Plot loci on LGs.
  geom_segment(data = polySNPs,                                                 # Only use polymorphic loci on known LGs.
               aes(x = num-0.1, xend = num+0.1, 
                   y = Position, yend = Position), size = 0.1) 

plot(loci) + mytheme +                                                          # Tidy plot.     
  xlab("") + ylab("Position (bp)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Note that the above plot is simply polymorphic loci. 
     # Upon further filtering, redo above visualization. 

#### 3. PLINK FORMATTING ####

#The following steps are variable and depend on the type of analysis being conducted. 

#### 3a. REMOVE LOCI IN LDE ####
snps_in_LD <- read_tsv('plink.prune.out', col_names = F) #gives snps that are in LDE and are removed. 
dim(snps_in_LD) #34 loci in LDE. 






