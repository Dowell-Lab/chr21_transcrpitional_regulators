library(tidyverse)
library(ggplot2)
library(ggridges)
library(ComplexHeatmap)
library(dplyr)
library(matrixStats)



indir = "/Shares/down/public/HTP/RNAseq/"
whichoutput="T21vsD21_noDNAdosagecorrection"
comorbidfileroot = "Patient_mondo_logical.csv"
comorbidfile=paste0("/Shares/down/public/HTP/RNAseq/outputdata/", comorbidfileroot, sep="")
allgenesfile = paste(indir, "outputdata/", whichoutput, "/allgeneswithgenenames_res_not_unique.csv", sep="")
normcountsfile =paste(indir, "outputdata/",whichoutput,"/normcounts.csv", sep="")
outdir="/scratch/Shares/dowell/temp/ChrisO/t21_reviewarticle/ouput/T21vsD21_noDNAdosagecorrection/furtherprocessing/"
min_n_gene_to_high=10
min_count_sum=100
min_people_with_comorbid=10

#read in who has which comorbidities, get "complete_trisomy_21" vs "mosaic" vs etc...
comorbidf <-  read.csv(comorbidfile, row.names = 1)
comorbidf
colnames(comorbidf) <- paste('comorbid', colnames(comorbidf), sep = '_')
comorbidf[] <- lapply(comorbidf, as.logical)
comorbidf$Patient <- gsub("-",".",rownames(comorbidf))

complete_trisomy_21 <- comorbidf %>% dplyr::filter(comorbid_complete_trisomy_21==1)
dim(complete_trisomy_21)
mosaic_trisomy_21 <-  comorbidf %>% dplyr::filter(comorbid_mosaic_trisomy_21==1)
dim(mosaic_trisomy_21)
translocation_Down_syndrome <-  comorbidf %>% dplyr::filter(comorbid_translocation_Down_syndrome==1)
dim(translocation_Down_syndrome)
mosaic_translocation_Down_syndrome <-  comorbidf %>% dplyr::filter(comorbid_mosaic_translocation_Down_syndrome==1)
dim(mosaic_translocation_Down_syndrome)

#read in genes and their chromosomes
geneinfo = "/Shares/down/public/HTP/RNAseq/selfannotated/genes.csv"
genedf = read.csv(geneinfo)
genedf$ENSEMBL = genedf$gene_id 
genedf$SYMBOL = genedf$gene_name 
genedf$CHR = genedf$seqnames 

genesymbolchr <- genedf %>% dplyr::select(ENSEMBL, CHR, SYMBOL) %>% unique()
genesymbolchr_just21 <- genesymbolchr %>% filter(CHR==21)
genesymbolchr_just22 <- genesymbolchr %>% filter(CHR==22)
#ensb gene to symbol not 1 to 1, ??? CDO looks like it is to me...
dim(genesymbolchr)
length(unique(genesymbolchr$ENSEMBL))
dim(genesymbolchr_just21)
length(unique(genesymbolchr_just21$ENSEMBL))


#read in which patients are at all T21 
mani <- read.csv(paste(indir, "inputdata/manifest_20230829_152940.csv", sep=""))
manimini <- mani %>% dplyr::select(Patient, sample_type)

manimini <- manimini %>%
  dplyr::mutate(Patient = Patient %>% 
                  str_replace_all("-", "."))

dim(manimini)
head(manimini)

maniminiD21 <- manimini %>% filter(sample_type=="D21")
maniminiT21 <- manimini %>% filter(sample_type=="T21")


complete_trisomy_21_with_RNA <- complete_trisomy_21$Patient[complete_trisomy_21$Patient %in% maniminiT21$Patient]
mosaic_trisomy_21_with_RNA <- mosaic_trisomy_21$Patient[mosaic_trisomy_21$Patient %in% maniminiT21$Patient]
translocation_Down_syndrome_with_RNA <- translocation_Down_syndrome$Patient[translocation_Down_syndrome$Patient %in% maniminiT21$Patient]
mosaic_translocation_Down_syndrome_with_RNA <- mosaic_translocation_Down_syndrome$Patient[mosaic_translocation_Down_syndrome$Patient %in% maniminiT21$Patient]


#read in counts per gene and make a D21/T21 count data set
ncdf <- read.csv(normcountsfile)
head(ncdf)
rownames(ncdf)<-ncdf$X
ncdf <- ncdf %>% dplyr::select(-X)
ncdfchr21<- ncdf %>% dplyr::filter(rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL)
ncdfchr21_D21 <- ncdf %>% dplyr::filter(rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL) %>% dplyr::select(maniminiD21$Patient)
ncdfnotchr21_D21 <- ncdf %>% dplyr::filter(!rownames(ncdf) %in% genesymbolchr_just21$ENSEMBL) %>% dplyr::select(maniminiD21$Patient)
ncdf_T21<-ncdf %>% dplyr::select(complete_trisomy_21_with_RNA) #254 people
ncdf_D21<-ncdf %>% dplyr::select(maniminiD21$Patient) #96 people
ncdfchr21_T21sim <- ncdfchr21_D21*1.5
ncdf_T21sim <- rbind(ncdfchr21_T21sim, ncdfnotchr21_D21)

###CDO ridge plot code###
#ncdf_T21 and ncdf_D21 which has 96 D21 individuals and 254 T21 individuals (excluding the mosaics/translocation patients)

#need to get gene name back, and only table with the 34 genes we want

ncdf_T21$ENSEMBL <- rownames(ncdf_T21) #make a new column for ENSEMBL names
ncdf_T21 <- ncdf_T21 %>% select(ENSEMBL, everything()) #move it to first column

ncdf_T21_genename <- left_join(genesymbolchr_just21, ncdf_T21, by = "ENSEMBL") ##H2BC12L has the symbol H2BFS
subset_ncdf_T21 <- ncdf_T21_genename %>% filter(SYMBOL %in% c("GABPA", "RUNX1",
                                                              "AIRE", "PKNOX1",
                                                              "ZBTB21", "PRDM15",
                                                              "ETS2", "ERG", "SIM2", 
                                                              "OLIG1", "OLIG2", 
                                                              "BACH1", "TCP10L", 
                                                              "PAXBP1", "PCBP3",
                                                              "NRIP1", "RIPPLY3",
                                                              "BTG3", "PTTG1IP", "SCAF4",
                                                              "PRMT2", "DNMT3L",
                                                              "HMGN1", "SETD4",
                                                              "N6AMT1", "MORC3",
                                                              "BRWD1", "CHAF1B",
                                                              "H2BFS", #new genes, H2BFS is H2BC12L
                                                              "U2AF1", "SON", "RBM11",
                                                              "SIK1", "DYRK1A", "SOD1", "RRP1B", "HSF2BP" ))
#remove first two columns
subset_ncdf_T21 <- subset_ncdf_T21[, -c(1:2)]
# transpose so that we can plot it
long_normcounts_T21 <- pivot_longer(subset_ncdf_T21, cols = -SYMBOL, names_to = "Patient", values_to = "Norm_counts")
# add ploidy column
long_normcounts_T21_ploidy <- long_normcounts_T21 %>% mutate(ploidy = "T21")

ggplot(long_normcounts_T21, aes(x = Norm_counts, y = SYMBOL, fill = SYMBOL)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
##works! But I think I need to Zscore so we can see it better##
##should be able to fill by T21 vs. D21, Just stick the other df on the bottom?? would that work...
#zscore T21
rownames(subset_ncdf_T21)<-subset_ncdf_T21$SYMBOL# need to get rid of first column again
subset_ncdf_T21 <- subset_ncdf_T21 %>% dplyr::select(-SYMBOL)

rowSigmaT21 <- apply(subset_ncdf_T21, 1, sd, na.rm = TRUE) 
rowMuT21 <- rowMeans(subset_ncdf_T21, na.rm = TRUE)
#zscore_subset_ncdf_T21 <- (subset_ncdf_T21 - rowMuT21) / rowSigmaT21
# I think I need to Zscore from the Disomics. Or make a giant DF with counts for both disomics and trisomics...
# Transform the d21 df, join the two, remove the excess rows, use those means and std deviations to zscore the two smaller dfs

ncdf_D21$ENSEMBL <- rownames(ncdf_D21) #make a new column for ENSEMBL names
ncdf_D21 <- ncdf_D21 %>% select(ENSEMBL, everything()) #move it to first column

ncdf_D21_genename <- left_join(genesymbolchr_just21, ncdf_D21, by = "ENSEMBL") ##H2BC12L has the symbol H2BFS
subset_ncdf_D21 <- ncdf_D21_genename %>% filter(SYMBOL %in% c("GABPA", "RUNX1",
                                                              "AIRE", "PKNOX1",
                                                              "ZBTB21", "PRDM15",
                                                              "ETS2", "ERG", "SIM2", 
                                                              "OLIG1", "OLIG2", 
                                                              "BACH1", "TCP10L", 
                                                              "PAXBP1", "PCBP3",
                                                              "NRIP1", "RIPPLY3",
                                                              "BTG3", "PTTG1IP", "SCAF4",
                                                              "PRMT2", "DNMT3L",
                                                              "HMGN1", "SETD4",
                                                              "N6AMT1", "MORC3",
                                                              "BRWD1", "CHAF1B",
                                                              "SOD1", "RRP1B",
                                                              "H2BFS", #new genes, H2BFS is H2BC12L
                                                              "U2AF1", "SON", "RBM11",
                                                              "SIK1", "DYRK1A", "HSF2BP"))
#remove first two columns
subset_ncdf_D21 <- subset_ncdf_D21[, -c(1:2)]
# transpose so that we can plot it
long_normcounts_D21 <- pivot_longer(subset_ncdf_D21, cols = -SYMBOL, names_to = "Patient", values_to = "Norm_counts")
# add ploidy column
long_normcounts_D21_ploidy <- long_normcounts_D21 %>% mutate(ploidy = "D21")

rownames(subset_ncdf_D21)<-subset_ncdf_D21$SYMBOL #probably good to get those rownames set for when I inevitably have to delete the column again

#join to make the big df for averages and error to zscore smaller ones
subset_ncdf_T21$SYMBOL <- rownames(subset_ncdf_T21)
subset_ncdf_T21 <- subset_ncdf_T21 %>% select(SYMBOL, everything())

subset_ncdf_SYMBOL_D21_T21 <- left_join(subset_ncdf_T21, subset_ncdf_D21, by = "SYMBOL")
#remove label column
rownames(subset_ncdf_SYMBOL_D21_T21) <- subset_ncdf_SYMBOL_D21_T21$SYMBOL #column names again
subset_ncdf_SYMBOL_D21_T21 <- subset_ncdf_SYMBOL_D21_T21[, -c(1)]

#now get zscore values from large df
rowSigmaALL <- apply(subset_ncdf_SYMBOL_D21_T21, 1, sd, na.rm = TRUE) 
rowMuALL <- rowMeans(subset_ncdf_SYMBOL_D21_T21, na.rm = TRUE)

#Zscore D21 and T21 subsets
subset_ncdf_T21 <- subset_ncdf_T21[, -c(1)]
zscore_subset_ncdf_T21 <- (subset_ncdf_T21 - rowMuALL) / rowSigmaALL
subset_ncdf_D21 <- subset_ncdf_D21[, -c(1)]
zscore_subset_ncdf_D21 <- (subset_ncdf_D21 - rowMuALL) / rowSigmaALL

zscore_subset_ncdf_D21$SYMBOL <- rownames(zscore_subset_ncdf_D21) #get back SYMBOL column, by default it goes to the end of the df
zscore_subset_ncdf_T21$SYMBOL <- rownames(zscore_subset_ncdf_T21)

# transpose so that we can plot it
long_Zscore_D21 <- pivot_longer(zscore_subset_ncdf_D21, cols = -SYMBOL, names_to = "Patient", values_to = "zscore")
long_Zscore_T21 <- pivot_longer(zscore_subset_ncdf_T21, cols = -SYMBOL, names_to = "Patient", values_to = "zscore")

# add ploidy column
long_Zscore_D21_ploidy <- long_Zscore_D21 %>% mutate(ploidy = "D21")
long_Zscore_T21_ploidy <- long_Zscore_T21 %>% mutate(ploidy = "T21")
# stick them together
long_Zscore_D21_T21_ploidy <- rbind(long_Zscore_D21_ploidy, long_Zscore_T21_ploidy)


#plot
ggplot(long_Zscore_D21_T21_ploidy, aes(x = zscore, y = SYMBOL, fill = ploidy)) +
  geom_density_ridges(alpha = 0.6, scale = 0.99) +
  theme_ridges() + 
  theme(legend.position = "none") +
  xlim(c(-2.6, 2.6)) +
  scale_fill_manual(values = alpha(c("red1", "blue1"), 0.6)) +
  theme(axis.text.y = element_text(size = 9, color = "black", face = "bold"))  # Customize y-axis text appearance

##
# Get TPM average values for each gene in D21 and T21
head(subset_ncdf_D21)
head(subset_ncdf_T21)

subset_ncdf_D21_avg <- as_tibble(subset_ncdf_D21, rownames = "SYMBOL")
subset_ncdf_D21_avg <- subset_ncdf_D21_avg %>%
  mutate(average = rowMeans(select(., where(is.numeric)), na.rm = TRUE))
subset_ncdf_T21_avg <- as_tibble(subset_ncdf_T21, rownames = "SYMBOL")
subset_ncdf_T21_avg <- subset_ncdf_T21_avg %>%
  mutate(average = rowMeans(select(., where(is.numeric)), na.rm = TRUE))
subset_ncdf_D21_red <- subset_ncdf_D21_avg %>% select(1, ncol(subset_ncdf_D21_avg)) #remove other columns
subset_ncdf_T21_red <- subset_ncdf_T21_avg %>% select(1, ncol(subset_ncdf_T21_avg)) #remove other columns

result <- long_Zscore_D21_ploidy %>%
  left_join(subset_ncdf_D21_red, by = "SYMBOL") #Left join by the SYMBOL column
result2 <- long_Zscore_T21_ploidy %>%
  left_join(subset_ncdf_T21_red, by = "SYMBOL") #Left join by the SYMBOL column

long_Zscore_D21_T21_ploidy_TpmAvg <- rbind(result, result2) #stick them together (unsorted)
#get a sorted data frame also
result_sort <- result %>%
  arrange(average) # used arrange(desc()) the first time but it was the opposite order on the plot
result2_sort <- result2 %>%
  arrange(average)
long_Zscore_D21_T21_ploidy_TpmAvg_sort <- rbind(result_sort, result2_sort) #stick them together (sorted)


##PLOT with TPMS
#plot with sorted TPMs - use text_data_sort and long_Zscore_D21_T21_ploidy_TpmAvg_sort...
#otherwise, use text_data and long_Zscore_D21_T21_ploidy_TpmAvg

long_Zscore_D21_T21_ploidy_TpmAvg_sort_factor <- long_Zscore_D21_T21_ploidy_TpmAvg_sort %>%
  mutate(SYMBOL = factor(SYMBOL, levels = unique(SYMBOL))) #changed to a factor to preserve order of the dataframe for plotting
structure(long_Zscore_D21_T21_ploidy_TpmAvg_sort_factor)

# Prepare data for text annotations, including average by SYMBOL and ploidy dataframe --> so you can plot the TPM values in the second layer, and they are still linked to the gene (SYMBOL) and ploidy
text_data_sort <- long_Zscore_D21_T21_ploidy_TpmAvg_sort_factor %>%
  group_by(SYMBOL, ploidy) %>%
  summarise(average = mean(average, na.rm = TRUE), .groups = 'drop')

# Create the base plot
p <- ggplot(long_Zscore_D21_T21_ploidy_TpmAvg_sort_factor, aes(x = zscore, y = SYMBOL, fill = ploidy)) +
  geom_density_ridges(alpha = 0.6, scale = 0.99) +
  theme_ridges() + 
  theme(legend.position = "none") +
  xlim(c(-2.6, 2.6)) +
  scale_fill_manual(values = alpha(c("red1", "blue1"), 0.6)) +
  theme(axis.text.y = element_text(size = 9, color = "black", face = "bold")) +
  ylab("Genes") +
  theme(axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
  xlab("Z-scored normalized counts") +
  theme(axis.title.x = element_text(size = 10, color = "black", face = "bold"))
  

# Add average values for D21
p + geom_text(data = text_data_sort %>% filter(ploidy == "D21"), 
              aes(x = -2.3, y = SYMBOL, label = sprintf("D21 %.1f", average)), 
              family = "Fira Sans", size = 2.0, hjust = 0.7, vjust = -0.9, color = "deepskyblue4", fontface = "bold") +
  geom_text(data = text_data_sort %>% filter(ploidy == "T21"), 
            aes(x = -2.3, y = SYMBOL, label = sprintf("T21 %.1f", average)), 
            family = "Fira Sans", size = 2.0, hjust = 0.7, vjust = 1.5, color = "deepskyblue4", fontface = "bold") +
  labs(title = "Avg TPM") +   # Add plot title
  theme(plot.title = element_text(size = 8, face = "bold", color = "deepskyblue4", family = "Fira Sans"))


## heatmap code

heatmap_data <- text_data_sort %>%
  pivot_wider(names_from = ploidy, values_from = average, values_fill = NA) %>%
  arrange(SYMBOL)


# Base ridgeline plot
p_ridge <- ggplot(long_Zscore_D21_T21_ploidy_TpmAvg_sort_factor, aes(x = zscore, y = SYMBOL, fill = ploidy)) +
  geom_density_ridges(alpha = 0.6, scale = 0.98) +
  theme_ridges() + 
  xlim(c(-2.6, 2.6)) +
  scale_fill_manual(
    values = c("D21" = "red1", "T21" = "blue1"),  # Define colors for D21 and T21
    name = "Ploidy"                              # Add a legend title
  ) +
  theme(
    axis.text.y = element_text(size = 8, color = "black", face = "bold"),
    axis.title.y = element_text(size = 9, color = "black", face = "bold"),
    axis.title.x = element_text(size = 9, color = "black", face = "bold"),
    legend.position = "right"                   # Place legend on the right
  ) +
  ylab("Genes involved in transcriptional regulation") +
  xlab("Z-scored normalized counts")

# Heatmap for average TPM values
p_heatmap <- ggplot(heatmap_data, aes(x = 1, y = SYMBOL)) +
  geom_tile(aes(fill = D21), color = "white") +  # Heatmap for D21
  geom_tile(aes(x = 2, fill = T21), color = "white") +  # Heatmap for T21
  scale_fill_gradient(
    low = "white", high = "grey30", 
    limits = c(0, 3000),         # Set the color range from 0 to 3000
    oob = scales::squish,        # Values above 3000 are squished to the max color
    na.value = "grey90"          # Set NA values to a light grey
  ) +
  scale_x_continuous(breaks = c(1, 2), labels = c("D21", "T21")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add a black border
  ) +
  labs(fill = "Avg TPM", x = NULL)

# Combine plots
final_plot <- p_ridge + p_heatmap +
  plot_layout(widths = c(4, 1))  # Adjust width ratio between the ridgeline plot and the heatmap

# Display the plot
final_plot
