###########################
#                         #
#     Get climate data    #
#                         #
###########################

library(sp)
library(raster)
library(geodata)
library(terra)

samples = read.table("accession_data.txt", header = T)
head(samples)
lon<-samples$lon
lat<-samples$lat

# Extract coordinate data
xy <- samples[, c("lon", "lat")]
str(xy)
# Load BioClim data. The following command checks if the data is present at the 
# specified path. If the data is not present, it will be downloaded. Don't Foget to change the path to a folder where you stored your climate data
biodata = worldclim_global(var = "bio", res = 10, "C:/Users/flo19/OneDrive/Uni/3. Semester/305 - Population & Quantitative Genetics/Project_Cologne")
biodata # inspect the data

# Read UV_index data as ascii, downloaded from: https://www.ufz.de/gluv/index.php?en=32367
# Don't Foget to change the path to a folder where you stored your climate data
uvindex = raster("C:/Users/flo19/OneDrive/Uni/3. Semester/305 - Population & Quantitative Genetics/Project_Cologne/get_data_climate_variables/get_data_climate_variables/56459_UVB1_Annual_Mean_UV-B.asc")
# plot UV index
plot(uvindex, main = "Global UV index")
# Extract UV index using a xy coordinates dataframe
uv_index = extract(uvindex, xy, df = T)
# Write UV data as a text file
write.table(uv_index, file = "uv_index.txt", row.names = FALSE, col.names = FALSE)

#########################
#                       #
#   Terminal work       #
#                       #
#########################



#########################
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data <- read.table("sample_project_GWAS2.assoc.txt", header = TRUE)

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = (0.05/dim(gwas_data)[[1]])


#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data, chr="chr", bp="ps", snp="rs", p="p_wald", 
          suggestiveline = -log10(pval_bonf), genomewideline = FALSE, 
          annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


## Alternative P-value correction (FDR)


# Make the Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold, you can adjust this if needed

# Add the adjusted p-values as a new column to the data frame
gwas_data$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE, 
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"))


### --> No result









##########################################
#          R script for running         ##
#             Gene Enrichment           ##
#                Analysis               ##
##########################################

install.packages("BiocManager“)
library(BiocManager)
BiocManager::install("") 


library(BiocManager)
library(KEGGREST)
library(org.At.tair.db)
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)
library(scales)



