### Load packages----

library(vcfR) # reading VCF
library(adegenet) # genind object for PCA
library(tibble)
library(ggplot2) # graphing

# Set the R working directory to the location where you have stored your indexed 'input.vcf' file 
# and the 'genomic_positions.bed' file." Don't forget to change the path to working directory in your system

setwd("~/Downloads/Group2_data/")

### Principle Component Analysis (PCA)----

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

vcf_file <- read.vcfR("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf")
df <- read.csv("accession_data.csv")

# Convert VCF to genind object
genind_vcf <- vcfR2genind(vcf_file)


# Scale genind object for PCA
genind_vcf_scaled = scaleGen(genind_vcf, NA.method = "mean")

# Perform PCA
pca <- dudi.pca(genind_vcf_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

# Check PCA dimensions
axis_all = pca$eig * 100 / sum(pca$eig)
barplot(axis_all[1:10], main = "PCA eigenvalues")

str(pca)
summary(pca)

#ind <- pca$rank
pca_axis <- pca$li
pca_eigenval <- pca$eig[1:10]
str(pca_axis)

getwd()

# set names
ind_names <- df$accession

str(ind_names)

# Add a new column "ind" to pca_axis using the values from ind
pca_axis$ind <- ind_names

# You can also directly provide a list of names
population_labels <- df$country
subpopulation_labels <- df$subpopulation

# Add a new column named "Population" to your data frame
pca_axis$Population <- population_labels
pca_axis$Subpopulation <- subpopulation_labels

pop <- pca_axis$Population
subpop <- pca_axis$Subpopulation

# remake data.frame
pca_2 <- as.tibble(data.frame(pca_axis, population_labels))
pca_3 <- as.tibble(data.frame(pca_axis, subpopulation_labels))

n <- length(pca_eigenval)
n # use this number PC=1:n

# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = pca_eigenval/sum(pca_eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

str(pca_2)
str(pca_3)

# plot pca PC1 and PC2 (by country)
b <- ggplot(pca_2, aes(Axis1, Axis2, col = pop)) + 
  geom_point(size = 6) +
  scale_colour_manual(breaks = c("ITA", "GER", "SWE"), values = c("red","green","blue")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(col="Countries") +
  ggtitle("A") +
  theme(plot.title = element_text(size = rel(4)),
        axis.text.x = element_text(size = rel(1.5)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)), 
        axis.title.y = element_text(size = rel(1.5), angle = 90))

# Print the plot
print(b)

# plot pca PC1 and PC2 (by subpopulation)
c <- ggplot(pca_3, aes(Axis1, Axis2, col = subpop)) + 
  geom_point(size = 6) +
  scale_colour_manual(breaks = c("Italy_South","Italy_Middle","Italy_North","Germany_South","Germany_Middle","Germany_North","Sweden_South","Sweden_Middle","Sweden_North"), 
                      values = c("red4", "red3", "red", "green4", "green3", "green", "navy", "blue3", "blue")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(col="Subpopulations") +
  ggtitle("B") +
  theme(plot.title = element_text(size = rel(4)),
        axis.text.x = element_text(size = rel(1.5)), 
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)), 
        axis.title.y = element_text(size = rel(1.5), angle = 90))

# Print the plot
print(c)
