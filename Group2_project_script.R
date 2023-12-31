################################################################################################## 
##  Welcome to this R script designed to explore genetic variation within specific genomic      ##
##  regions, focusing on key metrics like FST, Tajima's D, Pi, and Neutrality stats.However,    ##
##  it is always good practice to perform some exploratory analyses such as PCA, Admixture.     ##
##  Our analysis will zoom in on particular genomic regions, e.g., genes associated with plant  ##
##  responses to pathogens, shedding light on the intricate genetic landscape governing these   ##
##  critical biological processes. Let's embark on a journey to unravel the secrets encoded     ##
##  in the genome and understand the fine genetic variations within these crucial regions.      ##
##################################################################################################

###
# This script uses multiple R packages. The first step is to install all required packages if not installed earlier.
# Note: The R package "PopGenome" is not available in the R repository. 
# Download and install the most recent version of PopGenome ("PopGenome_2.7.5.tar.gz") using 
# the following link: https://cran.r-project.org/src/contrib/Archive/PopGenome/
# 
###

BiocManager::install("PopGenome")

# load all required packages

library(PopGenome) 
library(dplyr)
library(ggplot2)
library (readr)
library(tibble)
library(vcfR)
library(adegenet)
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(ggrepel)
library(gplots)
library(StAMPP)
library(RColorBrewer)

# Set the R working directory to the location where you have stored your indexed 'input.vcf' file 
# and the 'genomic_positions.bed' file." Don't forget to change the path to working directory in your system

setwd("~/Downloads/Group2_data/")



########################################################################
##  Let' conduct dimentionality reduction directly by using VCF file  ##
##              Principle Component Analysis (PCA)                    ##
########################################################################

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

vcf_file <- read.vcfR("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")

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
ind_names<- read.table("group_2_accession_names.txt", header= FALSE)
str(ind_names)

# Add a new column "ind" to pca_axis using the values from ind
pca_axis$ind <- ind_names$V1

# You can also directly provide a list of names
population_labels <- c("ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S", "ITA_S",
                       "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M", "ITA_M",
                       "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", "ITA_N", 
                       "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", "GER_S", 
                       "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", "GER_M", 
                       "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N", "GER_N",
                       "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", "SWE_S", 
                       "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", "SWE_M", 
                       "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N", "SWE_N") # Replace ... if you have different population labels

# Add a new column named "Population" to your data frame
pca_axis$Population <- population_labels

pop<-pca_axis$Population

# remake data.frame
pca_2 <- as.tibble(data.frame(pca_axis, population_labels))

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

# plot pca PC1 and PC2
b <- ggplot(pca_2, aes(Axis1, Axis2, col = pop)) + 
  geom_point(size = 4) +
  scale_colour_manual(values = c("green", "green3", "green4", "red", "red3", "red4", "blue", "blue3", "navy")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Italy (ITA), Germany (GER) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        plot.margin = margin(20, 20, 20, 20))   # Adjust margin for title

# Print the plot
print(b)

# plot pca PC1 and PC3

b <- ggplot(pca_2, aes(Axis1, Axis3, col = pop)) + 
  geom_point(size = 4) +
  scale_colour_manual(values = c("green", "green3", "green4", "red", "red3", "red4", "blue", "blue3", "navy")) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Italy (ITA), Germany (GER) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5),  # Center title
        plot.margin = margin(20, 20, 20, 20))   # Adjust margin for title

# Print the plot
print(b)

# Load required libraries
library(ggplot2)
library(ggrepel)  # for label repulsion

# Assuming you have a data frame called 'pca' with columns 'PC1', 'PC2', 'pop', and 'ind'
# 'pop' contains population labels and 'ind' contains individual names

# Create a scatterplot

b <- ggplot(pca_2, aes(Axis1, Axis2, col = pop, label = ind)) + 
  geom_point(size = 2) + # size of points
  geom_text_repel(show.legend = FALSE, size = 2) +  # , max.overlaps = Inf >> Use is option to plot all labels
  scale_colour_manual(values = c("green", "green3", "green4", "red", "red3", "red4", "blue", "blue3", "navy")) + 
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 2), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 2), "%)")) +
  ggtitle("Arabidopsis thaliana Accessions from Italy (ITA), Germany (GER) and Sweden (SWE)") +
  theme(plot.title = element_text(hjust = 0.5))  # Adjust the hjust value for centering

# Print the plot
print(b)

##
#PCA Visualization using plotly
#https://plotly.com/r/pca-visualization/ #, "green", "yellow3", "yellow","orange","purple"
library(plotly)
#3D plots
fig <- plot_ly(pca_2, x = ~Axis1, y = ~Axis2, z = ~Axis3, color =pop, colors = c("green", "green3", "green4", "red", "red3", "red4", "blue", "blue3", "navy") ) %>%
  add_markers(size = 12)

fig <- fig %>%
  layout(
    title = "Arabidopsis thaliana Accessions from Italy (ITA), Germany (GER) and Sweden(SWE)",
    scene = list(bgcolor = "#e5ecf6")
  )

fig

# make plot background transparent

fig <- fig %>%
  layout(title = "Arabidopsis thaliana Accessions from Italy (ITA), Germany (GER) and Sweden(SWE)") %>% 
  layout(plot_bgcolor = 'rgb(0, 0, 0,0)') %>% 
  layout(paper_bgcolor = 'rgb(0, 0, 0,0)') %>% 
  layout(margin = list(l = 80, r = 80, b = 80, t = 80)) %>%  # Adjust margin values as needed
  layout(legend = list(bgcolor = 'rgba(0,0,0,0)'))  # Set legend background color to transparent

fig


###
#Visualize a subset of the principal components
cumsum(pve$pve)

tit = 'Total Explained Variance = 51.8 %'

axis = list(showline=FALSE,
            zeroline=FALSE,
            gridcolor='#ffff',
            ticklen=4)

fig <- pca_2 %>%
  plot_ly() %>%
  add_trace(
    type = 'splom',
    dimensions = list(
      list(label='PC1 (18.4%)', values=~Axis1),
      list(label='PC2 (13.2%)', values=~Axis2),
      list(label='PC3 (10.5%)', values=~Axis3),
      list(label='PC4 (9.7%)', values=~Axis4)
    ),
    color = ~pop, colors = c("green", "green3", "green4", "red", "red3", "red4", "blue", "blue3", "navy"), #, "green", "yellow3", "yellow","orange","purple"),
    marker = list(
      size = 7
    )
  ) %>% style(diagonal = list(visible = F)) %>%
  layout(
    title= tit,
    hovermode='closest',
    dragmode= 'select',
    plot_bgcolor='rgba(240,240,240, 0.95)',
    xaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    yaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    xaxis2=axis,
    xaxis3=axis,
    xaxis4=axis,
    yaxis2=axis,
    yaxis3=axis,
    yaxis4=axis
  )
options(warn=-1)
fig


########################################################################
##            Stampp to calculate FST between populations             ##
##                         and individuals                            ##
########################################################################


# Convert VCF to genlight object
genlight_vcf <- vcfR2genlight(vcf_file)

# Read population information from file
pop <- read.table("pop2.txt", header = TRUE)
str(pop)

# Extract population data
pop1 <- pop$pop
pop2 = as.factor(pop1)
genlight_vcf$pop = pop2

# Convert genlight to stampp object
stampp_vcf = stamppConvert(genlight_vcf, type = "genlight")

# Calculate FST between populations
stamppFst = stamppFst(stampp_vcf, nboots = 100, percent = 95, nclusters = 8)
stamppFst_matrix = as.matrix(stamppFst$Fsts)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppFst_matrix) <- 0
stamppFst_matrix[upper.tri(stamppFst_matrix)]  <- t(stamppFst_matrix)[upper.tri(stamppFst_matrix)]
# Optional: order the names
stamppFst_matrix = stamppFst_matrix[order(row.names(stamppFst_matrix)), order(colnames(stamppFst_matrix))]

# Make an FST heatmap
heatmap(stamppFst_matrix,
        symm = TRUE,
        margins = c(10, 10),
        main = "Genetic Divergence (FST) b/w A. thaliana Pop from Italy, Germany & Sweden ")

#####################################################################
########## calculate genetic distance between individuals ###########
#####################################################################

stamppNeisD = stamppNeisD(stampp_vcf, pop = FALSE)
stamppNeisD_matrix = as.matrix(stamppNeisD)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppNeisD_matrix) <- 0
stamppNeisD_matrix[upper.tri(stamppNeisD_matrix)]  <- t(stamppNeisD_matrix)[upper.tri(stamppNeisD_matrix)]
heatmap(stamppNeisD_matrix)

# add row names
colnames(stamppNeisD_matrix) <- rownames(stamppNeisD_matrix)

# Create a heatmap with symmetric color scale
heatmap(stamppNeisD_matrix,
        symm = TRUE,
        main = "Genetic Divergence (FST) b/w A. thaliana individuals from Italy, Germany & Sweden ")

########################################################################
##                           Plot admixture                           ##
##       using input produced using admixture.sh in terminal          ##
########################################################################

library(stringr)
library(ggplot2)
library(dplyr)

#set working directory
setwd("~/Downloads/Group2_data/admixture")
getwd()

cv <- read.table("cross_validation.txt")

# Analyze the cross-validation results Then, add a K-cluster column indicating the number of K you test and select only two columns of interest, CV and K.

cv$K <-gsub("[\\(\\)]", "", regmatches(cv$V3, gregexpr("\\(.*?\\)", cv$V3)))
CV <- select(cv, V4,K)
CV$K <- as.numeric(sub("K=", "", CV$K))

# Rename your two columns CV and K-cluster
colnames(CV) <- c("CV","K")

# Do a graph showing the cross validation results. Then select the optimal number of clusters regarding :
# the lowest cross validation error
# when the cross-validation error decrease the most

graph_title="Cross-Validation plot"
x_title="K"
y_title="Cross-validation error"
graph_1<-ggplot(CV,aes(x=K,y=CV))
graph_1+geom_line()+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10))+
  labs(title=graph_title)+
  labs(x=x_title)+
  labs(y=y_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))

#Save the graph
#ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
#dev.off()

# One can modify the color palette to include the colors yellow2, blue, orange, and red2 and its shades 
# using the palette function and the colors parameter. You can define a custom palette
# with these colors. Here's how you can do it:
# Define a custom color palette
my_colors <- c("red","blue","yellow2","orange", "red2", "yellow3", "steelblue","royalblue", "yellowgreen", "goldenrod","darkblue")

# Set the custom palette
palette(my_colors)

#load admxiture function
source("admixFun.R") 

#list Q files and sort for K
files <- list.files("~/Downloads/Group2_data/admixture", full = TRUE, pattern = "Q")
files <- files[order(as.numeric(sub('.*a.thaliana_admix.(\\d+)\\.Q$', '\\1', files)))]

#population file
pop <- scan("~/Downloads/Group2_data/admixture/pop.txt",what="df",na="")
table(pop) 

# possible K
Kall <- 1:10

## read Qs
allQ <- list()
for(K in Kall)
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))

#exact ordering (slow for large K)
plotMulti(allQ,Kall=2:10,reorder=1,pop,fast=T,lwd=1,lty=1)

#make smaller line and change type to solid line
plotMulti(allQ,Kall=2:10,reorder=1,pop,fast=T,lwd=1,lty=1)

plotMulti(allQ,Kall=2:10,reorder=1,pop,fast=T,lwd=0) 


########################################################################
##         Estimate and plot Fst and Tajima'D and Neutrality          ##
##                       stats using PopGenome                        ##
########################################################################

setwd("~/Downloads/Group2_data")

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

at.VCF <- read.vcfR("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")

#get start and stop positions of your vcf file
head(getFIX(at.VCF))

tail(getFIX(at.VCF))

#Reading in the Arabidopsis vcf
#https://rdrr.io/cran/PopGenome/man/readVCF.html

At_Chr <-readVCF("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz",numcols=89,tid="1", frompos=306735, topos=25089536, include.unknown = TRUE) 

#To the class of object At_Chr
class(At_Chr)

At_Chr ###this is your genome.class data. You can push all genomics analysis into one object. 

####Examining the variant data
#Remember, you can look the data we have read in using the following command:
get.sum.data(At_Chr)

#From the n.biallelic.sites we can see there are  428 bilallelic SNPs and from n.polyallelic.sites,
#there are 0 positions with more than two alleles. So in total we have:
At_Chr@n.biallelic.sites + At_Chr@n.polyallelic.sites

#To see what slots in Genome class
show.slots(At_Chr)#

#To check total number of sites
At_Chr@n.sites

#To check starting position and last position of genome class
At_Chr@region.names

#####Define populations in your dataset####

#If you look at this, you will only see a blank list. So we need to supply our population data to
#the ch4 object. To make naming our populations simple, we will read in some external data. 
At_Chr@populations # check for population data

library (readr)

###population data is stored in data.frame that has two columns, one column for individual name, one column for pop

population_info <- read_delim("pop1.txt", delim = "\t")
# now get the data for the populations
populations <- split(population_info$sample, population_info$pop)

# now set 
At_Chr <- set.populations(At_Chr, populations, diploid = T)
##check if it worked
At_Chr@populations

####Setting up sliding windows###

#Per-SNP estimates of statistics such as ?? can often be extremely noisy when you are calculating them on
#very large numbers of markers. As well as this, there are issues with the fact that SNP positions in close
#proximity are not always independent due to recombination - this is a theme we will return too shortly. 
#So for this reason, it is often better to use a sliding-window approach - i.e. split the genome into
#windows of a particular size and then calculate the mean for a statistic within that window.

#We know already that chromosome 4 is 18584000 bp long, so we can get an idea of how many sliding windows
#we would generate by using some R code. We'll set our sliding window to be 100 bp wide.
#We will also set a step or jump for our window of 50 bp.
# set chromosome size
chr <- 26222517

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start <- seq(from = 1, to = chr, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size

# no windows start before the end of chromosome 4
sum(window_start > chr)
# but some window stop positions do occur past the final point
sum(window_stop > chr)

# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < chr)]
window_stop <- window_stop[which(window_stop < chr)]

chr - window_stop[length(window_stop)]

# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

#https://rdrr.io/cran/PopGenome/man/sliding.window.transform-methods.html
# make a sliding window dataset
At_sw <- sliding.window.transform(At_Chr, width = 100, jump = 50, type = 2)


#######Calculating sliding window estimates of nucleotide diversity and differentiation#####
#Now that we have set up the data, the population information and the sliding windows, it is quite
#straightforward for us to calculate some statistics we are interested in. In this case, we are going
#to calculate nucleotide diversity (i.e. ??) and FST. We will also generate a third statistic, d_XY_,
#which is the absolute nucleotide divergence between two populations.

#First we will calculate ??. Handily, the following command also sets up what we need for d_XY_.

# calculate diversity statistics
At_sw <- diversity.stats(At_sw, pi = TRUE)


#Next we will calculate FST, which again is very straight forward with a single command.

### calculate diversity statistics
At_sw <- F_ST.stats(At_sw, mode = "nucleotide")

#Note that here we use mode = "nucleotide" to specify we want it to be calculated sliding averages
#of nucleotides, rather than using haplotype data, which is the alternative. And that's it for 
#calculating the statistics! As you will see in the next section, extracting them from the 
#At_sw object is actually more difficult than generating them

####Extracting statistics for visualization####

#Since we ran our analysis on a sliding-window basis, we should have estimates of ??, FST and d_XY_ for
#each window. What we want to do now is extract all our statistics and place them in a single data.frame
#for easier downstream visualisation - this will let us identify how these statistics are interrelated.

#First of all, we will get the nucleotide diversity data.

# extract nucleotide diversity and correct for window size
nd <- At_sw@nuc.diversity.within/100

#This is straightforward, but remember also that our estimates need to be corrected for 
#window size - so we divide them by 100 bp here. We should also add the population names
#to each of them, since they are not already set.

# make population name vector
pops <- c("ITA","GER","SWE")
# set population names
colnames(nd) <- paste0(pops, "_pi")

# extract fst values
fst <- t(At_sw@nuc.F_ST.pairwise)
#Note that here, we need to use t() to transpose the F_ST matrix so that each column is a pairwise
#comparison and each row is an estimate for a genome window. Since F_ST is pairwise, the column
#names are also quite different and will also be the same for d_XY_, which is also a pairwise measure.

#So now we are ready to extract our final statistic, d_XY_. We can do this in a similar way to how we
#handled the FST data.

# extract dxy - pairwise absolute nucleotide diversity
dxy <- get.diversity(At_sw, between = T)[[2]]/100
#As with nucleotide diversity, we also corrected d_XY_ for the window size.

#Now we sort out the column names for our FST and d_XY_ data. This is where our R skills come in use!
#We will need to use some R-based string manipulation. The column names are identical for both datasets,
#so we will take the first one and use the sub function to replace the population names.

# get column names 
x <- colnames(fst)
# replace all occurrences of pop1 with ITA
x <- sub("pop1", "ITA", x)
# does the same thing as above but by indexing the pops vector
x <- sub("pop1", pops[1], x)
# look at x to confirm the replacement has occurred
x


# get column names 
x <- colnames(fst)
# does the same thing as above but by indexing the pops vector
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)

# replace forward slash
x <- sub("/", "_", x)
# look at x to confirm the replacement has occurred
x

#Now all that we need to do is make clear these names are for either FST or d_XY_. 
#The best way to do this is to append a suffix to our vector of pairwise comparison names.
#We can do this using paste0
paste0(x, "_fst")
paste0(x, "_dxy")

#So this function allows us to join strings together in a character vector. Very useful. 
#Next we will actually change the column names of our two data matrices, before we put 
#everything together in our final dataset.
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#Ok so now that our ??, FST and d_XY_ datasets are ready, we can combine them all together
#with our windows information from earlier into a big dataset.
library(tibble)
At_data <- as.tibble(data.frame(windows, nd, fst, dxy))


#####Visualizing the data - distributions#####

#For the purposes of this session, we will focus mainly on the difference between Italian, German and Swedish
#Arabidopsis pop.
#For example, let's say we want to look at mean nucleotide diversity, we can do that like so:

# select nucleotide diversity data and calculate means
At_data %>% select(contains("pi")) %>% summarise_all(mean)
#we used select and contains to select columns from our main dataset that contain 
#pi - i.e. nucleotide diversity columns. We then used summarise_all and mean to calculate
#the mean value for all of the columns we selected.

# To plot this we need to use "gather" on the data
library(tidyr)
pi_g <- At_data %>% select(contains("pi")) %>% gather(key = "populations", value = "pi")

# make a boxplot
library(ggplot2)
a <- ggplot(pi_g, aes(populations, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
a

# Taking the logarithm (in this case, log base 10) of the values can be useful when dealing with data
# that spans several orders of magnitude. This transformation helps in visually emphasizing relative
# differences in the data, especially when there are large variations in scale.

pi_g$log_pi <- log10(pi_g$pi)

a <- ggplot(pi_g, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "blue")) +  # Box fill colors
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)")

a

#This makes it much clearer how nucleotide diversity differs among the populations.

#When comparing two boxplots to determine if they are statistically different, one can perform
# statistical tests such as the t-test or Wilcoxon rank-sum test. 
# Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(log_pi ~ populations, data = pi_g)

# Print the test result
print(wilcox_test_result)

# Add p-value to the plot
a + annotate("text", x = 1.5, y = max(pi_g$log_pi), label = paste("p =", format.pval(wilcox_test_result$p.value, digits = 3)))

#####Visualizing patterns along the chromosome ####
#Let's have a look at how FST between Spanish and Swedish populations varies along chromosomes.
#We can do this very simply with ggplot.

a <- ggplot(At_data, aes(mid/10^6, ESP_SWE_fst)) + geom_line(colour = "red")
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a + theme_light()


#to plot ??, FST and d_XY_ to examine how they co-vary along the genome. 
#This requires a bit of data manipulation, but is relatively straightforward. We will break it down into steps.
# select data of interest
hs <- At_data %>% select(mid, ESP_pi, SWE_pi, ESP_SWE_fst, ESP_SWE_dxy)

# To set Fst values smaller than zero to zero in the specified columns of a data frame using dplyr and
# the pipe operator %>%, you can use the mutate function along with across

hs <- At_data %>%
  select(mid, ESP_pi, SWE_pi, ESP_SWE_fst, ESP_SWE_dxy) %>%
  mutate(across(c(ESP_SWE_fst, ESP_SWE_dxy), ~ ifelse(. < 0, 0, .)))

# use gather to rearrange everything
hs_g <- gather(hs, -mid, key = "stat", value = "value")

#All "gather" function does is collapse everything so we can plot them efficiently. We use -mid to tell 
#the function we want to leave this out of the gathering and use key = stat to make it clear
#we are arranging our data by the statistics we have calculated, value = value is just a name
#for the values of each of our statistics.

#Now we can easily plot everything together like so:
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + xlab("Position (Mb)")
a + theme_light()

# To take the logarithm of the value variable in your ggplot code, you can use the log10() function
# within the aes() mapping.
hs_g$log_value <- log10(hs_g$value)

a <- ggplot(hs_g, aes(mid/10^6, log_value, colour = stat)) + geom_line()
a <- a + xlab("Position (Mb)")
a + theme_light()


#OK so it should be immediately obvious that this plot is really unhelpful. We see the FST data again,
#but since that is on such a different scale to estimates of ?? and d_XY_, we can't see anything! Instead,
#it would make a lot more sense to split our plot into facets - i.e. a plot panel for each statistic. 
#This is simple with the ggplot function facet_grid. We will construct our plot first and then breakdown

###
#what facet_grid actually does?
###
# construct a plot with facets
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")

#The facet_grid function allows us to split our data across panels for quick and easy visualization.
#In this case, we split our data by the stat variable - we used stat~. to specify we want this done
#by rows (compare with .~stat for the column equivalent). We also specified that we wanted the scales
#on our y-axes to vary with scales = free_y.

#However, before we examine our plot in detail, it would also be easier if we rearranged everything 
#so FST came at the top, ?? beneath it and then finally, d_XY_. How can we do that? Well we need to
#reorder the stat factor in our hs_g dataset.
# first make a factor
x <- factor(hs_g$stat)
# then reorder the levels
x <- factor(x, levels(x)[c(3, 1, 4, 2)])

# add to data.frame
hs_g$stat <- x
#This looks a little complicated, but in the square brackets above we simply rearranged what order
#our facets are displayed. We can replot our figure to demonstrate this:

# construct a plot with facets
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")

#### calculate neutrality statistics####
At_sw <- neutrality.stats(At_sw)

get.neutrality(At_sw)

#Let's look at the first population [[1]].
get.neutrality(At_sw)[[1]]

#Let's look at the first population [[2]].
get.neutrality(At_sw)[[2]]

#extract Tajma's D
td <- At_sw@Tajima.D/100

# set population names
colnames(td) <- paste0(pops, "_td")

####


###Delimitate windows on chromosome

# set chromosome start and end position
chri<-1005
chrl <- 26222517

library(tibble)

#as_tibble: Coerce lists and matrices to data frames
ara_data <- as.tibble(data.frame(windows, td,nd))
nrow(windows)
nrow(nd)
(chrl-chri)/50
nrow(ara_data)
head(ara_data)
ara_data %>% select(contains("pi")) %>% summarise_all(mean)##mean pi across all windows is for Spanish pop 0.000790

### load selected positions from chromosome e.g., gene 4 5kb upstream and down stream of Defense related genes
# 
bed<-read.table("At_defense_only.bed")
head(bed)

colnames(bed)<-c("chr", "begin","end")
DF<-vector(length=nrow(ara_data))
ara_data <- as.tibble(data.frame(windows, nd, DF))###if you only want to look at pi
ara_data <- as.tibble(data.frame(windows, nd, td, DF))##if you want to look at tajima D and nucleotide diversity

for (i in 2:nrow(bed)){ara_data$DF[which(ara_data$start>bed$begin[i]&ara_data$stop<bed$end[i]) ]<-"DF"}##each window that overlaps a DF is tagged
ara_data$DF<-as.factor(ara_data$DF)
summary(ara_data)

####
# Spain
###

##Kolmogorov smirnov test - compare the distributions of Pi
sub1<-ara_data$ESP_pi[ara_data$DF=="DF"]
sub2<-ara_data$ESP_pi[ara_data$DF!="DF"]
ks.test(sub1, sub2)###difference is very significant if windows are small, otherwise not. 

#Draw Density plot "Pi"
plot(density(log(ara_data$ESP_pi)), main="Distribution log Pi")
lines(density(log(ara_data$ESP_pi[ara_data$DF=="DF"])), col="red")

##Kolmogorov smirnov test - compare the distributions of Tajima's D
sub1<-ara_data$ESP_td[ara_data$DF=="DF"]
sub2<-ara_data$ESP_td[ara_data$DF!="DF"]
ks.test(sub1, sub2)

# Draw Density plots "Tajima's D"
plot(density((ara_data$ESP_td), na.rm=T), main="Distribution Tajima D")
lines(density((ara_data$ESP_td[ara_data$DF=="DF"]), na.rm = T), col="red")

p<-ggplot(ara_data, aes(x=ESP_td, fill=DF))
p+geom_density(alpha=0.4)

# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_data$ESP_td, na.rm = TRUE)
lowest_x

p <- ggplot(ara_data, aes(x = ESP_td, fill = DF)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.03, max(ara_data$ESP_td, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution Tajima D")
# Center the title
p <- p + ggtitle("Distribution Tajima D") +
  theme(plot.title = element_text(hjust = 0.5))  # Adjust the hjust value for centering
# plot distribution
p


##Plot along chromosome using ggplot function
sub1<-(ara_data[ara_data$DF=="DF",])
sub2<-ara_data[ara_data$DF!="DF",]
p<-ggplot(sub2, aes(mid,ESP_pi))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)

p <- ggplot(sub2, aes(mid,ESP_td))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)

####
# Sweden
###

##Kolmogorov smirnov test - compare the distributions of Pi
sub1<-ara_data$SWE_pi[ara_data$DF=="DF"]
sub2<-ara_data$SWE_pi[ara_data$DF!="DF"]
ks.test(sub1, sub2)###difference is very significant if windows are small, otherwise not. 

#Draw Density plot "Pi"
plot(density(log(ara_data$SWE_pi)), main="Distribution log Pi")
lines(density(log(ara_data$SWE_pi[ara_data$DF=="DF"])), col="red")

##Kolmogorov smirnov test - compare the distributions of Tajima's D
sub1<-ara_data$SWE_td[ara_data$DF=="DF"]
sub2<-ara_data$SWE_td[ara_data$DF!="DF"]
ks.test(sub1, sub2)

# Draw Density plots "Tajima's D"
plot(density((ara_data$SWE_td), na.rm=T), main="Distribution Tajima D")
lines(density((ara_data$SWE_td[ara_data$DF=="DF"]), na.rm = T), col="red")

p<-ggplot(ara_data, aes(x=SWE_td, fill=DF))
p+geom_density(alpha=0.4)

##
# Base R plot
plot(density(ara_data$SWE_td, na.rm = TRUE), main = "Distribution Tajima D")
lines(density(ara_data$SWE_td[ara_data$DF == "DF"], na.rm = TRUE), col = "red")

# ggplot version

# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_data$SWE_td, na.rm = TRUE)
lowest_x

p <- ggplot(ara_data, aes(x = SWE_td, fill = DF)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.010, max(ara_data$SWE_td, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution Tajima D")
# plot distribution
p

##Plot along chromosome using ggplot function
sub1<-(ara_data[ara_data$DF=="DF",])
sub2<-ara_data[ara_data$DF!="DF",]
p<-ggplot(sub2, aes(mid,SWE_pi))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)

p<-ggplot(sub2, aes(mid,SWE_td))
p+geom_point(size=2)+geom_point(data=sub1, color="red", size=3)


##########################
# FST ESP_SWE           ##
#                       ##
##########################

#as_tibble: Coerce lists and matrices to data frames
ara_data2 <- as.tibble(data.frame(windows, fst))
nrow(windows)
nrow(fst)
(chrl-chri)/50
nrow(ara_data2)
head(ara_data2)
ara_data2 %>% select(contains("fst")) %>% summarise_all(mean)

### load selected positions from chromosome -> flowering time genes ####

bed2<-read.table("At_defense_only.bed")
head(bed2)

colnames(bed2)<-c("chr", "begin","end")
DF2<-vector(length=nrow(ara_data2))
ara_data2 <- as.tibble(data.frame(windows, fst, DF2))

for (i in 2:nrow(bed2)){
  ara_data2$DF22 <- "all"
} #horrible but works

for (i in 2:nrow(bed2)){
  ara_data2$DF2[which(ara_data2$start>bed2$begin[i]&ara_data2$stop<bed2$end[i])]<-"DF2"
}

DF2<-vector(length=nrow(ara_data2))

ara_data2$DF2<-as.factor(ara_data2$DF2)
summary(ara_data2)


Defense <- ara_data2 %>% filter(DF2 == "DF2")
Defense_fst <- Defense %>% filter(ESP_SWE_fst >= 0)

a <- ggplot(Defense_fst, aes(mid/10^6, ESP_SWE_fst)) + geom_line(colour = "red")
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
a + theme_light()

# plot with ara_data FST (<=0 values not removed) and FLOWER_fst (all flowering time FSTs also with <=0 values not removed)
tip <- ggplot() + 
  geom_line(data=ara_data2, aes(mid/10^6, ESP_SWE_fst), colour = "blue") + 
  geom_line(data=Defense, aes(mid/10^6, ESP_SWE_fst), colour="pink")
tip <- tip + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
tip + theme_light()

# remove fst values <=0 
ara_d2 <- ara_data2 %>% filter(ESP_SWE_fst >= 0)
ara_d2

# calculate means
mean_fst <- mean(ara_d2$ESP_SWE_fst)
mean_defense <- mean(Defense_fst$ESP_SWE_fst)

ks.test(ara_d2$ESP_SWE_fst, Defense_fst$ESP_SWE_fst) #p-value: 1

#outliers 95% quantile
threshold_95 <- quantile(Defense_fst$ESP_SWE_fst[Defense_fst$DF2=="DF2"], 0.975, na.rm = T)
Defense_fst <- Defense_fst %>% mutate(outlier_95 = ifelse(Defense_fst$ESP_SWE_fst > threshold_95, "outlier", "background"))

#outliers 99% quantile
threshold_99 <- quantile(Defense_fst$ESP_SWE_fst[Defense_fst$DF2=="DF2"], 0.995, na.rm = T)
Defense_fst <- Defense_fst %>% mutate(outlier_99 = ifelse(Defense_fst$ESP_SWE_fst > threshold_99, "outlier", "background"))

# plot with ara data and Defense_fst (all fst values below 0 removed)
top <- ggplot() + 
  geom_point(data=ara_d2, aes(mid/10^6, ESP_SWE_fst), colour = "lightblue") + 
  geom_point(data=Defense_fst, aes(mid/10^6, ESP_SWE_fst), colour="blue") +
  geom_point(data=Defense_fst[Defense_fst$outlier_95 == "outlier",], aes(mid/10^6, ESP_SWE_fst), color="orange") +
  geom_point(data=Defense_fst[Defense_fst$outlier_99 == "outlier",], aes(mid/10^6, ESP_SWE_fst), color="red") +
  geom_hline(yintercept = mean_fst) +
  geom_hline(yintercept = mean_defense, colour="orange")

top <- top + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
top + theme_light()


# Draw Density plots

plot(density((ara_d2$ESP_SWE_fst), na.rm=T), main="Distribution FST", )
lines(density((Defense_fst$ESP_SWE_fst[Defense_fst$DF2=="DF2"]), na.rm = T), col="red")

p<-ggplot(ara_d2, aes(x=ESP_SWE_fst, fill=DF2))
p+geom_density(alpha=0.4)


# To estimate the lowest value of the x-axis for a density plot
lowest_x <- min(ara_d2$ESP_SWE_fst, na.rm = TRUE)
lowest_x

p <- ggplot(ara_d2, aes(x = ESP_SWE_fst, fill = DF2)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(limits = c(-0.05, max(ara_d2$ESP_SWE_fst, na.rm = TRUE))) +  # Set x-axis limit
  ggtitle("Distribution ESP-SWE FST")
# plot distribution
p

# ggplot version with log scale for x-axis
p <- ggplot(ara_d2, aes(x = ESP_SWE_fst, fill = DF2)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +  # Set log scale for x-axis
  ggtitle("Distribution ESP-SWE FST")

# plot distribution
p

############################################################
### To check where these FST outliers are in the genome ####
###                                                     ####
############################################################

# Outliers 95% quantile whole genome
threshold_95 <- quantile(ara_d2$ESP_SWE_fst, 0.95, na.rm = T)
ara_d2 <- ara_d2 %>% mutate(outlier_95 = ifelse(ara_d2$ESP_SWE_fst > threshold_95, "outlier", "background"))

# Outliers 99% quantile whole genome
threshold_99 <- quantile(ara_d2$ESP_SWE_fst, 0.99, na.rm = T)
ara_d2 <- ara_d2 %>% mutate(outlier_99 = ifelse(ara_d2$ESP_SWE_fst > threshold_99, "outlier", "background"))

out <- ara_d2 %>% filter(outlier_95 == "outlier")
out2 <- out %>% filter(DF2 == "DF2")
out2 #no 99% outliers but 2 95% outliers which both correspond to sucrose synthase 3

print(out2,n=50)

