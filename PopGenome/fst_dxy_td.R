########################################################################
##         Estimate and plot Fst and Tajima'D and Neutrality          ##
##                       stats using PopGenome                        ##
########################################################################

library(PopGenome) 
library(vcfR)
library(VariantAnnotation)
library (readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dunn.test)

setwd("~/Downloads/Group2_data/")

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

at.VCF <- read.vcfR("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")

vcf <- readVcf("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz")
# Extract chromosome names
chromosomes <- seqlevels(vcf)
# Initialize an empty data frame to store results
chromosome_ranges <- data.frame(CHROM = character(), Start = numeric(), End = numeric(), stringsAsFactors = FALSE)

# Iterate over chromosomes
for (chrom in chromosomes) {
  # Extract positions for the current chromosome
  positions <- start(vcf[seqnames(vcf) == chrom])
  
  # Append results to the data frame
  chromosome_ranges <- rbind(chromosome_ranges, data.frame(CHROM = chrom, Start = min(positions), End = max(positions)))
}

# Display the results
print(chromosome_ranges)
#CHROM    Start      End
#1     1   306735 29974585
#2     2   627319 12632121
#3     3 15661067 23218238
#4     4  1110477 17030281
#5     5   631554 25089536 

chr1start <- 306735
chr1end <- 29974585

chr3start <- 15661067
chr3end <- 23218238

chr5start <- 631554
chr5end <- 25089536

# chromosome 1
At_Chr <- readVCF("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz", numcols=189, tid="1", frompos = chr1start, topos = chr1end, include.unknown = TRUE)

# chromosome 3
At_Chr <- readVCF("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz", numcols=189, tid="3", frompos = chr3start, topos = chr3end, include.unknown = TRUE)

# chromosome 5
At_Chr <- readVCF("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz", numcols=189, tid="5", frompos = chr5start, topos = chr5end, include.unknown = TRUE)


## Examining the variant data
#Remember, you can look the data we have read in using the following command:
#get.sum.data(At_Chr)

#From the n.biallelic.sites we can see there are  568 bilallelic SNPs and from n.polyallelic.sites,
#there are 0 positions with more than two alleles. So in total we have:
#At_Chr@n.biallelic.sites + At_Chr@n.polyallelic.sites # 568

#To check total number of sites
#At_Chr@n.sites

#To check starting position and last position of genome class
#At_Chr@region.names

## Define populations in your dataset

## population data is stored in data.frame that has two columns, one column for individual name, one column for pop
population_info <- read_delim("./PopGenome/pop1.txt", delim = "\t")
# now get the data for the populations
populations <- split(population_info$sample, population_info$pop)

# now set 
At_Chr <- set.populations(At_Chr, populations, diploid = T)
##check if it worked
#At_Chr@populations

## Setting up sliding windows

#We know already that chromosome 4 is 18584000 bp long, so we can get an idea of how many sliding windows
#we would generate by using some R code. We'll set our sliding window to be 100 bp wide.
#We will also set a step or jump for our window of 50 bp.
# set chromosome size
chr <- 29974600 # chr 1

chr <- 23218238 # chr 3

chr <- 25089536 # chr 5

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start <- seq(from = chr3start, to = chr, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size

# no windows start before the end of chromosome
sum(window_start > chr)
# but some window stop positions do occur past the final point
sum(window_stop > chr)

# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < chr)]
window_stop <- window_stop[which(window_stop < chr)]

# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

# make a sliding window dataset
At_sw <- sliding.window.transform(At_Chr, width = 100, jump = 50, type = 2)


## Calculating sliding window estimates of nucleotide diversity and differentiation
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

## Extracting statistics for visualization

# extract nucleotide diversity and correct for window size
nd <- At_sw@nuc.diversity.within/100

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

# Assuming you have a vector pops with four population names
pops <- c("ITA", "GER", "SWE")

# get column names 
x <- colnames(fst)

# Loop through each population and replace the corresponding population name in the column names
for (i in 1:length(pops)) {
  pattern <- paste0("pop", i)
  x <- sub(pattern, pops[i], x)
}

# replace forward slash
x <- sub("/", "_", x)

#Next we will actually change the column names of our two data matrices, before we put 
#everything together in our final dataset.
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#Ok so now that our ??, FST and d_XY_ datasets are ready, we can combine them all together
#with our windows information from earlier into a big dataset.
At_data <- as.tibble(data.frame(windows, nd, fst, dxy))

## Visualizing the data - distributions

# select nucleotide diversity data and calculate means
At_data %>% select(contains("pi")) %>% summarise_all(mean)


# To plot this we need to use "gather" on the data
pi_g <- At_data %>% select(contains("pi")) %>% gather(key = "populations", value = "pi")

pi_g$log_pi <- log10(pi_g$pi)

pi_g$populations <- factor(pi_g$populations,
                           levels = c("ITA_pi", "GER_pi", "SWE_pi"),ordered = TRUE)

a <- ggplot(pi_g, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(breaks = c("ITA_pi", "GER_pi", "SWE_pi"), values = c("red","green","blue")) +  # Box fill colors
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)") +
  labs(fill="Countries")

# When you have more than two groups (three in this case), and you want to compare the 
# central tendency of the distributions, you might consider using analysis of variance (ANOVA) or
# its non-parametric counterpart, the Kruskal-Wallis test.

# Kruskal-Wallis test
kruskal_test_result <- kruskal.test(log_pi ~ populations, data = pi_g)

# Add Kruskal-Wallis test p-value to the plot
a + annotate("text", x = 1.5, y = max(pi_g$log_pi), label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_result$p.value, digits = 3)))
#This makes it much clearer how nucleotide diversity differs among the populations.

#########################################################################################################################
#####Visualizing patterns along the chromosome ####

### Pi

# Select data of interest
hs_pi <- At_data %>%
  select(mid, ITA_pi, GER_pi, SWE_pi)

# Use gather to rearrange everything
hs_pi_g <- gather(hs_pi, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_pi_g$stat <- factor(hs_pi_g$stat, levels = c("ITA_pi", "GER_pi", "SWE_pi"))

# Take the logarithm of the value variable
hs_pi_g$log_value <- log10(hs_pi_g$value)

# Construct a plot with facets
a_pi <- ggplot(hs_pi_g, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_pi

### FST

# Select data of interest
hs_fst <- At_data %>%
  select(mid, ITA_GER_fst, ITA_SWE_fst, GER_SWE_fst)

# Mutate to set FST and dXY values smaller than zero to zero
hs_fst <- hs_fst %>%
  mutate(across(c(ITA_GER_fst, ITA_SWE_fst, GER_SWE_fst), 
                ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_fst_g <- gather(hs_fst, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_fst_g$stat <- factor(hs_fst_g$stat, levels = c("ITA_GER_fst", "ITA_SWE_fst", "GER_SWE_fst"))

# Take the logarithm of the value variable
hs_fst_g$log_value <- log10(hs_fst_g$value)

# Construct a plot with facets
a_fst <- ggplot(hs_fst_g, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_fst

### Dxy

# Select data of interest
hs_dxy <- At_data %>%
  select(mid, ITA_GER_dxy, ITA_SWE_dxy, GER_SWE_dxy)

# Mutate to set dXY values smaller than zero to zero
hs_dxy <- hs_dxy %>%
  mutate(across(c(ITA_GER_dxy, ITA_SWE_dxy, GER_SWE_dxy), 
                ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_dxy_g <- gather(hs_dxy, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_dxy_g$stat <- factor(hs_dxy_g$stat, levels = c("ITA_GER_dxy", "ITA_SWE_dxy", "GER_SWE_dxy"))

# Take the logarithm of the value variable
hs_dxy_g$log_value <- log10(hs_dxy_g$value)

# Construct a plot with facets
a_dxy <- ggplot(hs_dxy_g, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_dxy

### calculate neutrality statistics
At_sw <- neutrality.stats(At_sw)

get.neutrality(At_sw)

#Let's look at the first population [[1]].
#get.neutrality(At_sw)[[1]]

#Let's look at the first population [[2]].
#get.neutrality(At_sw)[[2]]

#extract Tajma's D
td <- At_sw@Tajima.D/100

# set population names
colnames(td) <- paste0(pops, "_td")
###Delimitate windows on chromosome

#as_tibble: Coerce lists and matrices to data frames
ara_data <- as.tibble(data.frame(windows, td,nd))
nrow(windows)
nrow(nd)
(chr3end-chr3start)/50
nrow(ara_data)
head(ara_data)
ara_data %>% select(contains("pi")) %>% summarise_all(mean)

### load selected positions from chromosome
bed<-read.table("./PopGenome/group_2_At_UV_stress_only.bed")

colnames(bed) <- c("chr", "begin","end")
DF<-vector(length=nrow(ara_data))
ara_data <- as.tibble(data.frame(windows, nd, DF))###if you only want to look at pi
ara_data <- as.tibble(data.frame(windows, nd, td, DF))##if you want to look at tajima D and nucleotide diversity

for (i in 2:nrow(bed)){ara_data$DF[which(ara_data$start>bed$begin[i]&ara_data$stop<bed$end[i]) ]<-"DF"}##each window that overlaps a DF is tagged
ara_data$DF<-as.factor(ara_data$DF)
#summary(ara_data)

# To check and print values in DF other than FALSE, you can use the following approach:
# Check the frequency of each unique value in DF
df_counts <- table(ara_data$DF)
# Print values other than FALSE
#print(df_counts[df_counts > 0 & names(df_counts) != "FALSE"])

# To see the other value(s) in the DF column besides FALSE, you can use the unique() function in 
# combination with a logical condition. Here's how you can do it:

# Extract unique values from the DF column excluding FALSE
other_values <- unique(ara_data$DF[ara_data$DF != "FALSE"])
# Print the other value(s)
#print(other_values)


#### Compare Inter population distributions of Pi and Tajima' D ----

# Compare Pi density distributions of populations 
# modify code change xlim to see complete distribution
plot(density(log(ara_data$ITA_pi)), col= "red", main="Distribution log Pi", xlim=c(-9, max(log(ara_data$GER_pi))))

# Add density plots for other populations with different colors
lines(density(log(ara_data$GER_pi)), col="green")
lines(density(log(ara_data$SWE_pi)), col="blue")

# Add legend
legend("topright", legend=c("ITA", "GER","SWE"),
       col=c("red", "green", "blue"), lty=1,
       title="Population")

# Perform Kruskal-Wallis test
pi_data <- list(
  GER = log(ara_data$GER_pi),
  SWE = log(ara_data$SWE_pi),
  ITA = log(ara_data$ITA_pi)
)

kruskal_result <- kruskal.test(pi_data)

# Perform post-hoc Dunn's test if Kruskal-Wallis is significant
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunn.test(pi_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}

# The post-hoc Dunn's test (dunn.test) will only be performed and printed if the p-value from the 
# Kruskal-Wallis test (kruskal_result$p.value) is less than 0.05.
# If kruskal_result$p.value is greater than or equal to 0.05, the if condition is not satisfied,
# and the code inside the if block, including the dunn.test and print(dunn_result), will not be executed.

# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_result$p.value, digits=4)), 
       bty="n", 
       cex=1)

### Tajima D

#let's say we want to look at mean Tajima'D, we can do that like so:                           
# Select columns containing "td" in the column name and gather them
td_g <- ara_data %>% 
  select(contains("td")) %>% 
  gather(key = "populations", value = "td")

# Remove rows with missing values
td_g2 <- na.omit(td_g)

td_g2$log_td <- log10(td_g2$td)

# Make a boxplot using ggplot2
a_td2 <- ggplot(td_g2, aes(x = populations, y = log_td, fill = populations)) +
  geom_boxplot() +
  theme_light() +
  xlab(NULL) +
  ylab("Tajima's D") +
  ggtitle("Distribution Tajima D")

# Display the plot
print(a_td2)


# Kruskal-Wallis test
kruskal_test_result <- kruskal.test(log_td ~ populations, data = td_g2)

# Add Kruskal-Wallis test p-value to the plot
a_td2 + annotate("text", x = 1.5, y = -1.5, label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_result$p.value, digits = 3)))


### Compare "Tajima's D" density distributions of populations

subset_ITA_td <- log(ara_data$ITA_td[!is.na(ara_data$ITA_td) & !is.nan(log(ara_data$ITA_td))])
subset_GER_td <- log(ara_data$GER_td[!is.na(ara_data$GER_td) & !is.nan(log(ara_data$GER_td))])
subset_SWE_td <- log(ara_data$SWE_td[!is.na(ara_data$SWE_td) & !is.nan(log(ara_data$SWE_td))])

# Compare Pi density distributions of populations from different Altitudes across Spain Spain
# modify code change xlim to see complete distribution
plot(density(subset_ITA_td), col = "red", main = "Distribution log TajimaD", xlim = c(-10,-2))


# Add density plots for other populations with different colors
lines(density(subset_GER_td), col="green")
lines(density(subset_SWE_td), col="blue")


# Add legend
legend("topright", legend=c("ITA", "GER","SWE"),
       col=c("red", "green", "blue"), lty=1,
       title="Population")

# Perform Kruskal-Wallis test
td_data <- list(
  GER = subset_GER_td,
  SWE = subset_SWE_td,
  ITA = subset_ITA_td
)

kruskal_result <- kruskal.test(td_data)

# Perform post-hoc Dunn's test if Kruskal-Wallis is significant
if (kruskal_result$p.value < 0.05) {
  dunn_result <- dunn.test(pi_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}

# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_result$p.value, digits=4)), 
       bty="n", 
       cex=1)
