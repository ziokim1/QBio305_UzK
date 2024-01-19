############################
# Prepare input file for   #
# Gene Enrichment Analysis #
############################

####################
# Create input file "1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing_80SwEs_DefenseOnly.p.fst_SWE-ESP.tsv"
# by running bash script "stacks.sh" in Cheops terminal
###################

######
# login to cheops and create this bash script file using following code
# and run the as "sbatch stacks.sh"
######
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --partition=devel-rh7
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln
#SBATCH --error=/home/group.kurse/qbcbp009/Group2/-%x-%j.err
#SBATCH --output=/home/group.kurse/qbcbp009/Group2/-%x-%j.out

module load stacks/2.65

# Input directory # don'T forget to change this directory where you have saved 
# "input.vcf" and "pop.txt" files

input_dir="/home/group.kurse/qbcbp009/Group2/"

# Output directory, don't forget to change this to your stacks output directory
output_dir="/home/group.kurse/qbcbp009/Group2/stacks"

populations -V "${input_dir}group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz" --popmap "${input_dir}pop.txt" --fstats -t 10 --structure -O "${output_dir}"

############
# Rscript to process stacks output file "*.p.fst_ITA_SWE.tsv", "*.p.fst_ITA_GER.tsv" and "*.p.fst_GER_SWE.tsv"
############
# Load required libraries
library(readr)
library(GenomicRanges)
library(dplyr)

# set working directory
setwd("~/Downloads/Group2_data/GO/")

######
# read --fstats output file from stacks
######
#fstats <- read.table("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_ITA-GER.tsv", header=FALSE)
#fstats <- read.table("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_ITA-SWE.tsv", header=FALSE)
fstats <- read.table("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_GER-SWE.tsv", header=FALSE)

# Set the column names manually and replace spaces with underscores
colnames(fstats) <- gsub(" ", "_", c(
  "Locus ID", "Pop 1 ID", "Pop 2 ID", "Chr", "BP", "Column", "Fishers P",
  "Odds Ratio", "CI Low", "CI High", "LOD", "AMOVA Fst", "Smoothed AMOVA Fst",
  "Smoothed AMOVA Fst P-value"
))

names(fstats)
# Create a subset of the dataframe with selected columns
fstats_fst <- data.frame(
  Chr = fstats$Chr,
  BP = fstats$BP,
  Fishers_P = fstats$Fishers_P,
  AMOVA_Fst = fstats$AMOVA_Fst
)

names(fstats_fst)

head(fstats_fst)
print(head(fstats_fst, n=1))

#############
#Read annotation file
#############
annot<-read.delim("group_2_At_UV_stress_only.gff", header = FALSE)
head(annot)
names(annot)

print(head(annot, n=1))

# Set the column names manually and replace spaces with underscores
colnames(annot) <- c(
  "gene_id", "chr", "tair_version", "type", "start", "end", "empty1",
  "strand", "empty2", "gene_id2", "gene_id3", "type_name", "short_description", "curation"
)

print(head(annot, n=5))

# Create a subset of the dataframe with selected columns
annot_subset <- data.frame(
  gene_id = annot$gene_id,
  chr = annot$chr,
  start = annot$start,
  end = annot$end
)  

print(head(annot_subset, n=5))

#remove "Chr" from each entry in the annot_subset$chr column to make it similar as in fstats_fst$Chr:
annot_subset$chr <- sub("Chr", "", annot_subset$chr)


# Convert data frames to GRanges objects
gr_fst <- GRanges(seqnames = fstats_fst$Chr, 
                  ranges = IRanges(start = fstats_fst$BP, end = fstats_fst$BP),
                  fst = fstats_fst$AMOVA_Fst,
                  p_value = fstats_fst$Fishers_P)

length(fstats_fst$Chr)
print(head(gr_fst, n=1))

gr_annotation <- GRanges(seqnames = annot_subset$chr, 
                         ranges = IRanges(start = annot_subset$start, end = annot_subset$end),
                         gene_id = annot_subset$gene_id)

# Find overlaps between fst_data and annotation_data
ovl <- findOverlaps(gr_fst, gr_annotation, type = "any", select = "all", ignore.strand = TRUE)

# Extract Gene IDs based on overlaps
overlapping_genes <- as.data.frame(gr_annotation[subjectHits(ovl)])
overlapping_fst <- as.data.frame(gr_fst[queryHits(ovl)])

# Print or further process the results
print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

# To change the name of the seqnames column to chrom and remove the strand column in both
# dataframes (overlapping_genes and overlapping_fst)

# Rename 'seqnames' to 'chrom'
colnames(overlapping_genes)[colnames(overlapping_genes) == "seqnames"] <- "chrom"
colnames(overlapping_fst)[colnames(overlapping_fst) == "seqnames"] <- "chrom"

# Remove the 'strand' column
overlapping_genes <- overlapping_genes[, !(names(overlapping_genes) %in% "strand")]
overlapping_fst <- overlapping_fst[, !(names(overlapping_fst) %in% "strand")]

# Print or further process the results
print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

# Write both overlapping_genes and overlapping_fst as text files
#write.table(overlapping_genes, file = "overlapping_genes_ITA_GER.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(overlapping_fst, file = "overlapping_fst_ITA_GER.txt", sep = "\t", row.names = FALSE, quote = FALSE)  

#write.table(overlapping_genes, file = "overlapping_genes_ITA_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(overlapping_fst, file = "overlapping_fst_ITA_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)  

write.table(overlapping_genes, file = "overlapping_genes_GER_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(overlapping_fst, file = "overlapping_fst_GER_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)  


####################
# Run this part in terminal
###################

# Copy these file to Cheops and run bedtools intersect to get loci which are within
# genomic boundaries of defense related genes
$ module load bedtools/
$ bedtools intersect -wa -wb -a overlapping_genes_ITA_GER.txt -b overlapping_fst_ITA_GER.txt > overlapping_genes_fst_ITA_GER.txt

# bedtools intersect -wa -wb -a overlapping_genes_ITA_SWE.txt -b overlapping_fst_ITA_SWE.txt > overlapping_genes_fst_ITA_SWE.txt
# bedtools intersect -wa -wb -a overlapping_genes_GER_SWE.txt -b overlapping_fst_GER_SWE.txt > overlapping_genes_fst_GER_SWE.txt

##################
# Process "overlapping_genes_fst.txt"
# in R
#################

#Read annotation file

#ovl_genes_fst<-read.delim("overlapping_genes_fst_ITA_GER.txt", header = FALSE)
#ovl_genes_fst<-read.delim("overlapping_genes_fst_ITA_SWE.txt", header = FALSE)
ovl_genes_fst<-read.delim("overlapping_genes_fst_GER_SWE.txt", header = FALSE)

head(ovl_genes_fst)
names(ovl_genes_fst)

print(head(ovl_genes_fst, n=1))

# Set the column names
colnames(ovl_genes_fst) <- c(
  "chrom", "gene_start", "gene_end", "gene_width", "gene_id","chrom", "snp_start", "snp_end", "snp_width",
  "fst", "p_value")

print(head(ovl_genes_fst, n=5))

# subset ovl_genes_fst
# Create a subset of the "ovl_genes_fst" dataframe with selected columns
ovl_genes_fst_sub <- data.frame(
  gene_id = ovl_genes_fst$gene_id,
  fst = ovl_genes_fst$fst,
  p_value = ovl_genes_fst$p_value
)  

print(head(ovl_genes_fst_sub, n=5))

#Drop Integer After Decimal Point in gene_id:
ovl_genes_fst_sub$gene_id <- sub("\\.\\d+", "", ovl_genes_fst_sub$gene_id)

# Order by fst (in descending order) and then by p_value (in ascending order) within each gene_id:

sorted_data <- ovl_genes_fst_sub %>%
  arrange(gene_id, desc(fst), p_value)

# Keep Unique gene_id with Highest fst but Lowest p_value:
unique_data <- sorted_data %>%
  group_by(gene_id) %>%
  filter(row_number() == 1)

print(head(unique_data, n=10))

########################

# Write both overlapping_genes and overlapping_fst as text files
#write.table(unique_data, file = "GO_input_ITA_GER.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(unique_data, file = "GO_input_ITA_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(unique_data, file = "GO_input_GER_SWE.txt", sep = "\t", row.names = FALSE, quote = FALSE)
