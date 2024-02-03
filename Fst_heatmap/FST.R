### Load packages----
library(vcfR) # reading VCF
library(StAMPP) # fst
library(ComplexHeatmap)
library("circlize")
library("RColorBrewer")

### FST between populations and individuals----

setwd("~/Downloads/Group2_data/")
vcf_file <- read.vcfR("group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf")
# Convert VCF to genlight object
genlight_vcf <- vcfR2genlight(vcf_file)

# Read population information from file
df <- read.csv("accession_data.csv")
pop <- df$country
subpop <- df$subpopulation

# Extract population data
pop2 = as.factor(pop)
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

# Make an FST heatmap (by countries)
Heatmap(stamppFst_matrix, name = "Fst", column_title = "A", column_title_gp = gpar(fontsize = 50), column_title_side = "top",
        col = colorRamp2(c(0, 0.02, 0.04, 0.06, 0.08), brewer.pal(n=5, name="Reds")))

# Extract subpopulation data
subpop2 = as.factor(subpop)
genlight_vcf$pop = subpop2

# Convert genlight to stampp object
stampp_vcf = stamppConvert(genlight_vcf, type = "genlight")

# Calculate FST between subpopulations
stamppFst = stamppFst(stampp_vcf, nboots = 100, percent = 95, nclusters = 8)
stamppFst_matrix = as.matrix(stamppFst$Fsts)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppFst_matrix) <- 0
stamppFst_matrix[upper.tri(stamppFst_matrix)]  <- t(stamppFst_matrix)[upper.tri(stamppFst_matrix)]
# Optional: order the names
stamppFst_matrix = stamppFst_matrix[order(row.names(stamppFst_matrix)), order(colnames(stamppFst_matrix))]

# Make an FST heatmap (by regions)
Heatmap(stamppFst_matrix, name = "Fst", column_title = "B", column_title_gp = gpar(fontsize = 50), column_title_side = "top",
        col = colorRamp2(c(-0.05, 0, 0.05, 0.10, 0.15, 0.2), brewer.pal(n=6, name="Reds")))

writeLines(capture.output(sessionInfo()), "sessionInfo_fst_heatmap.txt")
