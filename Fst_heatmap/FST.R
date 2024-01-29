### Load packages----
library(vcfR) # reading VCF
library(StAMPP) # fst

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
heatmap(stamppFst_matrix,
        symm = TRUE,
        margins = c(10, 10))
title(main = "A", adj=0.15, cex.main = 3)

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
heatmap(stamppFst_matrix,
        symm = TRUE,
        margins = c(10, 10),
        cexRow = 1.5,
        cexCol = 1.5)
title(main = "B", adj=0.15, cex.main = 3)

#writeLines(capture.output(sessionInfo()), "sessionInfo_fst_heatmap.txt")
