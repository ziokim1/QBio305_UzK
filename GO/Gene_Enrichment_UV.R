# Load necessary libraries
library(BiocManager)
library(KEGGREST)
library(org.At.tair.db)
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)
library(scales)

# Set the working directory. Make sure to adjust this path to your specific location.
setwd("~/Downloads/Group2_data/GO/")

# Read the input file (gene expression values and significance)
#gene_list <- read.csv("output.csv") 

#gene_list <- read.table("GO_input_ITA_GER.txt", header=TRUE) # output.csv
#gene_list <- read.table("GO_input_ITA_SWE.txt", header=TRUE) # output.csv
gene_list <- read.table("GO_input_GER_SWE.txt", header=TRUE) # output.csv


# Create a universe (list of all genes) for enrichments based on significance
univ <- gene_list[, 3]         # Select the third column with p-values
names(univ) <- gene_list[, 1]  # Give names to genes from the first column
univ <- univ[!is.na(univ)]     # Ignore NA values

# Filter out genes not in the database
#univ_clean <- gsub("\\..*", "", names(univ))
#univ_clean <- univ_clean[univ_clean %in% keys(org.At.tair.db)]

# Define a function to return TRUE/FALSE for p-values < 0.05
selection <- function(allScore) { return(allScore < 0.05) }

# Prepare topGO data
tGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", geneSel = selection, allGenes = univ, nodeSize = 3, mapping = "org.At.tair.db", annot = annFUN.org)

# Run Gene Enrichment Analysis
results.ks <- runTest(tGOdata, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichment <- GenTable(tGOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS < 0.05, ]
goEnrichment <- goEnrichment[, c("GO.ID", "Term", "KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")

class(goEnrichment$KS)
# If the class is not numeric, you'll need to convert it to numeric:
goEnrichment$KS <- as.numeric(goEnrichment$KS)

# check for NAs
any(is.na(goEnrichment$KS))
head(goEnrichment[is.na(goEnrichment$KS), ])

# remove NAs
goEnrichment <- goEnrichment[!is.na(goEnrichment$KS), ]

# Save the enriched GO terms table to a CSV file
#write.csv(goEnrichment, "sample_project_ITA_GER.csv")
#write.csv(goEnrichment, "sample_project_ITA_SWE.csv")
write.csv(goEnrichment, "sample_project_GER_SWE.csv")

# Read the exported file
#UV <- read.csv("export_file.csv")

# Select the top N terms
#ntop <- 10
#ggdata <- UV[1:ntop, ]

# Select the top N terms
ntop <- 10
ggdata <- goEnrichment[1:ntop, ]

ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))  # Fix order

# Plotting using ggplot2
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5, 12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'Pathways Enriched in Upregulated Genes in Stress (GER-SWE)',
    subtitle = 'Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
  ) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),colour = c("black", "black", "black"),size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12), axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'), axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),legend.key.size = unit(1, "cm"),legend.text = element_text(size = 16),
    title = element_text(size = 12)) + coord_flip()

# Display and save the plot
print(gg1)
#ggsave("GO_results_UV_ITA_GER.png", plot = gg1, width = 11, height = 10, dpi = 300)
#ggsave("GO_results_UV_ITA_SWE.png", plot = gg1, width = 11, height = 10, dpi = 300)
ggsave("GO_results_UV_GER_SWE.png", plot = gg1, width = 11, height = 10, dpi = 300)
