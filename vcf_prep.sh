###### Rscipt to Annotate GFF file ####
# Input files below are downloaded from: https://www.arabidopsis.org/
# TAIR10_functional_descriptions.txt
# TAIR10_GFF3_genes.gff
#######################################
#set working directory
setwd("C:/Data_tali/QBio/2023-24/AT_1001_genome/Annotate_GFF")

###
#Read GFF file with SNP positions
###

GFF<-read.delim("TAIR10_GFF3_genes.gff", header = FALSE)
head(GFF)
#Create a duplicated column 9
GFF$dup_colmn<-GFF$V9

#select all the rows with mRNA in cloumn 3
GFF_sub<-subset(GFF, GFF$V3 == "mRNA")
head(GFF_sub)
GFF_sub$dup_colmn<-gsub("ID=.*Name=","", GFF_sub$dup_colmn)
GFF_sub$dup_colmn<-gsub(";Index=1","", GFF_sub$dup_colmn)

write.table(GFF_sub, "GFF_final_dup_column.gff", sep = "\t", quote = FALSE)

#change name of column from duplicated to Model_name
GFF_sub$Model_name<-GFF_sub$dup_colmn

###
#Read annotation file
###
annot<-read.delim("TAIR10_functional_descriptions.txt", header = TRUE)
head(annot)
#Merge subset of GFF with annotation file using column Model_name (common in both data files)
GFF_sub_annot.m<-merge(GFF_sub, annot, by="Model_name")
head(GFF_sub_annot.m)
write.table(GFF_sub_annot.m, "GFF_final_with_annotation.gff", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

###
# Subset VCF file using annotated gff file | This part is implemented in the Linux terminal
# VCF file is downloaded using this command: $ wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
###

# apply quality filters on variants in the VCF file
$ vcftools --gzvcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz --minDP 10 --minGQ 20 --minQ 30 --max-missing 0.80 --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > /scratch/QBIO/VCF/1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz

# Keep only selected accesions from the VCF file
$ vcftools --gzvcf /scratch/QBIO/VCF/1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz --keep /scratch/QBIO/student_project_data/group_2_accession_names.txt --recode --recode-INFO-all --stdout | bgzip -c > /scratch/QBIO/student_project_data/group_2_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz

# grep required genes from annotated gff file with coordinates
$ grep -e "Glutathione S-transferase" -e "Ascorbate peroxidase" -e "photolyase 1" -e "RAD23B" -e "Y-family DNA polymerase H" -e "RAD4" -e "Rad23 UV excision repair protein" -e "single-stranded DNA endonuclease family protein" -e "5'-3' exonuclease family protein" GFF_final_with_annotation_2_gff > group_2_At_UV_stress_only.gff

# select columns from filtered gff and save as bed file
$ awk '{OFS="\t"; print $2, $5, $6}' group_2_At_UV_stress_only.gff > group_2_At_UV_stress_only.bed

# Filter VCF files to subset variants based on desired genomic regions (genes)
# Make sure the formate for CHROM column in the bed file and VCF file is the same
$ vcftools --gzvcf group_2_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz --bed group_2_At_UV_stress_only.bed --stdout --recode --keep-INFO-all | bgzip -c > group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz

# index VCF file using tabix
$ tabix -p vcf group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz
$ zgrep -v "^#" group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz | wc -l

### 